
extern crate log;

//use std::cmp::{Ordering, Reverse, max, min};
//use std::cmp::Reverse;
//use std::collections::{BinaryHeap, HashSet};
use std::collections::BinaryHeap;
use std::sync::Arc;
use triple_accel::{levenshtein_search,Match};

use crate::bv_bwt::BitVectorBWT;
use crate::dynamic_wfa::DynamicWFA;
use crate::stats_util;
use crate::string_util;

const VALID_CHARS: [u8; 4] = [1, 2, 3, 5];
const VALID_CHARS_LEN: usize = VALID_CHARS.len();

/// stores options for running the correction algorithms
pub struct CorrectionParameters {
    //use_fm_index: bool, //old parameter that toggled between bit vector and classic mode
    /// The k-mer sizes to use for correction, performs one pass per value in the array
    pub kmer_sizes: Vec<usize>,
    /// The absolute minimum k-mer count to be considered present
    pub min_count: u64,
    /// The maximum length of sequence to attempt to correct
    pub max_branch_attempt_length: usize,
    /// A multiplier on `k` to limit the number of branches explored
    pub branch_limit_factor: f64,
    /// A buffering factor to allow bridge corrections to be longer than the original sequence
    pub branch_buffer_factor: f64,
    //TODO: add a tail_truncate_factor that buts a bounding box around min length and max length
    /// A factor limiting partial corrections when bridges cannot be found
    pub midpoint_ed_factor: f64,
    /// A buffering factor to allow assembly corrections to be longer than the original sequence
    pub tail_buffer_factor: f64,
    /// A multiplier factor on the median for dynamic minimum k-mer count to be considered present
    pub frac: f64,
    //fm_bit_power: u8, //only matters for classic mode which isn't implemented currently
    /// Will calculate more stats if verbose is set to `true`
    pub verbose: bool
}

/// a basic struct for storing a correction to a sequence
pub struct Correction {
    start_pos: usize,
    end_pos: usize,
    seq: Vec<u8>
}

/// a struct for storing generic read
#[derive(Clone,Debug)]
pub struct LongReadFA {
    /// The index associated with the read
    pub read_index: u64,
    /// The read label/identifier
    pub label: String,
    /// The actual genomic sequence
    pub seq: String
}

/// a struct for storing the modified string
#[derive(Clone,Debug)]
pub struct CorrectionResults {
    /// The index associated with the read
    pub read_index: u64,
    /// The read label/identifier
    pub label: String,
    /// The original, uncorrected sequence
    pub original_seq: String,
    /// The modified, corrected sequence
    pub corrected_seq: String,
    /// If verbose is set, this will store the average k-mer count before correction
    pub avg_before: f64,
    /// If verbose is set, this will store the average k-mer count after correction
    pub avg_after: f64
}

/// This will run a correction "job" on a single long read using a shared BWT resource.
/// # Argument
/// * `arc_bwt` - the shared BitVectorBWT resource
/// * `long_read` - the read to correct
/// * `arc_params` - the shared parameters to use for performing the correction
pub fn correction_job(arc_bwt: Arc<BitVectorBWT>, long_read: LongReadFA, arc_params: Arc<CorrectionParameters>) -> CorrectionResults {
    //these clearly are not mutable, nor should they be
    let bwt: &BitVectorBWT = &*arc_bwt;
    let params: &CorrectionParameters = &*arc_params;

    //convert the input
    let mut seq_i: Vec<u8> = string_util::convert_stoi(&long_read.seq);

    //now do a correction pass for each k-mer in the correction params
    for k in params.kmer_sizes.iter() {
        seq_i = correction_pass(bwt, &seq_i, params, *k);
    }
    
    //convert back to a string
    let corrected_seq: String = string_util::convert_itos(&seq_i);
    
    let avg_before: f64 = if params.verbose {
        let orig_seq: Vec<u8> = string_util::convert_stoi(&long_read.seq);
        let counts_before: Vec<u64> = bwt.count_pileup(&orig_seq, params.kmer_sizes[0]);
        let sum_before: u64 = counts_before.iter().sum();
        sum_before as f64 / counts_before.len() as f64
    } else {
        0.0
    };
    let avg_after: f64 = if params.verbose {
        //TODO: is there a way to pre-serve this from earlier? this mode will always take longer, so it's not a huge deal
        let counts_after: Vec<u64> = bwt.count_pileup(&seq_i, params.kmer_sizes[0]);
        let sum_after: u64 = counts_after.iter().sum();
        sum_after as f64 / counts_after.len() as f64
    } else {
        0.0
    };

    //send it on back y'all
    CorrectionResults {
        read_index: long_read.read_index,
        label: long_read.label,
        original_seq: long_read.seq,
        corrected_seq,
        avg_before,
        avg_after
    }
}

/// This is the core correction function. It takes a BWT, a sequence, parameters, and a k-mer size and performs the
/// following high level steps: calculate expected k-mer counts, identify gaps, bridge/extend from solid k-mers, and finally
/// build the corrected sequence.
/// # Arguments
/// * `bwt` - the BitVectorBWT that contains all the implicit k-mer counts
/// * `seq_i` - the sequence to correct in integer form
/// * `params` - the correction parameters to use
/// * `kmer_size` - the length of k-mer (i.e. `k`) to use for this pass
pub fn correction_pass(bwt: &BitVectorBWT, seq_i: &[u8], params: &CorrectionParameters, kmer_size: usize) -> Vec<u8> {
    //count up the initial pileups
    let pileup: Vec<u64> = bwt.count_pileup(seq_i, kmer_size);
    let nz_med = stats_util::calculate_bounded_median(&pileup, params.min_count);
    if nz_med < params.min_count {
        return seq_i.to_owned();
    }
    let pileup_size: usize = pileup.len();
    
    //we can try some fixing, set up additional params for this pass
    let branch_limit = params.branch_limit_factor as usize*kmer_size;
    
    //try to dynamically set the threshold, but make sure its at least MIN_COUNT
    let threshold: u64 = std::cmp::max((params.frac * nz_med as f64) as u64, params.min_count);
    let upper_threshold: u64 = 10*nz_med;

    //prep for the actual corrections now
    let mut prev_found: isize = -1;
    let mut seed_kmer: Vec<u8> = vec![0; kmer_size];
    let mut target_kmer: Vec<u8> = vec![0; kmer_size];
    let mut max_branch_length: usize;

    let mut bridge_points: Vec<Vec<u8>>;
    let mut corrections_list: Vec<Correction> = Vec::<Correction>::new();
    let mut new_corr: Correction;
    
    let mut x: usize = 0;

    while x < pileup_size {
        if pileup[x] < threshold {
            prev_found = x as isize-1;
            
            //find the next index that is above the threshold
            while x < pileup_size && pileup[x] < threshold {
                x += 1;
            }

            if prev_found == -1 && x < pileup_size {
                //handle the head case
                max_branch_length = (params.tail_buffer_factor*(x+kmer_size) as f64) as usize;
                if max_branch_length < params.max_branch_attempt_length {
                    //copy and reverse complement since we're going backwards
                    seed_kmer[..].clone_from_slice(&seq_i[x..x+kmer_size]);
                    seed_kmer = string_util::reverse_complement_i(&seed_kmer);
                    
                    //now assemble out from it
                    bridge_points = assemble_from_kmer(bwt, &seed_kmer, threshold, branch_limit, max_branch_length);

                    //remember to rev comp this also
                    let orig: Vec<u8> = string_util::reverse_complement_i(&seq_i[0..x+kmer_size]);
                    //this is the dynamic form, but this is worse in our test?
                    //let max_edit_distance: u32 = ((x+kmer_size) as f64 * params.midpoint_ed_factor) as u32;
                    let max_edit_distance = 0xFFFFFFFF;
                    let best_match: Option<Vec<u8>> = pick_best_levenshtein_search(&orig, bridge_points, bwt, kmer_size, max_edit_distance);
                    if let Some(bm) = best_match {
                        //we found a best match, it will span from index 0 up to the start of the first found k-mer
                        let rev_comp_seq: Vec<u8> = string_util::reverse_complement_i(&bm);
                        
                        //now store the correction in range [0..x)
                        new_corr = Correction {
                            start_pos: 0,
                            end_pos: x+kmer_size,
                            seq: rev_comp_seq
                        };
                        corrections_list.push(new_corr);
                    }
                }
            }
            else if prev_found >= 0 && x < pileup_size {
                //handle a bridging case
                seed_kmer[..].clone_from_slice(&seq_i[prev_found as usize..prev_found as usize+kmer_size]);
                target_kmer[..].clone_from_slice(&seq_i[x..x+kmer_size]);
                max_branch_length = (params.branch_buffer_factor*(x-prev_found as usize+kmer_size) as f64) as usize;

                if max_branch_length < params.max_branch_attempt_length {
                    //try forward first
                    //bridge_points = bridge_kmers(bwt, &seed_kmer, &target_kmer, threshold, branch_limit, max_branch_length);
                    bridge_points = bridge_sequence(
                        bwt, 
                        &seq_i[prev_found as usize..x+kmer_size], 
                        kmer_size, 
                        threshold, 
                        upper_threshold,
                        branch_limit, 
                        max_branch_length
                    );
                    
                    //try reverse complement if we failed
                    if bridge_points.is_empty() {
                        bridge_points = bridge_sequence(
                            bwt, 
                            &string_util::reverse_complement_i(&seq_i[prev_found as usize..x+kmer_size]),
                            kmer_size,
                            threshold,
                            upper_threshold,
                            branch_limit,
                            max_branch_length
                        );
                        
                        //make sure to rev-comp the results here
                        for bp in &mut bridge_points {
                            *bp = string_util::reverse_complement_i(&bp);
                        }
                    }
                    
                    //we are doing a bridge, so use straight levenshtein
                    //let best_match: Option<Vec<u8>> = pick_best_levenshtein(&seq_i[prev_found as usize..x+kmer_size], bridge_points, bwt, kmer_size);
                    let best_match: Option<Vec<u8>> = pick_best_pileup(bridge_points, bwt, kmer_size);

                    //pick out the best result from the full levenshtein
                    if let Some(bm) = best_match {
                        //now store the correction in range [0..x)
                        new_corr = Correction {
                            start_pos: prev_found as usize,
                            end_pos: x+kmer_size,
                            seq: bm
                        };
                        corrections_list.push(new_corr);
                    } else if x - prev_found as usize > kmer_size {
                        //no bridges were found, try to extend into the midpoint (assuming no seed/target overlap)
                        let mid_point: usize = (prev_found as usize+x+kmer_size) / 2;
                        max_branch_length = (params.tail_buffer_factor*(mid_point - prev_found as usize) as f64) as usize;
                        let max_edit_distance: u32 = ((mid_point - prev_found as usize) as f64 * params.midpoint_ed_factor) as u32;

                        //extend left to midpoint
                        bridge_points = assemble_from_kmer(bwt, &seed_kmer, threshold, branch_limit, max_branch_length);
                        
                        //this and below did not work like I originally thought; it may need a slightly different
                        //approach that attempt to find the optimum edit distance while extending
                        /*
                        bridge_points = bridge_sequence(
                            bwt, 
                            &seq_i[prev_found as usize..mid_point], 
                            kmer_size, 
                            threshold, 
                            upper_threshold,
                            branch_limit, 
                            max_branch_length,
                            true
                        );
                        */

                        //let best_match: Option<Vec<u8>> = pick_best_pileup(bridge_points, bwt, kmer_size);
                        let best_match: Option<Vec<u8>> = pick_best_levenshtein_search(&seq_i[prev_found as usize..mid_point], bridge_points, bwt, kmer_size, max_edit_distance);
                        if let Some(bm) = best_match {
                            new_corr = Correction {
                                start_pos: prev_found as usize,
                                end_pos: mid_point,
                                seq: bm
                            };
                            corrections_list.push(new_corr);
                        }

                        //extend right to midpoint
                        let rev_target = string_util::reverse_complement_i(&target_kmer);
                        bridge_points = assemble_from_kmer(bwt, &rev_target, threshold, branch_limit, max_branch_length);
                        /*
                        bridge_points = bridge_sequence(
                            bwt, 
                            &string_util::reverse_complement_i(&seq_i[mid_point..x+kmer_size]),
                            kmer_size, 
                            threshold, 
                            upper_threshold,
                            branch_limit, 
                            max_branch_length,
                            true
                        );
                        */
                        
                        //remember to rev comp this also
                        let orig: Vec<u8> = string_util::reverse_complement_i(&seq_i[mid_point..x+kmer_size]);
                        let best_match: Option<Vec<u8>> = pick_best_levenshtein_search(&orig, bridge_points, bwt, kmer_size, max_edit_distance);
                        //let best_match: Option<Vec<u8>> = pick_best_pileup(bridge_points, bwt, kmer_size);
                        if let Some(bm) = best_match {
                            //need to rev-comp the result, and it goes from midpoint to the right side
                            let rev_comp_seq: Vec<u8> = string_util::reverse_complement_i(&bm);
                            new_corr = Correction {
                                start_pos: mid_point,
                                end_pos: x+kmer_size,
                                seq: rev_comp_seq
                            };
                            corrections_list.push(new_corr);
                        }
                    }
                }
            }
        }
        else {
            //the counts were okay, no correction needed here
            x += 1;
        }
    }

    //handle any tail sequences
    if prev_found >= 0 && pileup[pileup_size-1] < threshold {
        max_branch_length = (params.tail_buffer_factor*(seq_i.len()-prev_found as usize) as f64) as usize;
        //max_branch_length = (params.branch_buffer_factor*(pileup_size-1-prev_found as usize+kmer_size) as f64) as usize;
        if max_branch_length <= params.max_branch_attempt_length {
            //copy the last found k-mer and assemble outwards
            seed_kmer[..].clone_from_slice(&seq_i[prev_found as usize..prev_found as usize+kmer_size]);
            bridge_points = assemble_from_kmer(bwt, &seed_kmer, threshold, branch_limit, max_branch_length);
            //bridge_points = wave_correction(bwt, &seq_i[prev_found as usize..], kmer_size, threshold, branch_limit, max_branch_length, false);
            /*
            //this method doesn't seem to work as well, and it's a lot slower for some reason
            //probably need to consider how to do this in the one directional assembly approach
            bridge_points = bridge_sequence(
                bwt, 
                &seq_i[prev_found as usize..], 
                kmer_size, 
                threshold, 
                upper_threshold,
                branch_limit, 
                max_branch_length,
                true
            );*/

            //now get the best match
            //this is the dynamic form, but this is worse in our test?
            //let max_edit_distance: u32 = (params.midpoint_ed_factor*(seq_i.len()-prev_found as usize) as f64) as u32;
            let max_edit_distance: u32 = 0xFFFFFFFF;
            let best_match: Option<Vec<u8>> = pick_best_levenshtein_search(&seq_i[prev_found as usize..], bridge_points, bwt, kmer_size, max_edit_distance);
            if let Some(bm) = best_match {
                //now store the correction in range [0..x)
                new_corr = Correction {
                    start_pos: prev_found as usize,
                    end_pos: seq_i.len(),
                    seq: bm
                };
                corrections_list.push(new_corr);
            }
        }
    }

    let mut current_position: usize = 0;
    let mut corrected_seq: Vec<u8> = Vec::<u8>::new();
    for correction in corrections_list {
        if current_position > correction.start_pos {
            corrected_seq.truncate(corrected_seq.len() - (current_position - correction.start_pos));
        }
        else {
            //copy everything up to the correction
            corrected_seq.extend_from_slice(&seq_i[current_position..correction.start_pos]);
        }

        //copy in the correction
        corrected_seq.extend_from_slice(&correction.seq);
        
        //update our current position
        current_position = correction.end_pos as usize;
    }
    //add in anything through the end and send it back
    corrected_seq.extend_from_slice(&seq_i[current_position..]);
    corrected_seq
}

/// This function will take an original sequence and a candidate list and pick the best candidate to return
/// after performing some sanity checks on the mapping. It assumes that all matches must requires the first `kmer_size`
/// characters to exactly match between the original sequence and all candidates.
/// # Arguments
/// * `original` - the original sequence in integer format
/// * `candidates` - the Vec of candidates, each in integer format
/// * `bwt` - the BWT of counts (used for pileup tie-breaking)
/// * `kmer_size` - the k-mer size to use for pileup tie-breaking
/// * `max_ed` - the maximum allowed edit distance, set to 0xFFFFFFFF if anything goes
#[inline]
fn pick_best_levenshtein_search(original: &[u8], candidates: Vec<Vec<u8>>, bwt: &BitVectorBWT, kmer_size: usize, max_ed: u32) -> Option<Vec<u8>> {
    let mut ed_scores: Vec<Option<Match>> = Vec::<Option<Match>>::with_capacity(candidates.len());
    let mut min_score: u32 = max_ed;
    for candidate in candidates.iter() {
        //calculate the min distance
        let matches: Vec<Match> = levenshtein_search(&original, &candidate).collect();
        let mut best_match: Option<Match> = None;
        for m in matches {
            //make sure the start is index 0 and goes at least k long
            if m.start == 0 && m.end >= kmer_size {
                match &best_match {
                    Some(bm) => {
                        if m.k < bm.k || (m.k == bm.k && m.end < bm.end) {
                            best_match = Some(m);
                        }
                    },
                    None => {
                        best_match = Some(m);
                    }
                }
            }
        }

        //update the min score if we can
        match &best_match {
            Some(bm) => {
                if bm.k < min_score {
                    min_score = bm.k;
                }
            },
            None => {}
        }
        ed_scores.push(best_match);
    }
    
    //get everything with a good score
    let mut candidates_ed: Vec<&[u8]> = Vec::<&[u8]>::with_capacity(candidates.len());
    for y in 0..candidates.len() {
        match &ed_scores[y] {
            Some(eds) => {
                if eds.k == min_score {
                    //we have to truncate down to the match end
                    candidates_ed.push(&candidates[y][0..eds.end]);
                }
            },
            None => {}
        }
    }

    if candidates_ed.is_empty() {
        //do nothing, we didn't find anything good
        None
    }
    else if candidates_ed.len() == 1 {
        //only one with smallest edit distance
        Some(candidates_ed[0].to_vec())
    }
    else {
        //figure out which of the ones with equal edit distance has the most counts
        let mut max_counts: u64 = 0;
        let mut max_id: usize = 0;
        let mut ed_pu: Vec<u64>;
        let mut summation: u64;
        for (y, candidate) in candidates_ed.iter().enumerate() {
            ed_pu = bwt.count_pileup(candidate, kmer_size);
            summation = ed_pu.iter().sum();
            if summation > max_counts {
                max_id = y;
                max_counts = summation;
            }
        }

        //now return the best candidate
        Some(candidates_ed[max_id].to_vec())
    }
}

/// This function will take a candidate list and pick the best candidate to return based on the pileup support
/// in the BWT. 
/// # Arguments
/// * `candidates` - the Vec of candidates, each in integer format; for speed, we will remove the best candidate directly from this list
/// * `bwt` - the BWT of counts (used for pileup tie-breaking)
/// * `kmer_size` - the k-mer size to use for pileup tie-breaking
#[inline]
fn pick_best_pileup(mut candidates: Vec<Vec<u8>>, bwt: &BitVectorBWT, kmer_size: usize) -> Option<Vec<u8>> {
    //two short circuit points
    if candidates.is_empty() {
        //no valid bridges were found
        None
    }
    else if candidates.len() == 1 {
        //only one valid bridge, so short circuit here
        Some(candidates.remove(0))
    }
    else {
        //used to check ED here, but no longer necessary with WFA already handling it
        //TODO: is pileup or random better here?
        //figure out which of the ones with equal edit distance has the most counts
        let mut max_counts: u64 = 0;
        let mut max_id: usize = 0;
        let mut ed_pu: Vec<u64>;
        let mut summation: u64;
        for (y, candidate) in candidates.iter().enumerate() {
            ed_pu = bwt.count_pileup(&candidate, kmer_size);
            summation = ed_pu.iter().sum();
            if summation > max_counts {
                max_id = y;
                max_counts = summation;
            }
        }

        //now return the best candidate
        Some(candidates.remove(max_id))
    }
}

/// Traverses from one kmer to another using the BWT for querying counts
/// # Arguments
/// * `bwt` - the BWT that contains count data
/// * `seed_kmer` - the integer form seed k-mer for the bridge
/// * `target_kmer` - the integer form target k-mer for the bridge
/// * `min_count` - the minimum count required for a path to be consider solid
/// * `branch_limit` - the maximum number of branches to explore before giving up on all paths
/// * `max_branch_len` - the maximum branch length allowed before giving up on a path
pub fn bridge_kmers(
    bwt: &BitVectorBWT, seed_kmer: &[u8], target_kmer: &[u8], min_count: u64, branch_limit: usize, 
    max_branch_len: usize
) -> Vec<Vec<u8>> {
    //build up return
    let mut ret = Vec::<Vec<u8>>::new();

    //build some helper values we'll be looping through a lot
    let kmer_len = seed_kmer.len();
    assert_eq!(kmer_len, target_kmer.len());
    let mut counts: [u64; VALID_CHARS_LEN] = [0; VALID_CHARS_LEN];
    let mut fw_counts: [u64; VALID_CHARS_LEN] = [0; VALID_CHARS_LEN];
    let mut rev_counts: [u64; VALID_CHARS_LEN] = [0; VALID_CHARS_LEN];
    let mut num_branched: usize = 0;
    let mut max_pos: usize;

    //these buffers are used to create query slices
    let mut curr_buffer: Vec<u8> = vec![4; max_branch_len];
    let mut rev_buffer: Vec<u8> = vec![4; max_branch_len];

    //initialize the bridging with our seed k-mer
    let mut possible_bridges: Vec<Vec<u8>> = Vec::<Vec<u8>>::new();
    let mut seed_vec: Vec<u8> = vec![0; kmer_len];
    seed_vec.clone_from_slice(seed_kmer);
    possible_bridges.push(seed_vec);

    while !possible_bridges.is_empty() && num_branched < branch_limit {
        //get a bridge to extend
        let mut curr_bridge: Vec<u8> = possible_bridges.pop().unwrap();
        let mut curr_bridge_len = curr_bridge.len();
        let mut curr_offset: usize = 0;
        num_branched += 1;

        //TODO: replace this with copy slice? not sure how to do the rev comp in a one-line
        for x in 0..kmer_len {
            curr_buffer[x] = curr_bridge[curr_bridge_len-kmer_len+x];
            rev_buffer[x] = string_util::COMPLEMENT_INT[curr_buffer[x] as usize];
        }

        while curr_bridge_len < max_branch_len {
            //increase the offset into our buffers (i.e. shift the curr k-mer left and rev k-mer right)
            curr_offset += 1;
            
            //do all the k-mer counting, efficient on rev-comp, then forward queries are added in
            bwt.prefix_revkmer_noalloc_fixed(&rev_buffer[curr_offset..curr_offset+kmer_len-1], &mut rev_counts);
            
            //parallel forward query, doesn't seem to gain much if anything
            bwt.postfix_kmer_noalloc_fixed(&curr_buffer[curr_offset..curr_offset+kmer_len-1], &mut fw_counts);
            
            max_pos=0;
            for x in 0..VALID_CHARS_LEN {
                //change the last symbol, then do the counts    
                //curr_buffer[curr_offset+kmer_len-1] = VALID_CHARS[x];
                //counts[x] += bwt.count_kmer(&curr_buffer[curr_offset..curr_offset+kmer_len]);
                counts[x] = fw_counts[x] + rev_counts[x];
                if counts[x] > counts[max_pos] {
                    max_pos = x;
                }
            }
            
            //check if the best is good enough
            if counts[max_pos] < min_count {
                //its not, this is a dead-end bridge
                break;
            }

            //the best is good enough, time to see just how many are good enough
            curr_bridge.push(4);
            for x in 0..VALID_CHARS_LEN {
                if x != max_pos && counts[x] >= min_count {
                    curr_bridge[curr_bridge_len] = VALID_CHARS[x];
                    possible_bridges.push(curr_bridge.clone());
                }
            }

            //check if the branch limit has been hit
            if possible_bridges.len() >= branch_limit {
                return Vec::<Vec<u8>>::new();
            }

            //now finish out the main path
            curr_bridge[curr_bridge_len] = VALID_CHARS[max_pos];
            curr_bridge_len += 1;
            
            //extend k-mer buffers with the best character
            curr_buffer[curr_offset+kmer_len-1] = VALID_CHARS[max_pos];
            rev_buffer[curr_offset+kmer_len-1] = string_util::COMPLEMENT_INT[VALID_CHARS[max_pos] as usize];

            //check if we found the target
            if &curr_buffer[curr_offset..curr_offset+kmer_len] == target_kmer {
                //add this bridge
                ret.push(curr_bridge.clone());
                
                //check if we have too many
                if ret.len() >= branch_limit {
                    return Vec::<Vec<u8>>::new();
                }
            }
        }
    }

    //return the list of found bridges
    if num_branched < branch_limit {
        ret
    }
    else {
        Vec::<Vec<u8>>::new()
    }
}

/// Searches for a sequence that is supported by the `bwt` and requires the initial final k-mers to match.
/// # Arguments
/// * `bwt` - the BWT that contains count data
/// * `sequence` - the integer form sequence we want to replace
/// * `kmer_len` - the length of the k-mer to use for counting
/// * `min_count` - the minimum count required for a path to be considered solid
/// * `max_count` - the maximum count allowed for a path to be considered solit
/// * `branch_limit` - the maximum number of branches to explore before giving up on all paths
/// * `max_branch_len` - the maximum branch length allowed before giving up on a path
pub fn bridge_sequence(
    bwt: &BitVectorBWT, sequence: &[u8], kmer_len: usize, min_count: u64, max_count: u64,
    branch_limit: usize, max_branch_len: usize
) -> Vec<Vec<u8>> {
    //build up return
    let mut ret = Vec::<Vec<u8>>::new();

    //build some helper values we'll be looping through a lot
    let seed_kmer = &sequence[..kmer_len];
    let target_kmer = &sequence[sequence.len()-kmer_len..];
    assert_eq!(kmer_len, target_kmer.len());
    let mut counts: [u64; VALID_CHARS_LEN] = [0; VALID_CHARS_LEN];
    let mut fw_counts: [u64; VALID_CHARS_LEN] = [0; VALID_CHARS_LEN];
    let mut rev_counts: [u64; VALID_CHARS_LEN] = [0; VALID_CHARS_LEN];
    let mut num_branched: usize = 0;
    let mut total_found: usize;
    let mut only_pos: usize = 0;

    //these buffers are used to create query slices
    let mut curr_buffer: Vec<u8> = vec![4; max_branch_len];
    let mut rev_buffer: Vec<u8> = vec![4; max_branch_len];

    //initialize the bridging with our seed k-mer
    let mut possible_bridges = BinaryHeap::new();
    let mut seed_vec: Vec<u8> = vec![0; kmer_len];
    seed_vec.clone_from_slice(seed_kmer);
    let mut initial_bridge = DynamicWFA::new(sequence);
    for c in sequence[..kmer_len].iter() {
        initial_bridge.append(*c);
    }
    possible_bridges.push(initial_bridge);

    //let mut max_edit_distance: usize = (0.30*max_branch_len as f64) as usize; 
    let mut max_edit_distance: usize = max_branch_len - sequence.len();
    
    while !possible_bridges.is_empty() && num_branched < branch_limit {
        //get a bridge to extend
        //let Reverse((curr_ed, mut curr_bridge)) = possible_bridges.pop().unwrap();
        let mut curr_bridge: DynamicWFA = possible_bridges.pop().unwrap();
        
        //if this ED is already worse than our best, just skip this one
        if curr_bridge.get_edit_distance() > max_edit_distance {
            continue;
        }

        //get some helper values ready
        let cb_other_seq = curr_bridge.get_other_seq();
        let mut curr_bridge_len = cb_other_seq.len();
        let mut curr_offset: usize = 0;
        num_branched += 1;

        //TODO: replace this with copy slice? not sure how to do the rev comp in a one-line
        for x in 0..kmer_len {
            curr_buffer[x] = cb_other_seq[curr_bridge_len-kmer_len+x];
            rev_buffer[x] = string_util::COMPLEMENT_INT[curr_buffer[x] as usize];
        }

        while curr_bridge_len < max_branch_len && curr_bridge.get_edit_distance() <= max_edit_distance {
            //increase the offset into our buffers (i.e. shift the curr k-mer left and rev k-mer right)
            curr_offset += 1;
            
            //do all the k-mer counting, efficient on rev-comp, then forward queries are added in
            bwt.prefix_revkmer_noalloc_fixed(&rev_buffer[curr_offset..curr_offset+kmer_len-1], &mut rev_counts);
            
            //parallel forward query, doesn't seem to gain much if anything
            bwt.postfix_kmer_noalloc_fixed(&curr_buffer[curr_offset..curr_offset+kmer_len-1], &mut fw_counts);
            
            //figure out how many extensions exist
            total_found = 0;
            for x in 0..VALID_CHARS_LEN {
                //change the last symbol, then do the counts    
                counts[x] = fw_counts[x] + rev_counts[x];
                if counts[x] >= min_count && counts[x] <= max_count {
                    total_found += 1;
                    only_pos = x;
                }
            }
            
            match total_found {
                0 => {
                    //no solid extensions, break out of the loop
                    break;
                },
                1 => {
                    //only one extension, append it
                    curr_bridge.append(VALID_CHARS[only_pos]);
                    curr_bridge_len += 1;

                    //extend k-mer buffers with the best character
                    curr_buffer[curr_offset+kmer_len-1] = VALID_CHARS[only_pos];
                    rev_buffer[curr_offset+kmer_len-1] = string_util::COMPLEMENT_INT[VALID_CHARS[only_pos] as usize];

                    //check if we found the target
                    if &curr_buffer[curr_offset..curr_offset+kmer_len] == target_kmer {
                        if curr_bridge.get_max_distance() >= sequence.len() {
                            //the WFA reaches the end of the primary sequence, check the ED results
                            if curr_bridge.get_edit_distance() < max_edit_distance {
                                //bridge has lower ED, clear other results then push this one
                                max_edit_distance = curr_bridge.get_edit_distance();
                                ret.clear();
                                ret.push(curr_bridge.get_other_seq().clone());
                            } else if curr_bridge.get_edit_distance() == max_edit_distance {
                                //bridge has same ED, just add to list
                                ret.push(curr_bridge.get_other_seq().clone());
                            }
                        } else {
                            //WFA doesn't reach the end, clone it and finalize it
                            let mut finalized_bridge = curr_bridge.clone();
                            finalized_bridge.finalize();
                            assert!(finalized_bridge.get_max_distance() >= sequence.len());

                            //now check the finalized ED and do appropriate things
                            if finalized_bridge.get_edit_distance() < max_edit_distance {
                                //bridge has lower ED, clear other results then push this one
                                max_edit_distance = finalized_bridge.get_edit_distance();
                                ret.clear();
                                ret.push(finalized_bridge.get_other_seq().clone());
                            } else if finalized_bridge.get_edit_distance() == max_edit_distance {
                                //bridge has same ED, just add to list
                                ret.push(finalized_bridge.get_other_seq().clone());
                            }
                        }
                    }
                },
                _tf => {
                    //multiple are found, push all of them as new candidate for extension
                    for x in 0..VALID_CHARS_LEN {
                        if counts[x] >= min_count && counts[x] <= max_count {
                            //clone the bridge and then append the extension that was solid
                            let mut new_bridge = curr_bridge.clone();
                            new_bridge.append(VALID_CHARS[x]);
                            
                            //now we need to check if this extension matches our target and then push as appropriate
                            curr_buffer[curr_offset+kmer_len-1] = VALID_CHARS[x];
                            if &curr_buffer[curr_offset..curr_offset+kmer_len] == target_kmer {
                                if new_bridge.get_max_distance() >= sequence.len() {
                                    //the WFA reaches the end of the primary sequence, check the ED results
                                    if new_bridge.get_edit_distance() < max_edit_distance {
                                        //bridge has lower ED, clear other results then push this one
                                        max_edit_distance = new_bridge.get_edit_distance();
                                        ret.clear();
                                        ret.push(new_bridge.get_other_seq().clone());
                                    } else if new_bridge.get_edit_distance() == max_edit_distance {
                                        //bridge has same ED, just add to list
                                        ret.push(new_bridge.get_other_seq().clone());
                                    }
                                } else {
                                    //WFA doesn't reach the end, clone it and finalize it
                                    let mut finalized_bridge = new_bridge.clone();
                                    finalized_bridge.finalize();
                                    assert!(finalized_bridge.get_max_distance() >= sequence.len());
        
                                    //now check the finalized ED and do appropriate things
                                    if finalized_bridge.get_edit_distance() < max_edit_distance {
                                        //bridge has lower ED, clear other results then push this one
                                        max_edit_distance = finalized_bridge.get_edit_distance();
                                        ret.clear();
                                        ret.push(finalized_bridge.get_other_seq().clone());
                                    } else if finalized_bridge.get_edit_distance() == max_edit_distance {
                                        //bridge has same ED, just add to list
                                        ret.push(finalized_bridge.get_other_seq().clone());
                                    }
                                }
                            }

                            //finally add it to the extension list
                            possible_bridges.push(new_bridge);
                        }
                    }

                    //end this loop because we've added each one as an extension to the heap
                    break;
                }
            };
        }
    }

    //return the list of found bridges
    if num_branched < branch_limit {
        ret
    }
    else {
        Vec::<Vec<u8>>::new()
    }
}

/// Given an initial k-mer, this will extend outwards from that k-mer up to the max branch length.
/// # Arguments
/// * `bwt` - the BWT source for k-mer counts
/// * `seed_kmer` - the initial seed k-mer in integer form
/// * `min_count` - the minimum count to consider a k-mer as present
/// * `branch_limit` - the maximum allowed number of branches
/// * `max_branch_len` - the maximum length of any given branch
pub fn assemble_from_kmer(
    bwt: &BitVectorBWT, seed_kmer: &[u8], min_count: u64, branch_limit: usize, max_branch_len: usize
) -> Vec<Vec<u8>> {
    //build up return
    let mut ret = Vec::<Vec<u8>>::new();

    //build some helper values we'll be looping through a lot
    let kmer_len = seed_kmer.len();
    let mut counts: [u64; VALID_CHARS_LEN] = [0; VALID_CHARS_LEN];
    let mut fw_counts: [u64; VALID_CHARS_LEN] = [0; VALID_CHARS_LEN];
    let mut rev_counts: [u64; VALID_CHARS_LEN] = [0; VALID_CHARS_LEN];
    let mut num_branched: usize = 0;
    let mut max_pos: usize;

    //these buffers are used to create query slices
    let mut curr_buffer: Vec<u8> = vec![4; max_branch_len];
    let mut rev_buffer: Vec<u8> = vec![4; max_branch_len];

    //initialize the bridging with our seed k-mer
    let mut possible_bridges: Vec<Vec<u8>> = Vec::<Vec<u8>>::new();
    let mut seed_vec: Vec<u8> = vec![0; kmer_len];
    seed_vec.clone_from_slice(seed_kmer);
    possible_bridges.push(seed_vec);

    while !possible_bridges.is_empty() && num_branched < branch_limit {
        //get a bridge to extend
        let mut curr_bridge: Vec<u8> = possible_bridges.pop().unwrap();
        let mut curr_bridge_len = curr_bridge.len();
        let mut curr_offset: usize = 0;
        num_branched += 1;

        //TODO: replace this with copy slice? not sure how to do the rev comp in a one-line
        for x in 0..kmer_len {
            curr_buffer[x] = curr_bridge[curr_bridge_len-kmer_len+x];
            rev_buffer[x] = string_util::COMPLEMENT_INT[curr_buffer[x] as usize];
        }

        while curr_bridge_len < max_branch_len {
            //increment the offset
            curr_offset += 1;

            //do all the k-mer counting, efficient on rev-comp, then forward queries are added in
            bwt.prefix_revkmer_noalloc_fixed(&rev_buffer[curr_offset..curr_offset+kmer_len-1], &mut rev_counts);
            
            //parallel forward query, not much faster but there is some gain
            bwt.postfix_kmer_noalloc_fixed(&curr_buffer[curr_offset..curr_offset+kmer_len-1], &mut fw_counts);
            
            max_pos=0;
            for x in 0..VALID_CHARS_LEN {
                //change the last symbol, then do the counts    
                counts[x] = fw_counts[x] + rev_counts[x];
                if counts[x] > counts[max_pos] {
                    max_pos = x;
                }
            }
            
            //check if the best is good enough
            if counts[max_pos] < min_count {
                //its not, this is a dead-end bridge
                break;
            }

            //the best is good enough, time to see just how many are good enough
            curr_bridge.push(4);
            for x in 0..VALID_CHARS_LEN {
                if x != max_pos && counts[x] >= min_count {
                    curr_bridge[curr_bridge_len] = VALID_CHARS[x];
                    possible_bridges.push(curr_bridge.clone());
                }
            }

            //check if the branch limit has been hit
            if possible_bridges.len() >= branch_limit {
                return Vec::<Vec<u8>>::new();
            }

            //now finish out the main path
            curr_bridge[curr_bridge_len] = VALID_CHARS[max_pos];
            curr_bridge_len += 1;
            
            //extend k-mer buffers with the best character
            curr_buffer[curr_offset+kmer_len-1] = VALID_CHARS[max_pos];
            rev_buffer[curr_offset+kmer_len-1] = string_util::COMPLEMENT_INT[VALID_CHARS[max_pos] as usize];
        }

        //TODO: we should revisit this requirement to be as long as the max_branch_len, maybe add a min_branch_len?
        if curr_bridge_len == max_branch_len {
            //we hit the limit, save the result and return out
            ret.push(curr_bridge);
            if ret.len() >= branch_limit {
                return Vec::<Vec<u8>>::new();
            }
        }
    }

    //return the list of found bridges
    if num_branched < branch_limit {
        ret
    }
    else {
        Vec::<Vec<u8>>::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::bwt_converter::convert_to_vec;
    use crate::ropebwt2_util::create_bwt_from_strings;
    use crate::string_util::{convert_stoi, reverse_complement_i};
    
    fn get_constant_bwt() -> BitVectorBWT {
        //build the dataset
        let const_string = "AACGGATCAAGCTTACCAGTATTTACGT";
        let rep_count = 30;
        let mut data: Vec<&str> = vec![];
        for _i in 0..rep_count {
            data.push(&const_string);
        }

        //build our BWT
        let bwt_string = create_bwt_from_strings(&data).unwrap();
        let bwt_cursor = Cursor::new(bwt_string);
        let vec = convert_to_vec(bwt_cursor);
        let mut bv_bwt: BitVectorBWT = BitVectorBWT::new();
        bv_bwt.load_vector(vec);
        bv_bwt
    }

    fn get_diploid_bwt() -> BitVectorBWT {
        //build the dataset
        let const_string = "AACGGATCAAGCTTACCAGTATTTACGT";
        let rep_count = 30;
        let mut data: Vec<&str> = vec![];
        for _i in 0..rep_count {
            data.push(&const_string);
        }

        //this less string is reverse-complemented so we can verify RC counting & such is working from a correction standpoint
        //                   AACGGATCAAGCTTACCAGTATTTACGT
        //                      #         #            #
        let lesser_string = "AACTGATCAAGCTGACCAGTATTTACTT"; //3 basepair changes
        let lesser_string = string_util::convert_itos(
            &string_util::reverse_complement_i(
                &string_util::convert_stoi(&lesser_string)
            )
        );
        let lesser_count = 20;
        for _i in 0..lesser_count {
            data.push(&lesser_string);
        }

        //build our BWT
        let bwt_string = create_bwt_from_strings(&data).unwrap();
        let bwt_cursor = Cursor::new(bwt_string);
        let vec = convert_to_vec(bwt_cursor);
        let mut bv_bwt: BitVectorBWT = BitVectorBWT::new();
        bv_bwt.load_vector(vec);
        bv_bwt
    }

    fn get_expansion_bwt() -> BitVectorBWT {
        //build the dataset
        let const_string = "AACGGATACACACACACACACTTTACGT";
        let rep_count = 30;
        let mut data: Vec<&str> = vec![];
        for _i in 0..rep_count {
            data.push(&const_string);
        }

        //build our BWT
        let bwt_string = create_bwt_from_strings(&data).unwrap();
        let bwt_cursor = Cursor::new(bwt_string);
        let vec = convert_to_vec(bwt_cursor);
        let mut bv_bwt: BitVectorBWT = BitVectorBWT::new();
        bv_bwt.load_vector(vec);
        bv_bwt
    }

    #[test]
    fn test_bridge() {
        //build our test and verify it's fine
        let bwt = get_constant_bwt();
        let query = convert_stoi(&"AACGGATCAAGCTTACCAGTATTTACGT");
        assert_eq!(bwt.count_kmer(&query), 30);

        //okay, now test the actual bridging
        let seed = convert_stoi(&"AACGGAT");
        let target = convert_stoi(&"TTTACGT");
        let min_count = 15;
        let branch_lim = 20;
        let max_branch_len = 40;
        let bridges = bridge_kmers(&bwt, &seed, &target, min_count, branch_lim, max_branch_len);
        assert_eq!(bridges.len(), 1);
        assert_eq!(bridges[0], query);
        
        //now do it in reverse complement space
        let rev_seed = reverse_complement_i(&target);
        let rev_target = reverse_complement_i(&seed);
        let rev_bridges = bridge_kmers(&bwt, &rev_seed, &rev_target, min_count, branch_lim, max_branch_len);
        assert_eq!(rev_bridges.len(), 1);
        assert_eq!(rev_bridges[0], reverse_complement_i(&query));
    }

    #[test]
    fn test_bridge_sequence() {
        //build our test and verify it's fine
        let bwt = get_constant_bwt();
        let query =     convert_stoi(&"AACGGATCAAGCTTACCAGTATTTACGT");
        let error_seq = convert_stoi(&"AACGGATGAAGTTAACCAAGTTTTACGT");
        assert_eq!(bwt.count_kmer(&query), 30);

        //okay, now test the actual bridging
        let seed = convert_stoi(&"AACGGAT");
        let target = convert_stoi(&"TTTACGT");
        let min_count = 15;
        let max_count = 1000;
        let branch_lim = 20;
        let max_branch_len = 40;
        let bridges = bridge_sequence(&bwt, &error_seq, seed.len(), min_count, max_count, branch_lim, max_branch_len);
        assert_eq!(bridges.len(), 1);
        assert_eq!(bridges[0], query);
        
        //now do it in reverse complement space
        let rev_seed = reverse_complement_i(&target);
        let rev_bridges = bridge_sequence(&bwt, &reverse_complement_i(&error_seq), rev_seed.len(), min_count, max_count, branch_lim, max_branch_len);
        assert_eq!(rev_bridges.len(), 1);
        assert_eq!(rev_bridges[0], reverse_complement_i(&query));
    }

    #[test]
    fn test_expansion_bridge() {
        //build our test and verify it's fine
        let bwt = get_expansion_bwt();
        let query = convert_stoi(&"AACGGATACACACACACACACTTTACGT");
        assert_eq!(bwt.count_kmer(&query), 30);

        //okay, now test the actual bridging
        let seed = convert_stoi(&"AACGGAT");
        let target = convert_stoi(&"TTTACGT");
        let min_count = 15;
        let branch_lim = 20;
        let max_branch_len = query.len();
        let bridges = bridge_kmers(&bwt, &seed, &target, min_count, branch_lim, max_branch_len);

        //due to the repeat, you get a bunch of possible up to the length of the query
        assert_eq!(bridges.len(), 5);
        for bridge in bridges.iter() {
            //make sure each one starts and ends with the correct sequence
            assert_eq!(bridge[0..seed.len()], seed[..]);
            assert_eq!(bridge[bridge.len()-target.len()..bridge.len()], target[..]);
        }
    }

    #[test]
    fn test_assemble() {
        //build our test and verify it's fine
        let bwt = get_constant_bwt();
        let query = convert_stoi(&"AACGGATCAAGCTTACCAGTATTTACGT");
        assert_eq!(bwt.count_kmer(&query), 30);

        //okay, now test the actual bridging
        let seed = convert_stoi(&"AACGGAT");
        let min_count = 15;
        let branch_lim = 20;
        let max_branch_len = query.len(); //fix the max length to the strings in the BWT
        let bridges = assemble_from_kmer(&bwt, &seed, min_count, branch_lim, max_branch_len);
        assert_eq!(bridges.len(), 1);
        assert_eq!(bridges[0], query);
    }

    /*
    //originally for testing a version that allowed for one-way assembly, but it seemed to be very slow and give poor results
    #[test]
    fn test_assemble_bridge() {
        //build our test and verify it's fine
        let bwt = get_constant_bwt();
        let query =     convert_stoi(&"AACGGATCAAGCTTACCAGTATTTACGT");
        let error_seq = convert_stoi(&"AACGGATCAAAGCTACCACTATTACAT");
        //let query =     convert_stoi(&"AACGGATC");
        //let error_seq = convert_stoi(&"AACGGATT");
        assert_eq!(bwt.count_kmer(&query), 30);

        //okay, now test the actual bridging
        let min_count = 15;
        let max_count = 1000;
        let branch_lim = 20;
        let max_branch_len = query.len()+10; //fix the max length to the strings in the BWT
        let bridges = bridge_sequence(&bwt, &error_seq, 7, min_count, max_count, branch_lim, max_branch_len, true);
        assert_eq!(bridges.len(), 1);
        assert_eq!(bridges[0], query);
    }
    */

    #[test]
    fn test_triple_accel() {
        //this test is mostly just for my sanity of understanding the library
        let query = convert_stoi(&"AACGGATCAAGCTTACCAGTATTTACGT");
        let short_seq = convert_stoi(&"AACGGATGAA"); //one base change
        //let long_seq = convert_stoi(&"AACGGATCATAGCTTACCAGTATATACGT"); //one insert, one change

        //use this one when we want to compare the full strings and get an edit distance
        //let l_dist = levenshtein_exp(&long_seq, &query);
        //assert_eq!(l_dist, 2);

        //use this one when we need to do minimal head/tail corrections
        let matches: Vec<Match> = levenshtein_search(&short_seq, &query).collect();
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0], Match{start: 0, end: 10, k: 1})
    }

    #[test]
    fn test_correction_pass() {
        //build our test and verify it's fine
        let bwt = get_constant_bwt();
        let query = convert_stoi(&"AACGGATCAAGCTTACCAGTATTTACGT");
        assert_eq!(bwt.count_kmer(&query), 30);
        
        //now lets alter a sequence and test it
        let params = CorrectionParameters {
            kmer_sizes: vec![9],
            min_count: 5,
            max_branch_attempt_length: 10000,
            branch_limit_factor: 4.0,
            branch_buffer_factor: 1.3,
            midpoint_ed_factor: 0.4,
            tail_buffer_factor: 1.00, //normally - 1.05,
            frac: 0.1,
            verbose: true
        };

        //single head SNV change
        let change_head_seq = convert_stoi(&"AATGGATCAAGCTTACCAGTATTTACGT");
        let corrected = correction_pass(&bwt, &change_head_seq, &params, 9);
        assert_eq!(corrected, query);

        //truncate head test - removed 2 characters at start
        let change_head_seq = convert_stoi(&"TGGATCAAGCTTACCAGTATTTACGT");
        let corrected = correction_pass(&bwt, &change_head_seq, &params, 9);
        assert_eq!(corrected[..], query[2..]);

        //single bridge deletion
        let change_mid_seq = convert_stoi(&"AACGGATCAAGCTTCCAGTATTTACGT");
        let corrected = correction_pass(&bwt, &change_mid_seq, &params, 9);
        assert_eq!(corrected, query);
        
        //single tail <EVENT>
        let change_tail_seq = convert_stoi(&"AACGGATCAAGCTTACCAGTATTTGCGT");
        let corrected = correction_pass(&bwt, &change_tail_seq, &params, 9);
        assert_eq!(corrected, query);

        //truncated tail test - removed 3 characters at end
        let change_tail_seq = convert_stoi(&"AACGGATCAAGCTTACCAGTATTTG");
        let corrected = correction_pass(&bwt, &change_tail_seq, &params, 9);
        assert_eq!(corrected[..], query[..query.len()-3]);
    }

    #[test]
    fn test_diploid_correction_pass() {
        let bwt = get_diploid_bwt();
        //copied from the diploid bwt
        let const_string  = convert_stoi(&"AACGGATCAAGCTTACCAGTATTTACGT");
        let lesser_string = convert_stoi(&"AACTGATCAAGCTGACCAGTATTTACTT"); //3 basepair changes

        //now lets alter a sequence and test it
        let params = CorrectionParameters {
            kmer_sizes: vec![9],
            min_count: 5,
            max_branch_attempt_length: 10000,
            branch_limit_factor: 4.0,
            branch_buffer_factor: 1.3,
            midpoint_ed_factor: 0.4,
            tail_buffer_factor: 1.00, //normally - 1.05,
            frac: 0.1,
            verbose: true
        };
        
        //HEAD TESTS
        //should match const string
        //                                   AACGGATCAAGCTTACCAGTATTTACGT
        let change_head_seq = convert_stoi(&"ATCGGATCAAGCT");
        let corrected = correction_pass(&bwt, &change_head_seq, &params, 9);
        assert_eq!(&corrected[..], &const_string[0..corrected.len()]);

        //should match lesser string
        //                                   AACTGATCAAGCTGACCAGTATTTACTT
        let change_head_seq = convert_stoi(&"AAGTGATCAAGCT");
        let corrected = correction_pass(&bwt, &change_head_seq, &params, 9);
        assert_eq!(&corrected[..], &lesser_string[0..corrected.len()]);

        //could match either one, so it will match the greater
        let change_head_seq = convert_stoi(&"AACAGATCAAGCT");
        let corrected = correction_pass(&bwt, &change_head_seq, &params, 9);
        assert_eq!(&corrected[..], &const_string[0..corrected.len()]);

        //BRIDGE TESTS
        //should match const string
        let change_mid_seq = convert_stoi(&"AACGGATCAAGCTTAGCAGTATTTACGT");
        let corrected = correction_pass(&bwt, &change_mid_seq, &params, 9);
        assert_eq!(&corrected, &const_string);

        //should match less string
        let change_mid_seq = convert_stoi(&"AACTGATCAAGCTGAGCAGTATTTACTT");
        let corrected = correction_pass(&bwt, &change_mid_seq, &params, 9);
        assert_eq!(&corrected, &lesser_string);

        //the SNP event is changes to "A", so should match const string
        let change_mid_seq = convert_stoi(&"GATCAAGCTAACCAGTATT");
        let corrected = correction_pass(&bwt, &change_mid_seq, &params, 9);
        assert_eq!(&corrected, &convert_stoi(&"GATCAAGCTTACCAGTATT"));

        //TAIL TESTS
        //should match const string
        let change_tail_seq = convert_stoi(&"CAGTATTTAGGT");
        let corrected = correction_pass(&bwt, &change_tail_seq, &params, 9);
        assert_eq!(&corrected[..], &const_string[const_string.len()-corrected.len()..]);

        //should match less string
        let change_tail_seq = convert_stoi(&"CAGTATTTAGTT");
        let corrected = correction_pass(&bwt, &change_tail_seq, &params, 9);
        assert_eq!(&corrected[..], &lesser_string[lesser_string.len()-corrected.len()..]);
        
        //could match either, so should match const
        let change_tail_seq = convert_stoi(&"CAGTATTTACAT");
        let corrected = correction_pass(&bwt, &change_tail_seq, &params, 9);
        assert_eq!(&corrected[..], &const_string[const_string.len()-corrected.len()..]);
    }

    #[test]
    fn test_correction_job() {
        //shared bwt
        let bwt: BitVectorBWT = get_constant_bwt();
        let arc_bwt: Arc<BitVectorBWT> = Arc::new(bwt);

        //shared params
        let params = CorrectionParameters {
            kmer_sizes: vec![9],
            min_count: 5,
            max_branch_attempt_length: 10000,
            branch_limit_factor: 4.0,
            branch_buffer_factor: 1.3,
            midpoint_ed_factor: 0.4,
            tail_buffer_factor: 1.00, //normally - 1.05,
            frac: 0.1,
            verbose: true
        };
        let arc_params: Arc<CorrectionParameters> = Arc::new(params);

        //the one string in the bwt
        let const_string: String =  "AACGGATCAAGCTTACCAGTATTTACGT".to_string();
        
        let long_read_const_mod: LongReadFA = LongReadFA {
            read_index: 0,
            label: "test".to_string(),
            seq: "TACGGATCAAGCATACCAGTATGTACGT".to_string()
        };
        let corr_results = correction_job(arc_bwt, long_read_const_mod.clone(), arc_params);
        assert_eq!(corr_results.label, long_read_const_mod.label);
        assert_eq!(corr_results.original_seq, long_read_const_mod.seq);
        assert_eq!(corr_results.corrected_seq, const_string);
        assert_eq!(corr_results.avg_before, 6.0);
        assert_eq!(corr_results.avg_after, 30.0);
    }

    #[test]
    fn test_correction_deletion_special() {
        //shared bwt
        let bwt: BitVectorBWT = get_constant_bwt();
        let arc_bwt: Arc<BitVectorBWT> = Arc::new(bwt);

        //shared params
        let params = CorrectionParameters {
            kmer_sizes: vec![9],
            min_count: 5,
            max_branch_attempt_length: 10000,
            branch_limit_factor: 4.0,
            branch_buffer_factor: 1.3,
            midpoint_ed_factor: 0.4,
            tail_buffer_factor: 1.00, //normally - 1.05,
            frac: 0.1,
            verbose: true
        };
        let arc_params: Arc<CorrectionParameters> = Arc::new(params);

        //the one string in the bwt
        let const_string: String =  "AACGGATCAAGCTTACCAGTATTTACGT".to_string();
        
        //this specific error is cause by a poly-X insertion between two adjacent X-mers
        //here we just inserted a "T" in to the "TT" run and it triggers a panic due to the final string building
        let long_read_const_mod: LongReadFA = LongReadFA {
            read_index: 10,
            label: "test".to_string(),
            seq: "AACGGATCAAGCTTTACCAGTATTTACGT".to_string()
        };
        let corr_results = correction_job(arc_bwt, long_read_const_mod.clone(), arc_params);
        assert_eq!(corr_results.read_index, long_read_const_mod.read_index);
        assert_eq!(corr_results.label, long_read_const_mod.label);
        assert_eq!(corr_results.original_seq, long_read_const_mod.seq);
        assert_eq!(corr_results.corrected_seq, const_string);
        
    }
}