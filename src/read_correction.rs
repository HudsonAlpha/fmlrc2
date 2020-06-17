
use std::sync::Arc;
use triple_accel::{levenshtein,levenshtein_search,Match};

use crate::bv_bwt::BitVectorBWT;
use crate::stats_util;
use crate::string_util;

const VALID_CHARS: [u8; 4] = [1, 2, 3, 5];
const VALID_CHARS_LEN: usize = VALID_CHARS.len();

/// stores options for running the correction algorithms
pub struct CorrectionParameters {
    //use_fm_index: bool, //old parameter that toggled between bit vector and classic mode
    kmer_sizes: Vec<usize>,
    min_count: u64,
    max_branch_attempt_length: usize,
    branch_limit_factor: u64,
    branch_buffer_factor: f64,
    //TODO: add a tail_truncate_factor that buts a bounding box around min length and max length
    tail_buffer_factor: f64,
    frac: f64,
    //fm_bit_power: u8, //only matters for classic mode which isn't implemented currently
    verbose: bool
}

/// a basic struct for storing a correction to a sequence
pub struct Correction {
    start_pos: usize,
    end_pos: usize,
    seq: Vec<u8>
}

/// a struct for storing generic read
#[derive(Clone)]
pub struct LongReadFA {
    label: String,
    seq: String
}

/// a struct for storing the modified string
pub struct CorrectionResults {
    label: String,
    original_seq: String,
    corrected_seq: String,
    avg_before: f64,
    avg_after: f64
}

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

    let avg_before: f64;
    let avg_after: f64;
    if params.verbose {
        //TODO: is there a way to pre-serve this from earlier? this mode will always take longer, so it's not a huge deal
        let orig_seq: Vec<u8> = string_util::convert_stoi(&long_read.seq);
        let counts_before: Vec<u64> = bwt.count_pileup(&orig_seq, params.kmer_sizes[0]);
        let counts_after: Vec<u64> = bwt.count_pileup(&seq_i, params.kmer_sizes[0]);
        let sum_before: u64 = counts_before.iter().sum();
        let sum_after: u64 = counts_after.iter().sum();
        avg_before = sum_before as f64 / counts_before.len() as f64;
        avg_after = sum_after as f64 / counts_after.len() as f64;
    }
    else {
        avg_before = 0.0;
        avg_after = 0.0;
    }

    //send it on back y'all
    CorrectionResults {
        label: long_read.label,
        original_seq: long_read.seq,
        corrected_seq: corrected_seq,
        avg_before: avg_before,
        avg_after: avg_after
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
pub fn correction_pass(bwt: &BitVectorBWT, seq_i: &Vec<u8>, params: &CorrectionParameters, kmer_size: usize) -> Vec<u8> {
    //count up the initial pileups
    let pileup: Vec<u64> = bwt.count_pileup(seq_i, kmer_size);
    let nz_med = stats_util::calculate_bounded_median(&pileup, params.min_count);
    if nz_med < params.min_count {
        return seq_i.clone();
    }
    let pileup_size: usize = pileup.len();
    
    //we can try some fixing, set up additional params for this pass
    let branch_limit = params.branch_limit_factor as usize*kmer_size;
    
    //try to dynamically set the threshold, but make sure its at least MIN_COUNT
    let threshold: u64 = std::cmp::max((params.frac * nz_med as f64) as u64, params.min_count);

    //prep for the actual corrections now
    let mut prev_found: isize = -1;
    let mut seed_kmer: Vec<u8> = vec![0; kmer_size];
    let mut target_kmer: Vec<u8> = vec![0; kmer_size];
    let mut max_branch_length: usize;

    let mut bridge_points: Vec<Vec<u8>>;
    //let mut bridge_points_ed: Vec<Vec<u8>> = Vec::<Vec<u8>>::new();
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
                    let best_match: Option<Vec<u8>> = pick_best_levenshtein_search(&orig, bridge_points, bwt, kmer_size);

                    match best_match {
                        Some(bm) => {
                            //we found a best match, it will span from index 0 up to the start of the first found k-mer
                            let truncated_seq: Vec<u8> = string_util::reverse_complement_i(&bm)[0..bm.len()-kmer_size].to_vec();
                            
                            //now store the correction in range [0..x)
                            new_corr = Correction {
                                start_pos: 0,
                                end_pos: x,
                                seq: truncated_seq
                            };
                            corrections_list.push(new_corr);
                        },
                        None => {
                            //no match was found that works
                        }
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
                    bridge_points = bridge_kmers(bwt, &seed_kmer, &target_kmer, threshold, branch_limit, max_branch_length);

                    //try reverse complement if we failed
                    if bridge_points.len() == 0 {
                        bridge_points = bridge_kmers(
                            bwt, 
                            &string_util::reverse_complement_i(&target_kmer),
                            &string_util::reverse_complement_i(&seed_kmer),
                            threshold,
                            branch_limit,
                            max_branch_length
                        );
                        
                        //make sure to rev-comp the results here
                        for y in {0..bridge_points.len()} {
                            bridge_points[y] = string_util::reverse_complement_i(&bridge_points[y]);
                        }
                    }
                    
                    //we are doing a bridge, so use straight levenshtein
                    let best_match: Option<Vec<u8>> = pick_best_levenshtein(&seq_i[prev_found as usize..x+kmer_size], bridge_points, bwt, kmer_size);

                    //pick out the best result from the full levenshtein
                    match best_match {
                        Some(bm) => {
                            //we found a best match, it will span from k-mer to k-mer, but clipping "k" from each end
                            let truncated_seq: Vec<u8> = bm[kmer_size..bm.len()-kmer_size].to_vec();

                            //now store the correction in range [0..x)
                            new_corr = Correction {
                                start_pos: prev_found as usize+kmer_size,
                                end_pos: x,
                                seq: truncated_seq
                            };
                            corrections_list.push(new_corr);
                        },
                        None => {
                            //no match was found that works
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
        max_branch_length = (params.tail_buffer_factor*(pileup_size-1-prev_found as usize+kmer_size) as f64) as usize;
        if max_branch_length <= params.max_branch_attempt_length {
            //copy the last found k-mer and assemble outwards
            seed_kmer[..].clone_from_slice(&seq_i[prev_found as usize..prev_found as usize+kmer_size]);
            bridge_points = assemble_from_kmer(bwt, &seed_kmer, threshold, branch_limit, max_branch_length);
            
            //now get the best match
            let best_match: Option<Vec<u8>> = pick_best_levenshtein_search(&seq_i[prev_found as usize..], bridge_points, bwt, kmer_size);
            match best_match {
                Some(bm) => {
                    //we found a best match, it will span from index "prev_found+kmer_size" through the end
                    let truncated_seq: Vec<u8> = bm[kmer_size..].to_vec();
                    
                    //now store the correction in range [0..x)
                    new_corr = Correction {
                        start_pos: prev_found as usize+kmer_size,
                        end_pos: seq_i.len(),
                        seq: truncated_seq
                    };
                    corrections_list.push(new_corr);
                },
                None => {
                    //no match was found that works
                }
            }
        }
    }

    let mut current_position: usize = 0;
    let mut corrected_seq: Vec<u8> = Vec::<u8>::new();
    for correction in corrections_list {
        //TODO: is there a better way than extend_from_slice?
        //copy everything up to the correction
        corrected_seq.extend_from_slice(&seq_i[current_position..correction.start_pos]);

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
/// TODO: figure out if the checks occuring here are slowing us down, and identity ways to safely short-circuit
/// # Arguments
/// * `original` - the original sequence in integer format
/// * `candidates` - the Vec of candidates, each in integer format
/// * `bwt` - the BWT of counts (used for pileup tie-breaking)
/// * `kmer_size` - the k-mer size to use for pileup tie-breaking
fn pick_best_levenshtein_search(original: &[u8], candidates: Vec<Vec<u8>>, bwt: &BitVectorBWT, kmer_size: usize) -> Option<Vec<u8>> {
    let mut ed_scores: Vec<Option<Match>> = Vec::<Option<Match>>::with_capacity(candidates.len());
    let mut min_score: u32 = 0xFFFFFFFF;
    for y in {0..candidates.len()} {
        //calculate the min distance
        let matches: Vec<Match> = levenshtein_search(&original, &candidates[y]).collect();
        let mut best_match: Option<Match> = None;
        for m in matches {
            //make sure the start is index 0
            if m.start == 0 {
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
    let mut candidates_ed: Vec<Vec<u8>> = Vec::<Vec<u8>>::with_capacity(candidates.len());
    for y in {0..candidates.len()} {
        match &ed_scores[y] {
            Some(eds) => {
                if eds.k == min_score {
                    //we have to truncate down to the match end
                    let truncated: Vec<u8> = candidates[y][0..eds.end].to_vec();
                    candidates_ed.push(truncated);
                }
            },
            None => {}
        }
    }

    if candidates_ed.len() == 0 {
        //do nothing, we didn't find anything good
        return None;
    }
    else if candidates_ed.len() == 1 {
        //only one with smallest edit distance
        return Some(candidates_ed.remove(0));
    }
    else {
        //figure out which of the ones with equal edit distance has the most counts
        let mut max_counts: u64 = 0;
        let mut max_id: usize = 0;
        let mut ed_pu: Vec<u64>;
        let mut summation: u64;
        for y in {0..candidates_ed.len()} {
            ed_pu = bwt.count_pileup(&candidates_ed[y], kmer_size);
            summation = ed_pu.iter().sum();
            if summation > max_counts {
                max_id = y;
                max_counts = summation;
            }
        }

        //now return the best candidate
        return Some(candidates_ed.remove(max_id));
    }
}

/// This function will take an original sequence and a candidate list and pick the best candidate to return
/// after performing some sanity checks on the mapping. This performs a full comparison because it's a bridge.
/// TODO: figure out if the checks occuring here are slowing us down, and identity ways to safely short-circuit
/// # Arguments
/// * `original` - the original sequence in integer format
/// * `candidates` - the Vec of candidates, each in integer format
/// * `bwt` - the BWT of counts (used for pileup tie-breaking)
/// * `kmer_size` - the k-mer size to use for pileup tie-breaking
fn pick_best_levenshtein(original: &[u8], candidates: Vec<Vec<u8>>, bwt: &BitVectorBWT, kmer_size: usize) -> Option<Vec<u8>> {
    //two short circuit points
    if candidates.len() == 0 {
        //no valid bridges were found
        return None;
    }
    else if candidates.len() == 1 {
        //only one valid bridge, so short circuit here
        return Some(candidates[0].clone());
    }
    else {
        //we have multiple values, so check for edit distance
        let mut ed_scores: Vec<u32> = Vec::<u32>::with_capacity(candidates.len());
        for y in {0..candidates.len()} {
            //calculate the min distance
            let score: u32 = levenshtein(&original, &candidates[y]);
            ed_scores.push(score);
        }
        let min_score: u32 = *ed_scores.iter().min().unwrap();

        //get everything with a good score
        let mut candidates_ed: Vec<Vec<u8>> = Vec::<Vec<u8>>::with_capacity(candidates.len());
        for y in {0..candidates.len()} {
            if ed_scores[y] == min_score {
                candidates_ed.push(candidates[y].clone());
            }
        }

        if candidates_ed.len() == 1 {
            //only one with smallest edit distance
            return Some(candidates_ed.remove(0));
        }
        else {
            //figure out which of the ones with equal edit distance has the most counts
            let mut max_counts: u64 = 0;
            let mut max_id: usize = 0;
            let mut ed_pu: Vec<u64>;
            let mut summation: u64;
            for y in {0..candidates_ed.len()} {
                ed_pu = bwt.count_pileup(&candidates_ed[y], kmer_size);
                summation = ed_pu.iter().sum();
                if summation > max_counts {
                    max_id = y;
                    max_counts = summation;
                }
            }

            //now return the best candidate
            return Some(candidates_ed.remove(max_id));
        }
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
/// # Examples
/// TODO: See unit tests for current examples.
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
    let mut num_branched: usize = 0;
    let mut max_pos: usize;

    //the queries will populates these vectors
    let mut curr_kmer: Vec<u8> = vec![4; kmer_len];
    let mut rev_kmer: Vec<u8> = vec![4; kmer_len];

    //initialize the bridging with our seed k-mer
    let mut possible_bridges: Vec<Vec<u8>> = Vec::<Vec<u8>>::new();
    let mut seed_vec: Vec<u8> = vec![0; kmer_len];
    seed_vec.clone_from_slice(seed_kmer);
    possible_bridges.push(seed_vec);

    while possible_bridges.len() > 0 && num_branched < branch_limit {
        //get a bridge to extend
        let mut curr_bridge: Vec<u8> = possible_bridges.pop().unwrap();
        let mut curr_bridge_len = curr_bridge.len();
        num_branched += 1;

        //TODO: replace this with copy slice and call to string_util::reverse_complement_i?
        for x in {0..kmer_len} {
            curr_kmer[x] = curr_bridge[curr_bridge_len-kmer_len+x];
            rev_kmer[kmer_len-x-1] = string_util::COMPLEMENT_INT[curr_kmer[x] as usize];
        }

        while curr_bridge_len < max_branch_len {
            //shift the current k-mer over one in preparation for the last base toggle
            for x in {0..kmer_len-1} {
                curr_kmer[x] = curr_kmer[x+1];
                rev_kmer[kmer_len-x-1] = rev_kmer[kmer_len-x-2];
            }
            
            //do all the k-mer counting
            max_pos = 0;
            for x in {0..VALID_CHARS_LEN} {
                curr_kmer[kmer_len-1] = VALID_CHARS[x];
                rev_kmer[0] = string_util::COMPLEMENT_INT[VALID_CHARS[x] as usize];
                counts[x] = bwt.count_kmer(&curr_kmer) + bwt.count_kmer(&rev_kmer);
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
            for x in {0..VALID_CHARS_LEN} {
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
            
            //update k-mers
            curr_kmer[kmer_len-1] = VALID_CHARS[max_pos];
            rev_kmer[0] = string_util::COMPLEMENT_INT[VALID_CHARS[max_pos] as usize];

            //check if we found the target
            if curr_kmer == target_kmer {
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
        return ret;
    }
    else {
        return Vec::<Vec<u8>>::new();
    }
}

pub fn assemble_from_kmer(
    bwt: &BitVectorBWT, seed_kmer: &[u8], min_count: u64, branch_limit: usize, max_branch_len: usize
) -> Vec<Vec<u8>> {
    //build up return
    let mut ret = Vec::<Vec<u8>>::new();

    //build some helper values we'll be looping through a lot
    let kmer_len = seed_kmer.len();
    let mut counts: [u64; VALID_CHARS_LEN] = [0; VALID_CHARS_LEN];
    let mut num_branched: usize = 0;
    let mut max_pos: usize;

    //the queries will populates these vectors
    let mut curr_kmer: Vec<u8> = vec![4; kmer_len];
    let mut rev_kmer: Vec<u8> = vec![4; kmer_len];

    //initialize the bridging with our seed k-mer
    let mut possible_bridges: Vec<Vec<u8>> = Vec::<Vec<u8>>::new();
    let mut seed_vec: Vec<u8> = vec![0; kmer_len];
    seed_vec.clone_from_slice(seed_kmer);
    possible_bridges.push(seed_vec);

    while possible_bridges.len() > 0 && num_branched < branch_limit {
        //get a bridge to extend
        let mut curr_bridge: Vec<u8> = possible_bridges.pop().unwrap();
        let mut curr_bridge_len = curr_bridge.len();
        num_branched += 1;

        //TODO: replace this with copy slice and call to string_util::reverse_complement_i?
        for x in {0..kmer_len} {
            curr_kmer[x] = curr_bridge[curr_bridge_len-kmer_len+x];
            rev_kmer[kmer_len-x-1] = string_util::COMPLEMENT_INT[curr_kmer[x] as usize];
        }

        while curr_bridge_len < max_branch_len {
            //shift the current k-mer over one in preparation for the last base toggle
            for x in {0..kmer_len-1} {
                curr_kmer[x] = curr_kmer[x+1];
                rev_kmer[kmer_len-x-1] = rev_kmer[kmer_len-x-2];
            }
            
            //do all the k-mer counting
            max_pos = 0;
            for x in {0..VALID_CHARS_LEN} {
                curr_kmer[kmer_len-1] = VALID_CHARS[x];
                rev_kmer[0] = string_util::COMPLEMENT_INT[VALID_CHARS[x] as usize];
                counts[x] = bwt.count_kmer(&curr_kmer) + bwt.count_kmer(&rev_kmer);
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
            for x in {0..VALID_CHARS_LEN} {
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
            
            //update k-mers
            curr_kmer[kmer_len-1] = VALID_CHARS[max_pos];
            rev_kmer[0] = string_util::COMPLEMENT_INT[VALID_CHARS[max_pos] as usize];
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
        return ret;
    }
    else {
        return Vec::<Vec<u8>>::new();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::bwt_converter::convert_to_vec;
    use crate::ropebwt2_util::create_bwt_from_strings;
    use crate::string_util::convert_stoi;
    
    fn get_constant_bwt() -> BitVectorBWT {
        //build the dataset
        let const_string = "AACGGATCAAGCTTACCAGTATTTACGT";
        let rep_count = 30;
        let mut data: Vec<&str> = vec![];
        for _i in {0..rep_count} {
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
        for _i in {0..rep_count} {
            data.push(&const_string);
        }
        //                   AACGGATCAAGCTTACCAGTATTTACGT
        //                      #         #            #
        let lesser_string = "AACTGATCAAGCTGACCAGTATTTACTT"; //3 basepair changes
        let lesser_count = 20;
        for _i in {0..lesser_count} {
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
        for _i in {0..rep_count} {
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

    #[test]
    fn test_triple_accel() {
        //this test is mostly just for my sanity of understanding the library
        let query = convert_stoi(&"AACGGATCAAGCTTACCAGTATTTACGT");
        let short_seq = convert_stoi(&"AACGGATGAA"); //one base change
        let long_seq = convert_stoi(&"AACGGATCATAGCTTACCAGTATATACGT"); //one insert, one change

        //use this one when we want to compare the full strings and get an edit distance
        let l_dist = levenshtein(&long_seq, &query);
        assert_eq!(l_dist, 2);

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
            branch_limit_factor: 4,
            branch_buffer_factor: 1.3,
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
            branch_limit_factor: 4,
            branch_buffer_factor: 1.3,
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
            branch_limit_factor: 4,
            branch_buffer_factor: 1.3,
            tail_buffer_factor: 1.00, //normally - 1.05,
            frac: 0.1,
            verbose: true
        };
        let arc_params: Arc<CorrectionParameters> = Arc::new(params);

        //the two strings in the bwt
        let const_string: String =  "AACGGATCAAGCTTACCAGTATTTACGT".to_string();
        
        let long_read_const_mod: LongReadFA = LongReadFA {
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
}