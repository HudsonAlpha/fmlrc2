
extern crate log;

use log::{info, error};
use crate::bwt_converter::*;
use crate::indexed_bit_vec::IndexedBitVec;
use crate::string_util;

pub struct BitVectorBWT {
    bwt: Vec<u8>,
    //TODO: look into ArrayVec or raw array for speed, no reason to pre-mature optimize it though
    binary_bwt: Vec<IndexedBitVec>,
    total_counts: [u64; VC_LEN],
    start_index: [u64; VC_LEN],
    end_index: [u64; VC_LEN],
    total_size: u64
}

/// Basic struct for containing a range in a BWT
pub struct BWTRange {
    /// the lower bound, inclusive
    l: u64,
    /// the upper bound, exclusive
    h: u64
}

/// Implements a BWT and FM-index based on indexed bit vectors. The implementation is quite fast compared to traditional
/// sampled FM-index, but requires more space. In general, for `N` characters with alphabet size `C`, one can expect approximately
/// 2*C*N bits of storage for both the BWT and the corresponding index. 
impl BitVectorBWT {
    /// Allocation function for the BWT, look at `load_vector(...)` for initialization.
    /// # Examples
    /// ```rust
    /// use fmlrc::bv_bwt::BitVectorBWT;
    /// let mut bwt: BitVectorBWT = BitVectorBWT::new();
    /// ```
    pub fn new() -> Self {
        Self {
            bwt: Vec::<u8>::new(),
            //binary_bwt: [IndexedBitVec::with_capacity(0); VC_LEN],
            //binary_bwt: std::iter::repeat_with(|| IndexedBitVec::with_capacity(0)).take(VC_LEN).collect::<[IndexedBitVec; VC_LEN]>>().into_inner(),
            binary_bwt: Vec::<IndexedBitVec>::new(),
            total_counts: [0; VC_LEN],
            start_index: [0; VC_LEN],
            end_index: [0; VC_LEN],
            total_size: 0
        }
    }

    /// Initializes the BWT from a compressed BWT vector.
    /// # Arguments
    /// * `bwt` - the run-length encoded BWT stored in a Vec<u8> 
    /// # Examples
    /// ```rust
    /// use std::io::Cursor;
    /// use fmlrc::bv_bwt::BitVectorBWT;
    /// use fmlrc::bwt_converter::convert_to_vec;
    /// //strings "ACGT" and "CCGG"
    /// let seq = "TG$$CAGCCG";
    /// let seq = Cursor::new(seq);
    /// let vec = convert_to_vec(seq);
    /// let mut bwt = BitVectorBWT::new();
    /// bwt.load_vector(vec);
    /// ```
    pub fn load_vector(&mut self, bwt: Vec<u8>) {
        //i am the captain now
        self.bwt = bwt;
        info!("Loading BWT from vector of length {}", self.bwt.len());
        
        //first pass does a count so we can pre-allocate the indices correctly
        self.calculate_totals();

        //now we can construct the FM-index pieces in the binary storage format for rapid speed
        self.construct_fmindex(false);
    }

    /// This calculates the total number of each symbol in the BWT from the compressed representation.
    /// Running this is part of initialization and is a pre-requisite for FM-index construction.
    fn calculate_totals(&mut self) {
        let mut prev_char: u8 = 255;
        let mut current_char: u8;
        let mut power_multiple: u64 = 1;
        let mut current_count: u64;

        //go through each compressed block in the RLE encoded vector to calculate total character counts
        self.total_counts = [0; VC_LEN];
        for value in &self.bwt {
            current_char = value & MASK;
            if current_char == prev_char {
                power_multiple *= NUM_POWER as u64;
            }
            else {
                power_multiple = 1;
            }
            prev_char = current_char;
            current_count = (value >> LETTER_BITS) as u64 * power_multiple;
            self.total_counts[current_char as usize] += current_count;
        }

        //calculate start/end indices from the total
        self.start_index = [0; VC_LEN];
        self.end_index = [0; VC_LEN];
        let mut sum_offset: u64 = 0;
        for i in {0..VC_LEN} {
            self.start_index[i] = sum_offset;
            sum_offset += self.total_counts[i];
            self.end_index[i] = sum_offset;
        }
        self.total_size = self.end_index[VC_LEN-1];
        info!("Loaded BWT with symbol counts: {:?}", self.total_counts);
    }

    /// This function will build the indices that represent the FM-index.
    /// Should be called once as a part of initialization.
    /// # Arguments
    /// * `store_dollar` - if True, the index for the sentinel will also be constructed
    fn construct_fmindex(&mut self, store_dollar: bool) {
        //setup vectors
        self.binary_bwt = Vec::<IndexedBitVec>::with_capacity(VC_LEN);
        for x in {0..VC_LEN} {
            if x != 0 || store_dollar {
                self.binary_bwt.push(IndexedBitVec::with_capacity(self.total_size as usize));
            }
            else {
                self.binary_bwt.push(IndexedBitVec::with_capacity(1));
            }
        }
        
        //set up parse variables
        let mut prev_char: u8 = 0;
        let mut total_char_count: u64 = 0;
        let mut power_multiple: u64 = 1;
        let mut bwt_index: u64 = 0;
        let mut current_char: u8;

        //now we loop through the compressed BWT
        let num_bytes: usize = self.bwt.len();
        for x in {0..num_bytes} {
            current_char = self.bwt[x] & MASK;
            if current_char == prev_char {
                //same char, increment lengths
                total_char_count += (self.bwt[x] >> LETTER_BITS) as u64 * power_multiple;
                power_multiple *= NUM_POWER as u64;
            }
            else {
                //first save the current FM-index entry
                if prev_char != 0 || store_dollar {
                    for y in {bwt_index..bwt_index+total_char_count} {
                        unsafe {
                            self.binary_bwt[prev_char as usize].set_bit(y as usize);
                        }
                    }
                }
                
                //now add the previous
                bwt_index += total_char_count;
                prev_char = current_char;
                total_char_count = (self.bwt[x] >> LETTER_BITS) as u64;
                power_multiple = NUM_POWER as u64;
            }
        }
        
        //do the last block of characters
        if prev_char != 0 || store_dollar {
            for y in {bwt_index..bwt_index+total_char_count} {
                unsafe {
                    self.binary_bwt[prev_char as usize].set_bit(y as usize);
                }
            }
        }
        
        //build the indices for each vector
        for x in {0..VC_LEN} {
            if x != 0 || store_dollar {
                self.binary_bwt[x].build_index(self.start_index[x]);
            }
        }
    }

    /// Returns the total number of occurences of a given symbol
    /// # Arguments
    /// * `symbol` - the symbol in integer form
    /// # Examples
    /// ```rust
    /// # use std::io::Cursor;
    /// # use fmlrc::bv_bwt::BitVectorBWT;
    /// # use fmlrc::bwt_converter::convert_to_vec;
    /// # let seq = "TG$$CAGCCG";
    /// # let seq = Cursor::new(seq);
    /// # let vec = convert_to_vec(seq);
    /// # let mut bwt = BitVectorBWT::new();
    /// # bwt.load_vector(vec);
    /// let string_count = bwt.get_total_counts(0);
    /// assert_eq!(string_count, 2);
    /// ```
    #[inline]
    pub fn get_total_counts(&self, symbol: u8) -> u64 {
        self.total_counts[symbol as usize]
    }

    /// Performs a range constraint on a BWT range. This implicitly represents prepending a character `sym` to a k-mer
    /// represented by `input_range` to create a new range representing a (k+1)-mer.
    /// # Arguments
    /// * `sym` - the symbol to pre-pend in integer form
    /// * `input_range` - the range to pre-pend to
    #[inline]
    fn constrain_range(&self, sym: u8, input_range: BWTRange) -> BWTRange {
        let l = self.binary_bwt[sym as usize].rank(input_range.l as usize);
        let h = self.binary_bwt[sym as usize].rank(input_range.h as usize);
        BWTRange {
            l: l,
            h: h
        }
    }

    /// Returns the total number of occurrences of a given k-mer in the BWT.
    /// # Arguments
    /// * `kmer` - the integer-encoded kmer sequence to count
    /// # Examples
    /// ```rust
    /// # use std::io::Cursor;
    /// # use fmlrc::bv_bwt::BitVectorBWT;
    /// # use fmlrc::bwt_converter::convert_to_vec;
    /// # let seq = "TG$$CAGCCG";
    /// # let seq = Cursor::new(seq);
    /// # let vec = convert_to_vec(seq);
    /// # let mut bwt = BitVectorBWT::new();
    /// # bwt.load_vector(vec);
    /// let kmer_count = bwt.count_kmer(&vec![1, 2, 3, 5]); //ACGT
    /// assert_eq!(kmer_count, 1);
    /// ```
    #[inline]
    pub fn count_kmer(&self, kmer: &[u8]) -> u64 {
        //init to everything
        let mut ret: BWTRange = BWTRange {
            l: 0,
            h: self.total_size
        };
        
        //TODO: test if this is the fastest way to do this loop (with internal break & such)
        for i in {0..kmer.len()}.rev() {
            ret = self.constrain_range(kmer[i], ret);
            if ret.h == ret.l {
                return 0;
            }
        }

        //return the delta
        ret.h-ret.l
    }

    /// Returns the counts for all k-mers in a given sequence as a Vec<u64>. This function looks at both forward and
    /// reverse complemented counts. The total number of counts will be `max(seq.len() - kmer_size+1, 0)`.
    /// # Arguments
    /// * `seq` - a sequence in integer form
    /// * `kmer_size` - the size or `k` parameter for the k-mers
    /// ```rust
    /// # use std::io::Cursor;
    /// # use fmlrc::bv_bwt::BitVectorBWT;
    /// # use fmlrc::bwt_converter::convert_to_vec;
    /// # let seq = "TG$$CAGCCG";
    /// # let seq = Cursor::new(seq);
    /// # let vec = convert_to_vec(seq);
    /// # let mut bwt = BitVectorBWT::new();
    /// # bwt.load_vector(vec);
    /// let pu = bwt.count_pileup(&vec![2, 3, 3], 2); //get all 2-mer counts for "CGG"
    /// assert_eq!(pu, vec![4, 2]); //counts include reverse-complement
    /// ```
    #[inline]
    pub fn count_pileup(&self, seq: &Vec<u8>, kmer_size: usize) -> Vec<u64> {
        let seq_len = seq.len();
        if seq_len < kmer_size {
            return vec![0];
        }

        //get reverse complement
        let rev_seq = string_util::reverse_complement_i(&seq);
        
        //build up return and count the k-mers
        let num_counts = seq_len - kmer_size + 1;
        let mut ret = Vec::<u64>::with_capacity(num_counts);
        for x in {0..num_counts} {
            ret.push(
                self.count_kmer(&seq[x..x+kmer_size])+
                self.count_kmer(&rev_seq[seq_len-kmer_size-x..seq_len-x])
            );
        }
        ret
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_simple_index() {
        //strings "ACGT\nCCGG"
        let seq = "TG$$CAGCCG";
        let seq = Cursor::new(seq);
        let vec = convert_to_vec(seq);
        
        let mut bwt = BitVectorBWT::new();
        bwt.load_vector(vec);
        
        let expected_totals = vec![2, 1, 3, 3, 0, 1];
        let expected_starts = vec![0, 2, 3, 6, 9, 9];
        let expected_ends = vec![2, 3, 6, 9, 9, 10];

        for i in {0..6} {
            //make sure the total counts are correct
            assert_eq!(bwt.get_total_counts(i as u8), expected_totals[i]);
            assert_eq!(bwt.start_index[i], expected_starts[i]);
            assert_eq!(bwt.end_index[i], expected_ends[i]);

            //basic count of all 1-mers characters (excluding '$' of course)
            if i != 0 {
                assert_eq!(bwt.count_kmer(&vec![i as u8]), expected_totals[i]);
            }
        }

        //original strings "ACGT\nCCGG"
        //2-mers
        assert_eq!(bwt.count_kmer(&vec![1, 2]), 1); //AC
        assert_eq!(bwt.count_kmer(&vec![2, 3]), 2); //CG
        assert_eq!(bwt.count_kmer(&vec![3, 5]), 1); //GT
        assert_eq!(bwt.count_kmer(&vec![2, 2]), 1); //CC
        assert_eq!(bwt.count_kmer(&vec![3, 3]), 1); //GG
        assert_eq!(bwt.count_kmer(&vec![5, 3]), 0); //TG - not present

        //3-mers
        assert_eq!(bwt.count_kmer(&vec![1, 2, 3]), 1); //ACG
        assert_eq!(bwt.count_kmer(&vec![2, 3, 5]), 1); //CGT
        assert_eq!(bwt.count_kmer(&vec![2, 2, 3]), 1); //CCG
        assert_eq!(bwt.count_kmer(&vec![2, 3, 3]), 1); //CGG
        assert_eq!(bwt.count_kmer(&vec![2, 3, 1]), 0); //CGA - not present

        //4-mers
        assert_eq!(bwt.count_kmer(&vec![1, 2, 3, 5]), 1); //ACGT
        assert_eq!(bwt.count_kmer(&vec![2, 2, 3, 3]), 1); //CCGG
        assert_eq!(bwt.count_kmer(&vec![5, 5, 5, 5]), 0); //TTTT - not present
    }

    #[test]
    fn test_simple_pileup() {
        //strings "ACGT\nCCGG"
        let seq = "TG$$CAGCCG";
        let seq = Cursor::new(seq);
        let vec = convert_to_vec(seq);
        
        let mut bwt = BitVectorBWT::new();
        bwt.load_vector(vec);

        let pu = bwt.count_pileup(&vec![1, 2, 3, 5], 1);
        assert_eq!(pu, vec![2, 6, 6, 2]);
        
        let pu = bwt.count_pileup(&vec![2, 3, 3], 2);
        assert_eq!(pu, vec![4, 2]);

        let pu = bwt.count_pileup(&vec![1, 2, 3, 5], 3); //this string is it's own rev-comp, so we get two counts on the 3-mers
        assert_eq!(pu, vec![2, 2]);
    }
}