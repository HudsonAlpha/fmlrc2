
extern crate log;

use log::info;
use std::io::prelude::*;
use std::fs;

use crate::bwt_converter::*;
use crate::indexed_bit_vec::IndexedBitVec;
use crate::string_util;

const ZERO_COUNT_VEC: [u64; 256] = [0; 256];

// sweet spot seems to be 8, providing a definitely speed boost while using ~25MB,
// so negligible overhead compared to FMLRC as a whole
// 9 is ~153MB, and 10 is ~922MB; 9 does offer a boost, but I'd be hesistant to use >=10 due to mem reqs
const DEFAULT_CACHE_K: usize = 8; 

/// Structure that contains a bit vector-based BWT+FM-index implementation
pub struct BitVectorBWT {
    bwt: Vec<u8>,
    //TODO: look into ArrayVec or raw array for speed, no reason to pre-mature optimize it though
    binary_bwt: Vec<IndexedBitVec>,
    total_counts: [u64; VC_LEN],
    start_index: [u64; VC_LEN],
    end_index: [u64; VC_LEN],
    total_size: u64,
    fixed_init: [BWTRange; 4],
    cache_k: usize,
    kmer_cache: Vec<BWTRange>
}

/// Basic struct for containing a range in a BWT.
/// Only contains fields `l` and `h`, representing a range [l, h).
#[derive(Clone,Copy,Default,Debug,Eq,PartialEq)]
pub struct BWTRange {
    /// the lower bound, inclusive
    pub l: u64,
    /// the upper bound, exclusive
    pub h: u64
}

impl Default for BitVectorBWT {
    /// Default function
    /// # Examples
    /// ```rust
    /// use fmlrc::bv_bwt::BitVectorBWT;
    /// let mut bwt: BitVectorBWT = Default::default();
    /// ```    
    fn default() -> Self {
        Self {
            bwt: Vec::<u8>::new(),
            binary_bwt: Vec::<IndexedBitVec>::new(),
            total_counts: [0; VC_LEN],
            start_index: [0; VC_LEN],
            end_index: [0; VC_LEN],
            total_size: 0,
            fixed_init: [Default::default(); 4],
            cache_k: DEFAULT_CACHE_K,
            kmer_cache: Vec::<BWTRange>::new()
        }
    }
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
        Default::default()
    }

    /// Allocation function for the BWT with a setting for cache k-mer length.
    /// If using `new()`, the default cache uses `k=8` which leads to ~25MB of cache size.
    /// Each integer increment of `k` increases the cache size 6-fold, but also generally improves runtime.
    /// If memory is not a major issue, using `k=10` can provide a noticeable speedup while maintaining a cache size <1GB.
    /// Due to the minimal size of the 8-mer cache, we do not recommend using <8 as that will likely only cause a performance hit.
    /// # Arguments
    /// * `cache_k` - sets the k-mer size used in the cache (default is 8)
    /// # Examples
    /// ```rust
    /// // creates a BWT with ~1GB of memory used for caching
    /// use fmlrc::bv_bwt::BitVectorBWT;
    /// let mut bwt: BitVectorBWT = BitVectorBWT::with_cache_size(10);
    /// ```
    pub fn with_cache_size(cache_k: usize) -> Self {
        let mut ret: BitVectorBWT = Default::default();
        ret.cache_k = cache_k;
        ret
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
        
        //we copied it in, standard init now
        self.standard_init();
    }

    /// Initializes the BWT from the numpy file format for compressed BWTs
    /// # Arguments
    /// * `filename` - the name of the file to load into memory
    /// # Examples
    /// ```ignore
    /// use fmlrc::bv_bwt::BitVectorBWT;
    /// let mut bwt = BitVectorBWT::new();
    /// let filename: String = "/path/to/my/file/comp_msbwt.npy".to_string();
    /// bwt.load_numpy_file(&filename);
    /// ```
    pub fn load_numpy_file(&mut self, filename: &str) -> std::io::Result<()> {
        //read the numpy header: http://docs.scipy.org/doc/numpy-1.10.1/neps/npy-format.html
        //get the initial file size
        let file_metadata: fs::Metadata = fs::metadata(&filename)?;
        let full_file_size: u64 = file_metadata.len();

        //read the initial fixed header
        let mut file = fs::File::open(&filename)?;
        let mut init_header: Vec<u8> = vec![0; 16];
        let read_count: usize = file.read(&mut init_header[..])?;
        if read_count != 16 {
            panic!("Could not read initial 16 bytes of header for file {:?}", filename);
        }

        //read the dynamic header
        let header_len: usize = init_header[8] as usize + 256 * init_header[9] as usize;
        let mut skip_bytes: usize = 10+header_len;
        if skip_bytes % 16 != 0 {
            skip_bytes = ((skip_bytes / 16)+1)*16;
        }
        let mut skip_header: Vec<u8> = vec![0; skip_bytes-16];
        let read_count: usize = file.read(&mut skip_header[..])?;
        if read_count != skip_bytes-16 {
            panic!("Could not read bytes 16-{:?} of header for file {:?}", skip_bytes, filename);
        }

        //finally read in everything else
        let bwt_disk_size: u64 = full_file_size - skip_bytes as u64;
        self.bwt = vec![0; bwt_disk_size as usize];
        let read_count: usize = file.read_to_end(&mut self.bwt)?;
        if read_count as u64 != bwt_disk_size {
            panic!("Only read {:?} of {:?} bytes of BWT body for file {:?}", read_count, bwt_disk_size, filename);
        }

        //TODO: I imagine we want to use the info here somehow?
        //printf("loaded bwt with %lu compressed values\n", this->bwt.size());
        info!("Loading BWT with {:?} compressed values", bwt_disk_size);

        //we loaded the file into memory, standard init now
        self.standard_init();

        Ok(())
    }

    fn standard_init(&mut self) {
        //first pass does a count so we can pre-allocate the indices correctly
        self.calculate_totals();

        //now we can construct the FM-index pieces in the binary storage format for rapid speed
        self.construct_fmindex(false);

        //now do the fixed initialization
        let full_range: BWTRange = BWTRange {
            l: 0,
            h: self.total_size
        };
        unsafe {
            self.fixed_init[0] = self.constrain_range(1, &full_range);
            self.fixed_init[1] = self.constrain_range(2, &full_range);
            self.fixed_init[2] = self.constrain_range(3, &full_range);
            self.fixed_init[3] = self.constrain_range(5, &full_range);
        }

        self.populate_cache(false);

        info!("Finished BWT initialization.");
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
        for i in 0..VC_LEN {
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
        info!("Allocating binary vectors...");
        self.binary_bwt = Vec::<IndexedBitVec>::with_capacity(VC_LEN);
        for x in 0..VC_LEN {
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
        info!("Calculating binary vectors...");
        let num_bytes: usize = self.bwt.len();
        for x in 0..num_bytes {
            current_char = self.bwt[x] & MASK;
            if current_char == prev_char {
                //same char, increment lengths
                total_char_count += (self.bwt[x] >> LETTER_BITS) as u64 * power_multiple;
                power_multiple *= NUM_POWER as u64;
            }
            else {
                //first save the current FM-index entry
                if prev_char != 0 || store_dollar {
                    for y in bwt_index..bwt_index+total_char_count {
                        self.binary_bwt[prev_char as usize].set_bit(y as usize);
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
            for y in bwt_index..bwt_index+total_char_count {
                self.binary_bwt[prev_char as usize].set_bit(y as usize);
            }
        }
        
        //build the indices for each vector
        info!("Constructing FM-indices...");
        for x in 0..VC_LEN {
            if x != 0 || store_dollar {
                self.binary_bwt[x].build_index(self.start_index[x]);
            }
        }
    }

    /// This will pull the full cache index for a single in-order k-mer (i.e. forward order)
    #[inline]
    fn get_cache_index(&self, arr: &[u8]) -> usize {
        assert!(arr.len() >= self.cache_k);
        let mut ret: usize = 0;
        for i in 0..self.cache_k {
            ret *= VC_LEN;
            ret += arr[i] as usize;
        }
        ret
    }

    /// This will pull the full cache index for a single reversed k-mer
    #[inline]
    fn get_rev_cache_index(&self, arr: &[u8]) -> usize {
        assert!(arr.len() >= self.cache_k);
        let mut ret: usize = 0;
        for i in (0..self.cache_k).rev() {
            ret *= VC_LEN;
            ret += arr[i] as usize;
        }
        ret
    }

    /// This will pull all fixed extensions `c` for (k-1)-mer K, getting 4 cache keys `Kc`
    #[inline]
    fn get_fixed_cache_index(&self, arr: &[u8]) -> [usize; 4] {
        assert!(arr.len() >= self.cache_k-1);
        let mut initial_value: usize = 0;
        for i in 0..self.cache_k-1 {
            initial_value *= VC_LEN;
            initial_value += arr[i] as usize;
        }
        let mut ret: [usize; 4] = [initial_value; 4];
        for i in 0..4 {
            ret[i] *= VC_LEN;
        }
        ret[0] += 1;
        ret[1] += 2;
        ret[2] += 3;
        ret[3] += 5;
        ret
    }

    /// This will pre-build a cache of all k-mers of length `cache_k`.
    /// If `store_dollar` is set, it will build all k-mers with the `$` character also.
    /// Only needs to be called once after indexing is complete.
    /// # Arguments
    /// * `store_dollar` - if true, create cache entries for the k-mers with `$` characters
    fn populate_cache(&mut self, store_dollar: bool) {
        info!("Building {:?}-mer cache...", self.cache_k);
        let initial_value: u8 = if store_dollar {0} else {1};
        //let mut current_key: [u8; self.cache_k] = [initial_value; self.cache_k];
        let mut current_key: Vec<u8> = vec![initial_value; self.cache_k];
        let mut cache_index: usize;
        let mut ck_ind: usize;
        //let final_key: [u8; self.cache_k] = [5; self.cache_k];
        let final_key: Vec<u8> = vec![5; self.cache_k];
        self.kmer_cache = vec![Default::default(); VC_LEN.pow(self.cache_k as u32)];
        loop {
            //set the cache
            cache_index = self.get_cache_index(&current_key);
            self.kmer_cache[cache_index] = self.range_kmer(&current_key);
            
            //check if we're at the final key
            if current_key == final_key {
                break;
            }

            //increment our current_key
            ck_ind = self.cache_k - 1;
            current_key[ck_ind] += 1;
            while current_key[ck_ind] == 6 {
                current_key[ck_ind-1] += 1;
                current_key[ck_ind] = initial_value;
                ck_ind -= 1;
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
    unsafe fn constrain_range(&self, sym: u8, input_range: &BWTRange) -> BWTRange {
        BWTRange {
            l: self.binary_bwt.get_unchecked(sym as usize).rank_unchecked(input_range.l as usize),
            h: self.binary_bwt.get_unchecked(sym as usize).rank_unchecked(input_range.h as usize)
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
        let mut ret: BWTRange;
        let cut_kmer: &[u8];
        
        //check for cache entry
        if kmer.len() >= self.cache_k {
            ret = self.kmer_cache[self.get_cache_index(&kmer[kmer.len()-self.cache_k..])];
            cut_kmer = &kmer[..kmer.len()-self.cache_k];
        } else {
            cut_kmer = kmer;
            ret = BWTRange {
                l: 0,
                h: self.total_size
            };
        }
        
        //go through what remains in reverse
        for c in cut_kmer.iter().rev() {
            assert!(*c < VC_LEN as u8);
            unsafe {
                ret = self.constrain_range(*c, &ret);
            }
            if ret.h == ret.l {
                return 0;
            }
        }

        //return the delta
        ret.h-ret.l
    }

    /// This function is only used to initially build up the k-mer cache
    #[inline]
    fn range_kmer(&self, kmer: &[u8]) -> BWTRange {
        //init to everything
        let mut ret: BWTRange = BWTRange {
            l: 0,
            h: self.total_size
        };
        
        //must do full k-mer since cache doesn't exist yet
        for c in kmer.iter().rev() {
            assert!(*c < VC_LEN as u8);
            unsafe {
                ret = self.constrain_range(*c, &ret);
            }
        }
        ret
    }

    /// Returns the total number of occurrences of each given symbol before a given k-mer in the BWT.
    /// Functionally, for a k-mer `K` and a symbol set `C`, this is identical to calculating the occurence of `cK`
    /// for each `c` in `C`. This function reduces the work by re-using the shared components of the calculation.
    /// # Arguments
    /// * `kmer` - the integer-encoded kmer sequence to count
    /// * `symbols` - the integer-encoded symbols to pre-pend to the k-mer and count.
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
    /// let kmer_counts = bwt.prefix_kmer(&vec![2, 3], &vec![1, 2, 3, 5]); //count XCG
    /// assert_eq!(kmer_counts, vec![1, 1, 0, 0]); //the set has one occurrence each of ACG and CCG
    /// ```
    #[inline]
    pub fn prefix_kmer(&self, kmer: &[u8], symbols: &[u8]) -> Vec<u64> {
        let mut ret_counts: Vec<u64> = vec![0; symbols.len()];
        self.prefix_kmer_noalloc(kmer, symbols, &mut ret_counts[..]);
        ret_counts
    }

    /// Populates the total number of occurrences of each given symbol before a given k-mer in the BWT.
    /// Functionally, for a k-mer `K` and a symbol set `C`, this is identical to calculating the occurence of `cK`
    /// for each `c` in `C`. This function reduces the work by re-using the shared components of the calculation.
    /// This function does not allocate a count array, but populates the passed in value instead.
    /// # Arguments
    /// * `kmer` - the integer-encoded kmer sequence to count
    /// * `symbols` - the integer-encoded symbols to pre-pend to the k-mer and count.
    /// * `counts` - a mutable array to be populated with the outputs, must be at least the same length as `symbols`
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
    /// let mut kmer_counts = vec![0; 4];
    /// bwt.prefix_kmer_noalloc(&vec![2, 3], &vec![1, 2, 3, 5], &mut kmer_counts[..]); //count XCG
    /// assert_eq!(kmer_counts, vec![1, 1, 0, 0]); //the set has one occurrence each of ACG and CCG
    /// ```
    #[inline]
    pub fn prefix_kmer_noalloc(&self, kmer: &[u8], symbols: &[u8], counts: &mut [u64]) {
        //init to everything
        let mut ret: BWTRange;
        let cut_kmer: &[u8];
        
        //check for cache results
        if kmer.len() >= self.cache_k {
            ret = self.kmer_cache[self.get_cache_index(&kmer[kmer.len()-self.cache_k..])];
            cut_kmer = &kmer[..kmer.len()-self.cache_k];
        } else {
            ret = BWTRange {
                l: 0,
                h: self.total_size
            };
            cut_kmer = kmer;
        }
        
        //go through remaining sequence in reverse order
        for c in cut_kmer.iter().rev() {
            assert!(*c < VC_LEN as u8);
            unsafe {
                ret = self.constrain_range(*c, &ret);
            }
            if ret.h == ret.l {
                counts[0..symbols.len()].clone_from_slice(&ZERO_COUNT_VEC[0..symbols.len()]);
                return;
            }
        }

        //pre-pend the passed symbols as final counts
        for (i, c) in symbols.iter().enumerate() {
            assert!(*c < VC_LEN as u8);
            let subrange = unsafe { self.constrain_range(*c, &ret) };
            counts[i] = subrange.h-subrange.l
        }
    }

    /// This is a specialty function for fmlrc that accepts a *reversed* k-mer sequence along with a mutable
    /// counts array. It will then calculate the occurences of the forward k-mer if pre-pended with the fixed
    /// characters T, G, C, A.  This implicitly represents doing the reverse-complement counts of the (k+1)-mer sequence.
    /// Given a k-mer `K`, the results counts array will contain the number of occurences of: `[T+rev(K), G+rev(K), C+rev(K), A+rev(K)]`.
    /// In fmlrc, this is then added to forward counts for the rev-comp sequences to obtain total counts.
    /// This function is counter-intuitive, but efficient; make sure you understand it before use.
    /// # Arguments
    /// * `rev_kmer` - a k-mer sequence that will be traversed in the forward direction for counting (normal k-mer counting is from end to start)
    /// * `counts` - a mutable array that will be populated with counts for (k+1)-mers `[T+rev(K), G+rev(K), C+rev(K), A+rev(K)]`
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
    /// let mut kmer_counts = vec![0; 4];
    /// bwt.prefix_revkmer_noalloc_fixed(&vec![3, 2], &mut kmer_counts[..]); //count XCG (note the k-mer is revesed order GC=[3, 2])
    /// assert_eq!(kmer_counts, vec![0, 0, 1, 1]); //the set has one occurrence each of ACG and CCG; so complemented is [0, 0, 1, 1]
    /// ```
    #[inline]
    pub fn prefix_revkmer_noalloc_fixed(&self, rev_kmer: &[u8], counts: &mut [u64]) {
        //init to everything
        assert!(counts.len() >= 4);
        let mut ret: BWTRange;
        let cut_kmer: &[u8];

        //check cache
        if rev_kmer.len() >= self.cache_k {
            let cache_index: usize = self.get_rev_cache_index(&rev_kmer);
            ret = self.kmer_cache[cache_index];
            cut_kmer = &rev_kmer[self.cache_k..];
        } else {
            ret = BWTRange {
                l: 0,
                h: self.total_size
            };
            cut_kmer = rev_kmer;
        }

        //iterate forward
        for c in cut_kmer.iter() {
            assert!(*c < VC_LEN as u8);
            unsafe {
                ret = self.constrain_range(*c, &ret);
            }
            if ret.h == ret.l {
                //impossible to have any counts now, zero them out
                counts[0] = 0;
                counts[1] = 0;
                counts[2] = 0;
                counts[3] = 0;
                return;
            }
        }
        //reverse complemented, so T G C A
        let subrange = unsafe { self.constrain_range(5, &ret) };
        counts[0] = subrange.h-subrange.l;
        let subrange = unsafe { self.constrain_range(3, &ret) };
        counts[1] = subrange.h-subrange.l;
        let subrange = unsafe { self.constrain_range(2, &ret) };
        counts[2] = subrange.h-subrange.l;
        let subrange = unsafe { self.constrain_range(1, &ret) };
        counts[3] = subrange.h-subrange.l;
    }

    /// This is a specialty function for fmlrc that accepts a k-mer sequence along with a mutable counts array. 
    /// It will then calculate the counts for that k-mer sequence with each postfix [A, C, G, T].
    /// Given a k-mer `K`, the results counts array will contain the number of occurences of (k+1)-mers: `[KA, KC, KG, KT]`.
    /// In fmlrc, this is then added to reverse-complemented counts to obtain total counts.
    /// # Arguments
    /// * `kmer` - a k-mer sequence, `K`
    /// * `counts` - a mutable array that will be populated with counts for (k+1)-mers `[KA, KC, KG, KT]`
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
    /// let mut kmer_counts = vec![0; 4];
    /// bwt.postfix_kmer_noalloc_fixed(&vec![2, 3], &mut kmer_counts[..]); //count CGX
    /// assert_eq!(kmer_counts, vec![0, 0, 1, 1]); //the set has one occurrence each of CGG and CGT
    /// ```
    #[inline]
    pub fn postfix_kmer_noalloc_fixed(&self, kmer: &[u8], counts: &mut [u64]) {
        //init to everything
        assert!(counts.len() >= 4);
        let mut ranges: [BWTRange; 4];
        let cut_kmer: &[u8];
        if kmer.len() >= self.cache_k {
            let cache_indices: [usize; 4] = self.get_fixed_cache_index(&kmer[kmer.len()-self.cache_k+1..]);
            ranges = [
                self.kmer_cache[cache_indices[0]],
                self.kmer_cache[cache_indices[1]],
                self.kmer_cache[cache_indices[2]],
                self.kmer_cache[cache_indices[3]]
            ];
            cut_kmer = &kmer[0..kmer.len()-self.cache_k+1];
        } else {
            ranges = self.fixed_init.clone();
            cut_kmer = kmer;
        }

        //go through whatever remains in reverse; single iterator with checks is fastest
        for c in cut_kmer.iter().rev() {
            assert!(*c < VC_LEN as u8);
            unsafe {
                for x in 0..4 {
                    if ranges[x].l != ranges[x].h {
                        ranges[x] = self.constrain_range(*c, &ranges[x]);
                    }
                }
            }
        }
        
        //store the deltas
        for x in 0..4 {
            counts[x] = ranges[x].h - ranges[x].l;
        }
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
    pub fn count_pileup(&self, seq: &[u8], kmer_size: usize) -> Vec<u64> {
        let seq_len = seq.len();
        if seq_len < kmer_size {
            return vec![0];
        }

        //get reverse complement
        let rev_seq = string_util::reverse_complement_i(&seq);
        
        //build up return and count the k-mers
        let num_counts = seq_len - kmer_size + 1;
        let mut ret = Vec::<u64>::with_capacity(num_counts);
        for x in 0..num_counts {
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
    use crate::ropebwt2_util::stream_bwt_from_fastqs;
    use flate2::{Compression, GzBuilder};
    use std::io::Cursor;
    use tempfile::{Builder, NamedTempFile};

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

        for i in 0..6 {
            //make sure the total counts are correct
            assert_eq!(bwt.get_total_counts(i as u8), expected_totals[i]);
            assert_eq!(bwt.start_index[i], expected_starts[i]);
            assert_eq!(bwt.end_index[i], expected_ends[i]);

            //basic count of all 1-mers characters (excluding '$' of course)
            if i != 0 {
                assert_eq!(bwt.count_kmer(&vec![i as u8]), expected_totals[i]);
            }
        }

        //make sure the fixed inits are correct
        assert_eq!(bwt.fixed_init[0], BWTRange {l:expected_starts[1], h:expected_ends[1]});
        assert_eq!(bwt.fixed_init[1], BWTRange {l:expected_starts[2], h:expected_ends[2]});
        assert_eq!(bwt.fixed_init[2], BWTRange {l:expected_starts[3], h:expected_ends[3]});
        assert_eq!(bwt.fixed_init[3], BWTRange {l:expected_starts[5], h:expected_ends[5]});

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

        //prefix counters
        assert_eq!(bwt.prefix_kmer(&vec![2], &vec![1, 2, 3, 5]), vec![1, 1, 0, 0]); //XC
        assert_eq!(bwt.prefix_kmer(&vec![2, 3], &vec![1, 2, 3, 5]), vec![1, 1, 0, 0]); //XCG
        assert_eq!(bwt.prefix_kmer(&vec![5, 5, 5, 5], &vec![1, 2, 3, 5]), vec![0, 0, 0, 0]); //XTTTT - absent

        //prefix counts
        let mut counts: Vec<u64> = vec![0, 0, 0, 0];
        bwt.prefix_revkmer_noalloc_fixed(&vec![2], &mut counts);
        assert_eq!(counts, vec![0, 0, 1, 1]);//XC (but rev-comped)
        bwt.prefix_revkmer_noalloc_fixed(&vec![3, 2], &mut counts);
        assert_eq!(counts, vec![0, 0, 1, 1]);//XCG (but rev-comped)
        bwt.prefix_revkmer_noalloc_fixed(&vec![5, 5, 5, 5], &mut counts);
        assert_eq!(counts, vec![0, 0, 0, 0]);//XTTTT
        
        //postfix counts
        bwt.postfix_kmer_noalloc_fixed(&vec![2], &mut counts);
        assert_eq!(counts, vec![0, 1, 2, 0]); //CC and CG
        bwt.postfix_kmer_noalloc_fixed(&vec![2, 3], &mut counts);
        assert_eq!(counts, vec![0, 0, 1, 1]); //CGG and CGT
        bwt.postfix_kmer_noalloc_fixed(&vec![5, 5, 5, 5], &mut counts);
        assert_eq!(counts, vec![0, 0, 0, 0]); //none
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

    fn write_strings_to_fqgz(data: Vec<&str>) -> NamedTempFile {
        let file: NamedTempFile = Builder::new().prefix("temp_data_").suffix(".fq.gz").tempfile().unwrap();
        let mut gz = GzBuilder::new().write(file, Compression::default());
        let mut i: usize = 0;
        for s in data {
            writeln!(gz, "@seq_{}\n{}\n+\n{}", i, s, "F".repeat(s.len())).unwrap();
            i += 1;
        }

        //have to keep the file handle or everything blows up
        gz.finish().unwrap()
    }

    #[test]
    fn test_load_from_file() {
        //strings - "CCGT\nACG\nN"
        //build the BWT
        let data: Vec<&str> = vec!["CCGT", "N", "ACG"];
        let file = write_strings_to_fqgz(data);
        let fastq_filenames: Vec<&str> = vec![
            &file.path().to_str().unwrap(),
        ];

        //stream and compress the BWT
        let bwt_stream = stream_bwt_from_fastqs(&fastq_filenames).unwrap();
        let compressed_bwt = convert_to_vec(bwt_stream);
        
        //save the output to a temporary numpy file
        let bwt_file: NamedTempFile = Builder::new().prefix("temp_data_").suffix(".npy").tempfile().unwrap();
        let filename: String = bwt_file.path().to_str().unwrap().to_string();
        save_bwt_numpy(&compressed_bwt[..], &filename).unwrap();
        
        //load it back in and verify counts
        let mut bwt = BitVectorBWT::new();
        bwt.load_numpy_file(&filename).unwrap();

        let expected_totals = vec![3, 1, 3, 2, 1, 1];
        for i in 0..6 {
            //make sure the total counts are correct
            assert_eq!(bwt.get_total_counts(i as u8), expected_totals[i]);
        }
    }
}