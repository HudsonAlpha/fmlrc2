
use math::round;

/*
TODO: when rust allows loops in const functions, this will be preferred IMO
const fn build_low_shift() -> [u64; 65] {
    let mut ret: [u64; 65] = [0; 65];
    for set_count in 1..64 {
        ret[set_count] = 0xFFFF_FFFF_FFFF_FFFF >> (64 - set_count);
    }
    ret
}
const LOW_SET_FLAGS: [u64; 65] = build_low_shift();
*/

//interim hard set the low-set bits array, replace with above once fixed
const LOW_SET_FLAGS: [u64; 65] = [
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0001, 0x0000_0000_0000_0003, 0x0000_0000_0000_0007, 0x0000_0000_0000_000F,
    0x0000_0000_0000_001F, 0x0000_0000_0000_003F, 0x0000_0000_0000_007F, 0x0000_0000_0000_00FF,
    0x0000_0000_0000_01FF, 0x0000_0000_0000_03FF, 0x0000_0000_0000_07FF, 0x0000_0000_0000_0FFF,
    0x0000_0000_0000_1FFF, 0x0000_0000_0000_3FFF, 0x0000_0000_0000_7FFF, 0x0000_0000_0000_FFFF,
    0x0000_0000_0001_FFFF, 0x0000_0000_0003_FFFF, 0x0000_0000_0007_FFFF, 0x0000_0000_000F_FFFF,
    0x0000_0000_001F_FFFF, 0x0000_0000_003F_FFFF, 0x0000_0000_007F_FFFF, 0x0000_0000_00FF_FFFF,
    0x0000_0000_01FF_FFFF, 0x0000_0000_03FF_FFFF, 0x0000_0000_07FF_FFFF, 0x0000_0000_0FFF_FFFF,
    0x0000_0000_1FFF_FFFF, 0x0000_0000_3FFF_FFFF, 0x0000_0000_7FFF_FFFF, 0x0000_0000_FFFF_FFFF,
    0x0000_0001_FFFF_FFFF, 0x0000_0003_FFFF_FFFF, 0x0000_0007_FFFF_FFFF, 0x0000_000F_FFFF_FFFF,
    0x0000_001F_FFFF_FFFF, 0x0000_003F_FFFF_FFFF, 0x0000_007F_FFFF_FFFF, 0x0000_00FF_FFFF_FFFF,
    0x0000_01FF_FFFF_FFFF, 0x0000_03FF_FFFF_FFFF, 0x0000_07FF_FFFF_FFFF, 0x0000_0FFF_FFFF_FFFF,
    0x0000_1FFF_FFFF_FFFF, 0x0000_3FFF_FFFF_FFFF, 0x0000_7FFF_FFFF_FFFF, 0x0000_FFFF_FFFF_FFFF,
    0x0001_FFFF_FFFF_FFFF, 0x0003_FFFF_FFFF_FFFF, 0x0007_FFFF_FFFF_FFFF, 0x000F_FFFF_FFFF_FFFF,
    0x001F_FFFF_FFFF_FFFF, 0x003F_FFFF_FFFF_FFFF, 0x007F_FFFF_FFFF_FFFF, 0x00FF_FFFF_FFFF_FFFF,
    0x01FF_FFFF_FFFF_FFFF, 0x03FF_FFFF_FFFF_FFFF, 0x07FF_FFFF_FFFF_FFFF, 0x0FFF_FFFF_FFFF_FFFF,
    0x1FFF_FFFF_FFFF_FFFF, 0x3FFF_FFFF_FFFF_FFFF, 0x7FFF_FFFF_FFFF_FFFF, 0xFFFF_FFFF_FFFF_FFFF,
];

/// the number of bits to shift out when finding the index
const INDEX_SHIFT: usize = 6;
/// the flag for selecting
const INDEX_LOWER_MASK: u64 = LOW_SET_FLAGS[INDEX_SHIFT];

/// Represents a bit vector with an very fast index on top of it.
/// For a vector of length N, the index takes up O(N) bits.
pub struct IndexedBitVec {
    index_size: usize,
    bitvec: Vec<u64>,
    index: Vec<u64>
}

impl IndexedBitVec {
    /// Returns an IndexedBitVec with pre-allocated bit-vector and index space.  Note that both the bit-array and index are empty.
    /// # Arguments
    /// * `size` - the number of bits in the bit-vector
    /// # Examples
    /// ```rust
    /// use fmlrc::indexed_bit_vec::IndexedBitVec;
    /// let mut ibv = IndexedBitVec::with_capacity(128);
    /// ```
    #[inline]
	pub fn with_capacity(size: usize) -> Self {
        let index_size: usize = (round::ceil((size as f64) / 64.0, 0) as usize)+1;
        let bitvec = vec![0; index_size];
		Self {
            index_size,
            bitvec,
            index: vec![0; index_size]
        }
    }

    /// Sets a bit in the array. This MUST be used prior to indexing.
    /// # Arguments
    /// * `pos` - the index of the bit to set
    /// # Examples
    /// ```rust
    /// # use fmlrc::indexed_bit_vec::IndexedBitVec;
    /// # let mut ibv = IndexedBitVec::with_capacity(128);
    /// ibv.set_bit(64);
    /// ```
    #[inline]
    pub fn set_bit(&mut self, pos: usize) {
        self.bitvec[pos >> INDEX_SHIFT] |= 0x1 << (pos & INDEX_LOWER_MASK as usize);
    }

    /// Currently, this function is strictly for testing
    #[inline]
    pub fn get_bit(&mut self, pos: usize) -> bool {
        ((self.bitvec[pos >> INDEX_SHIFT] >> (pos & INDEX_LOWER_MASK as usize)) & 0x1) != 0
    }

    /// Builds an index for the array to perform rank queries.
    /// # Arguments
    /// * `initial_rank` - the base rank at bit 0, set to 0 for traditional ranking
    /// # Examples
    /// ```rust
    /// # use fmlrc::indexed_bit_vec::IndexedBitVec;
    /// # let mut ibv = IndexedBitVec::with_capacity(128);
    /// # ibv.set_bit(64);
    /// let initial_rank=0;
    /// ibv.build_index(initial_rank);
    /// ```
    pub fn build_index(&mut self, initial_rank: u64) {
        let mut current_rank = initial_rank;
        for (i, bv_val) in self.bitvec.iter().enumerate().take(self.index_size-1) {
            self.index[i] = current_rank;
            current_rank += bv_val.count_ones() as u64
        }
        self.index[self.index_size-1] = current_rank;
    }

    /// Performs a rank-1 query, returned the number of set bits up to but NOT including `pos`
    /// # Arguments
    /// * `pos` - the position to use for ranking
    /// # Examples
    /// ```rust
    /// # use fmlrc::indexed_bit_vec::IndexedBitVec;
    /// # let mut ibv = IndexedBitVec::with_capacity(128);
    /// # ibv.set_bit(64);
    /// # let initial_rank=0;
    /// # ibv.build_index(initial_rank);
    /// assert_eq!(ibv.rank(64), initial_rank);
    /// assert_eq!(ibv.rank(65), initial_rank+1);
    /// ```
    #[inline]
    pub fn rank(&self, pos: usize) -> u64 {
        let offset: usize = pos >> INDEX_SHIFT;
        self.index[offset] + (self.bitvec[offset] & LOW_SET_FLAGS[pos & INDEX_LOWER_MASK as usize]).count_ones() as u64
    }

    /// Performs a rank-1 query without bounds checking, returned the number of set bits up to but NOT including `pos`
    /// # Arguments
    /// * `pos` - the position to use for ranking
    /// # Safety
    /// If `pos` is outside the length of the bit vector, this will have undefined behavior.
    /// # Examples
    /// ```rust
    /// # use fmlrc::indexed_bit_vec::IndexedBitVec;
    /// # let mut ibv = IndexedBitVec::with_capacity(128);
    /// # ibv.set_bit(64);
    /// # let initial_rank=0;
    /// # ibv.build_index(initial_rank);
    /// unsafe {
    ///   assert_eq!(ibv.rank_unchecked(64), initial_rank);
    ///   assert_eq!(ibv.rank_unchecked(65), initial_rank+1);
    /// }
    /// ```
    #[inline]
    pub unsafe fn rank_unchecked(&self, pos:usize) -> u64 {
        let vec_offset: usize = pos >> INDEX_SHIFT;
        self.index.get_unchecked(vec_offset) + (self.bitvec.get_unchecked(vec_offset) & LOW_SET_FLAGS.get_unchecked(pos & INDEX_LOWER_MASK as usize)).count_ones() as u64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_allocation() {
        //just test index allocation at the corners
        let ibv = IndexedBitVec::with_capacity(63);
        assert_eq!(ibv.index.len(), 2);
        let ibv = IndexedBitVec::with_capacity(64);
        assert_eq!(ibv.index.len(), 2);
        let ibv = IndexedBitVec::with_capacity(65);
        assert_eq!(ibv.index.len(), 3);
    }

    #[test]
    fn test_set_bit() {
        //array with 128 "1"s, then 128 "0"s, then 1 "1"
        let mut ibv = IndexedBitVec::with_capacity(257);
        for i in 0..128 {
            ibv.set_bit(i);
        }
        ibv.set_bit(256);

        //check the sets
        for i in 0..128 {
            assert_eq!(ibv.get_bit(i), true);
        }
        for i in 128..256 {
            assert_eq!(ibv.get_bit(i), false);
        }
        assert_eq!(ibv.get_bit(256), true);

        //check that our slicing is correct also
        let slice: &[u64] = ibv.bitvec.as_slice();
        assert_eq!(slice.len(), 6);
        for i in 0..2 {
            assert_eq!(slice[i], 0xffffffffffffffff);
        }
        for i in 2..4 {
            assert_eq!(slice[i], 0x0);
        }
        assert_eq!(slice[4], 0x1)
    }

    #[test]
    fn test_indexing() {
        //array with 128 "1"s, then 128 "0"s, then 1 "1"
        let mut ibv = IndexedBitVec::with_capacity(257);
        for i in 0..128 {
            ibv.set_bit(i);
        }
        ibv.set_bit(256);

        ibv.build_index(0);
        let slice: &[u64] = ibv.index.as_slice();
        assert_eq!(slice.len(), 6);
        let counts = vec![0, 64, 128, 128, 128, 129];
        for i in 0..slice.len() {
            assert_eq!(counts[i], slice[i]);
        }

        //progressively adding bits
        let mut ibv = IndexedBitVec::with_capacity(64);
        for i in 0..64 {
            ibv.set_bit(i);
            ibv.build_index(0);
            let slice: &[u64] = ibv.index.as_slice();
            assert_eq!(slice[1], (i+1) as u64);
        }
    }

    #[test]
    fn test_offset() {
        //progressively adding bits
        let mut ibv = IndexedBitVec::with_capacity(64);
        for i in 0..64 {
            ibv.set_bit(i);
        }
        ibv.build_index(100);
        let slice: &[u64] = ibv.index.as_slice();
        assert_eq!(slice[1], 164);
    }

    #[test]
    fn test_rank() {
        //array with 128 "1"s, then 128 "0"s, then 1 "1"
        let mut ibv = IndexedBitVec::with_capacity(257);
        for i in 0..128 {
            ibv.set_bit(i);
        }
        ibv.set_bit(256);
        ibv.build_index(0);

        for i in 0..128 {
            assert_eq!(ibv.rank(i), i as u64);
        }
        for i in 128..257 {
            assert_eq!(ibv.rank(i), 128);
        }
        assert_eq!(ibv.rank(257), 129);
    }
}