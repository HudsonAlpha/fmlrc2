
use math::round;

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
            index_size: index_size,
            bitvec: bitvec,
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
    /// unsafe{ ibv.set_bit(64) };
    /// ```
    #[inline]
    pub fn set_bit(&mut self, pos: usize) {
        self.bitvec[pos >> 6] |= (0x1 << (pos & 0x3F));
    }

    /// Currently, this function is strictly for testing
    #[inline]
    fn get_bit(&mut self, pos: usize) -> bool {
        ((self.bitvec[pos >> 6] >> (pos & 0x3F)) & 0x1) != 0
    }

    /// Builds an index for the array to perform rank queries.
    /// # Arguments
    /// * `initial_rank` - the base rank at bit 0, set to 0 for traditional ranking
    /// # Examples
    /// ```rust
    /// # use fmlrc::indexed_bit_vec::IndexedBitVec;
    /// # let mut ibv = IndexedBitVec::with_capacity(128);
    /// # unsafe{ ibv.set_bit(64) };
    /// let initial_rank=0;
    /// ibv.build_index(initial_rank);
    /// ```
    pub fn build_index(&mut self, initial_rank: u64) {
        let slice: &[u64] = self.bitvec.as_slice();
        let mut current_rank = initial_rank;
        for i in {0..self.index_size-1} {
            self.index[i] = current_rank;
            current_rank += slice[i].count_ones() as u64
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
    /// # unsafe{ ibv.set_bit(64) };
    /// # let initial_rank=0;
    /// # ibv.build_index(initial_rank);
    /// assert_eq!(ibv.rank(64), initial_rank);
    /// assert_eq!(ibv.rank(65), initial_rank+1);
    /// ```
    #[inline]
    pub fn rank(&self, pos: usize) -> u64{
        self.index[pos >> 6] + ((self.bitvec.as_slice()[pos >> 6] << (!pos & 0x3f)) << 1).count_ones() as u64
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
        for i in {0..128} {
            unsafe {ibv.set_bit(i);}
        }
        unsafe {ibv.set_bit(256);}

        //check the sets
        for i in {0..128} {
            assert_eq!(ibv.get_bit(i), true);
        }
        for i in {128..25} {
            assert_eq!(ibv.get_bit(i), false);
        }
        assert_eq!(ibv.get_bit(256), true);

        //check that our slicing is correct also
        let slice: &[u64] = ibv.bitvec.as_slice();
        assert_eq!(slice.len(), 6);
        for i in {0..2} {
            assert_eq!(slice[i], 0xffffffffffffffff);
        }
        for i in {2..4} {
            assert_eq!(slice[i], 0x0);
        }
        assert_eq!(slice[4], 0x1)
    }

    #[test]
    fn test_indexing() {
        //array with 128 "1"s, then 128 "0"s, then 1 "1"
        let mut ibv = IndexedBitVec::with_capacity(257);
        for i in {0..128} {
            unsafe {ibv.set_bit(i);}
        }
        unsafe {ibv.set_bit(256);}

        ibv.build_index(0);
        let slice: &[u64] = ibv.index.as_slice();
        assert_eq!(slice.len(), 6);
        let counts = vec![0, 64, 128, 128, 128, 129];
        for i in {0..slice.len()} {
            assert_eq!(counts[i], slice[i]);
        }

        //progressively adding bits
        let mut ibv = IndexedBitVec::with_capacity(64);
        for i in {0..64} {
            unsafe {ibv.set_bit(i);}
            ibv.build_index(0);
            let slice: &[u64] = ibv.index.as_slice();
            assert_eq!(slice[1], (i+1) as u64);
        }
    }

    #[test]
    fn test_offset() {
        //progressively adding bits
        let mut ibv = IndexedBitVec::with_capacity(64);
        for i in {0..64} {
            unsafe {ibv.set_bit(i);}
        }
        ibv.build_index(100);
        let slice: &[u64] = ibv.index.as_slice();
        assert_eq!(slice[1], 164);
    }

    #[test]
    fn test_rank() {
        //array with 128 "1"s, then 128 "0"s, then 1 "1"
        let mut ibv = IndexedBitVec::with_capacity(257);
        for i in {0..128} {
            unsafe {ibv.set_bit(i);}
        }
        unsafe {ibv.set_bit(256);}
        ibv.build_index(0);

        for i in {0..128} {
            assert_eq!(ibv.rank(i), i as u64);
        }
        for i in {128..257} {
            assert_eq!(ibv.rank(i), 128);
        }
        assert_eq!(ibv.rank(257), 129);
    }
}