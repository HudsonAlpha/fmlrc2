
use std::cmp::{max,min};
use std::mem::swap;

/// Returns the edit distance between two u8 Vecs by doing the full grid calculation.
/// # Arguments
/// * `v1` - the first Vec
/// * `v2` - the second Vec
/// # Examples
/// ```rust
/// use fmlrc::align::edit_distance;
/// let v1: Vec<u8> = vec![0, 1, 2, 4, 5];
/// let v2: Vec<u8> = vec![0, 1, 3, 4, 5];
/// let v3: Vec<u8> = vec![1, 2, 3, 5];
/// assert_eq!(edit_distance(&v1, &v1), 0);
/// assert_eq!(edit_distance(&v1, &v2), 1);
/// assert_eq!(edit_distance(&v1, &v3), 2);
/// ```
pub fn edit_distance(v1: &[u8], v2: &[u8]) -> usize {
    let l2: usize = v2.len();

    let mut col: Vec<usize> = vec![0; l2+1];
    let mut prev_col: Vec<usize> = (0..l2+1).collect();
    for (i, &c1) in v1.iter().enumerate() {
        col[0] = i+1;
        for (j, &c2) in v2.iter().enumerate() {
            col[j+1] = min(
                min(
                    prev_col[j+1]+1, 
                    col[j]+1
                ), 
                prev_col[j]+({
                    if c1 == c2 {
                        0
                    } else {
                        1
                    }
                })
            );
        }

        swap(&mut col, &mut prev_col);
    }

    prev_col[l2]
}

/// Contains values for a partial match result
#[derive(Debug,PartialEq)]
pub struct MatchScore {
    /// The edit distance of the match
    pub score: usize,
    /// The length of the match with that edit distance
    pub match_length: usize
}

/// Returns the edit distance and index (y) between v1 and a slice of v2[0..y] such that the edit distance is minimized by doing the full grid calculation.
/// This used when you have a correction that is longer than the original sequence and you want to truncate the correction to the closest prefix string.
/// # Arguments
/// * `v1` - the first Vec, this one is always fully used in the edit distance calculation
/// * `v2` - the second Vec, this one will have a prefix used in the final calculation
/// # Examples
/// ```rust
/// use fmlrc::align::{edit_distance_minimize,MatchScore};
/// let v1: Vec<u8> = vec![0, 1, 2, 4, 5];
/// let v2: Vec<u8> = vec![0, 1, 3, 4, 5, 4, 3, 2, 1, 2, 3, 4, 5];
/// assert_eq!(edit_distance_minimize(&v1, &v1), MatchScore {
///     score: 0,
///     match_length: 5
/// });
/// //The best match is v2[0..5] and the edit distance from v1 is 1
/// assert_eq!(edit_distance_minimize(&v1, &v2), MatchScore {
///     score: 1,
///     match_length: 5
/// });
/// ```
pub fn edit_distance_minimize(v1: &[u8], v2: &[u8]) -> MatchScore {
    let l2: usize = v2.len();

    let mut col: Vec<usize> = vec![0; l2+1];
    let mut prev_col: Vec<usize> = (0..l2+1).collect();
    for (i, &c1) in v1.iter().enumerate() {
        col[0] = i+1;
        for (j, &c2) in v2.iter().enumerate() {
            col[j+1] = min(
                min(
                    prev_col[j+1]+1, 
                    col[j]+1
                ), 
                prev_col[j]+({
                    if c1 == c2 {
                        0
                    } else {
                        1
                    }
                })
            );
        }

        swap(&mut col, &mut prev_col);
    }

    let mut min_index = 0;
    for (i, &val) in prev_col.iter().enumerate() {
        if val <= prev_col[min_index] {
            min_index = i;
        }
    }

    MatchScore {
        score: prev_col[min_index],
        match_length: min_index
    }
}

/// Returns the edit distance between two u8 Vecs by using a version of WFA.
/// # Arguments
/// * `v1` - the first Vec
/// * `v2` - the second Vec
/// # Examples
/// ```rust
/// use fmlrc::align::wfa_ed;
/// let v1: Vec<u8> = vec![0, 1, 2, 4, 5];
/// let v2: Vec<u8> = vec![0, 1, 3, 4, 5];
/// let v3: Vec<u8> = vec![1, 2, 3, 5];
/// assert_eq!(wfa_ed(&v1, &v1), 0);
/// assert_eq!(wfa_ed(&v1, &v2), 1);
/// assert_eq!(wfa_ed(&v1, &v3), 2);
/// ```
pub fn wfa_ed(v1: &[u8], v2: &[u8]) -> usize {
    //we need the lengths to know where we are in the vecs
    let l1 = v1.len();
    let l2 = v2.len();

    //stores the next indices that should be compared
    let mut curr_wf: Vec<(usize, usize)> = vec![(0, 0)];
    let mut next_wf: Vec<(usize, usize)> = vec![(0, 0); 3];
    let mut edits = 0;

    //main idea is to iterate until we're at the end of BOTH vecs, this is guaranteed because i and j monotonically increase
    loop {
        //during each iteration, we go over all wavefronts; at iteration e, there are 2*e+1 current wavefronts that will generate 2*(e+1)+1 wavefronts
        //"e" in this context corresponds to the edit distance "edits"
        for (wf_index, &wf) in curr_wf.iter().enumerate() {
            let mut i = wf.0;
            let mut j = wf.1;

            //as long as the symbols match, keep moving along the diagonal
            while i < l1 && j < l2 && v1[i] == v2[j] {
                i += 1;
                j += 1;
            }
            
            if i == l1 && j == l2 {
                //we found the end, return the number of edits required to get here
                return edits;
            }
            else if i == l1 {
                //push the wavefront, but i cannot increase
                next_wf[wf_index] = max(next_wf[wf_index], (i, j));
                next_wf[wf_index+1] = max(next_wf[wf_index+1], (i, j+1));
                next_wf[wf_index+2] = max(next_wf[wf_index+2], (i, j+1));
            } else if j == l2 {
                //push the wavefront, but j cannot increase
                next_wf[wf_index] = max(next_wf[wf_index], (i+1, j));
                next_wf[wf_index+1] = max(next_wf[wf_index+1], (i+1, j));
                next_wf[wf_index+2] = max(next_wf[wf_index+2], (i, j));
            } else {
                //v1 and v2 do not match at i, j; add mismatch, insert, and del to the next wavefront
                next_wf[wf_index] = max(next_wf[wf_index], (i+1, j)); //v2 has a deletion relative to v1
                next_wf[wf_index+1] = max(next_wf[wf_index+1], (i+1, j+1)); //v2 has a mismatch relative to v1
                next_wf[wf_index+2] = max(next_wf[wf_index+2], (i, j+1)); //v2 has an insertion relative to v1
            }
        }

        //we finished this wave, increment the edit count and generate the buffer for the next wavefront
        edits += 1;
        curr_wf = next_wf;
        next_wf = vec![(0, 0); 3+2*edits];
    }
}

/// Returns the edit distance and index (y) between v1 and a slice of v2[0..y] such that the edit distance is minimized by using a version of WFA.
/// This used when you have a correction that is longer than the original sequence and you want to truncate the correction to the closest prefix string.
/// # Arguments
/// * `v1` - the first Vec, this one is always fully used in the edit distance calculation
/// * `v2` - the second Vec, this one will have a prefix used in the final calculation
/// # Examples
/// ```rust
/// use fmlrc::align::{edit_distance_minimize,MatchScore};
/// let v1: Vec<u8> = vec![0, 1, 2, 4, 5];
/// let v2: Vec<u8> = vec![0, 1, 3, 4, 5, 4, 3, 2, 1, 2, 3, 4, 5];
/// assert_eq!(edit_distance_minimize(&v1, &v1), MatchScore {
///     score: 0,
///     match_length: 5
/// });
/// //The best match is v2[0..5] and the edit distance from v1 is 1
/// assert_eq!(edit_distance_minimize(&v1, &v2), MatchScore {
///     score: 1,
///     match_length: 5
/// });
/// ```
pub fn wfa_minimize(v1: &[u8], v2: &[u8]) -> MatchScore {
    //we need the lengths to know where we are in the vecs
    let l1 = v1.len();
    let l2 = v2.len();

    //stores the next indices that should be compared
    let mut curr_wf: Vec<(usize, usize)> = vec![(0, 0)];
    let mut next_wf: Vec<(usize, usize)> = vec![(0, 0); 3];
    let mut edits = 0;

    //main idea is to iterate until we're at the end of v1, this is guaranteed because i abd j are monotonically increasing
    //the loop condition is triggered because we're somewhere in v2 at this point, and it will generally be > 0 (unless a user is doing something real dumb)
    let mut max_v2: usize = 0;
    while max_v2 == 0 {
        //during each iteration, we go over all wavefronts; at iteration e, there are 2*e+1 current wavefronts that will generate 2*(e+1)+1 wavefronts
        //"e" in this context corresponds to the edit distance "edits"
        for (wf_index, &wf) in curr_wf.iter().enumerate() {
            let mut i = wf.0;
            let mut j = wf.1;
            
            //as long as the symbols match, keep moving along the diagonal
            while i < l1 && j < l2 && v1[i] == v2[j] {
                i += 1;
                j += 1;
            }
            
            if i == l1 {
                //we reached the end of v1, see if we are farther along v2 than prior findings
                max_v2 = max(max_v2, j);
            } else if j == l2 {
                //we reached the end of v2, so j can no longer increase
                next_wf[wf_index] = max(next_wf[wf_index], (i+1, j));
                next_wf[wf_index+1] = max(next_wf[wf_index+1], (i+1, j));
                next_wf[wf_index+2] = max(next_wf[wf_index+2], (i, j));

            } else {
                //v1 and v2 do not match at i, j; add mismatch, insert, and del to the next wavefront
                next_wf[wf_index] = max(next_wf[wf_index], (i+1, j)); //v2 has a deletion relative to v1
                next_wf[wf_index+1] = max(next_wf[wf_index+1], (i+1, j+1)); //v2 has a mismatch relative to v1
                next_wf[wf_index+2] = max(next_wf[wf_index+2], (i, j+1)); //v2 has an insertion relative to v1
            }
        }

        //we finished this wave, increment the edit count and generate the buffer for the next wavefront
        edits += 1;
        swap(&mut curr_wf, &mut next_wf);
        next_wf = vec![(0, 0); 3+2*edits];
    }

    MatchScore {
        score: edits-1,
        match_length: max_v2
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_edit_distance() {
        let v1: Vec<u8> = vec![0, 1, 2, 4, 5];
        let v2: Vec<u8> = vec![0, 1, 3, 4, 5];
        let v3: Vec<u8> = vec![1, 2, 3, 5];

        assert_eq!(edit_distance(&v1, &v1), 0);
        assert_eq!(edit_distance(&v1, &v2), 1);
        assert_eq!(edit_distance(&v1, &v3), 2);

        assert_eq!(edit_distance(&v2, &v2), 0);
        assert_eq!(edit_distance(&v2, &v3), 3);
        
        assert_eq!(edit_distance(&v3, &v3), 0);
    }

    #[test]
    fn test_edit_distance_minimize() {
        let v1: Vec<u8> = vec![0, 1, 2, 4, 5];
        let v2: Vec<u8> = vec![0, 1, 3, 4, 5, 4, 3, 2, 1, 2, 3, 4, 5];

        assert_eq!(edit_distance_minimize(&v1, &v1), MatchScore {
            score: 0,
            match_length: 5
        });
        assert_eq!(edit_distance_minimize(&v1, &v2), MatchScore {
            score: 1,
            match_length: 5
        });
    }

    #[test]
    fn test_wfa_ed() {
        let v1: Vec<u8> = vec![0, 1, 2, 4, 5];
        let v2: Vec<u8> = vec![0, 1, 3, 4, 5];
        let v3: Vec<u8> = vec![1, 2, 3, 5];

        assert_eq!(wfa_ed(&v1, &v1), 0);
        assert_eq!(wfa_ed(&v1, &v2), 1);
        assert_eq!(wfa_ed(&v1, &v3), 2);

        assert_eq!(wfa_ed(&v2, &v2), 0);
        assert_eq!(wfa_ed(&v2, &v3), 3);
        
        assert_eq!(wfa_ed(&v3, &v3), 0);
    }

    #[test]
    fn test_wfa_minimize() {
        let v1: Vec<u8> = vec![0, 1, 2, 4, 5];
        let v2: Vec<u8> = vec![0, 1, 3, 4, 5, 4, 3, 2, 1, 2, 3, 4, 5];

        assert_eq!(wfa_minimize(&v1, &v1), MatchScore {
            score: 0,
            match_length: 5
        });
        assert_eq!(wfa_minimize(&v1, &v2), MatchScore {
            score: 1,
            match_length: 5
        });
    }
}