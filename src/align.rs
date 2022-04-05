
use std::cmp::min;
use std::mem::swap;

/// Returns the edit distance between two u8 Vecs.
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

/// Returns the edit distance and index (y) between v1 and a slice of v2[0..y] such that the edit distance is minimized.
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
}