
/// Returns the median of the array ignoring all values less than `min_value`. 
/// For an odd number of accepted values, the middle value is returned.
/// For an even number of accepted values, the lower value is used in place of the average.
/// If all values are smaller than `min_value`, this function returns 0.
/// Runs in O(n*log(n)) time.
/// # Arguments
/// * `in_array` - the array of values
/// * `min_value` - the minimum allowed value for use in the calculation
/// # Examples
/// ```rust
/// use fmlrc::stats_util::calculate_bounded_median;
/// let test: Vec<u64> = vec![0, 1, 2, 4, 5];
/// assert_eq!(calculate_bounded_median(&test, 0), 2);
/// assert_eq!(calculate_bounded_median(&test, 1), 2);
/// assert_eq!(calculate_bounded_median(&test, 5), 5);
/// //if all values are too small, returns 0
/// assert_eq!(calculate_bounded_median(&test, 6), 0);
/// ```
#[inline]
pub fn calculate_bounded_median(in_array: &[u64], min_value: u64) -> u64 {
    let array_len: u64 = in_array.len() as u64;
    let mut allowed_values: Vec<u64> = Vec::<u64>::with_capacity(array_len as usize);
    
    //filter down to only the values that are usable
    for val in in_array {
        if *val >= min_value {
            allowed_values.push(*val);
        }
    }
    
    match allowed_values.len() {
        0 => 0,
        1 => allowed_values[0],
        count => {
            allowed_values.sort_unstable();
            allowed_values[(count - 1) / 2]
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_median() {
        let test: Vec<u64> = vec![0, 1, 2, 4, 5];
        assert_eq!(calculate_bounded_median(&test, 0), 2);
        assert_eq!(calculate_bounded_median(&test, 1), 2);
        assert_eq!(calculate_bounded_median(&test, 2), 4);
        assert_eq!(calculate_bounded_median(&test, 3), 4);
        assert_eq!(calculate_bounded_median(&test, 4), 4);
        assert_eq!(calculate_bounded_median(&test, 5), 5);

        //if all values are too small, returns 0
        assert_eq!(calculate_bounded_median(&test, 6), 0);
    }
}