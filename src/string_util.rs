
/// contains ASCII to integer encoding
const STRING_TO_INT: [u8; 256] = build_stoi();

/// contains integer to ASCII encoding
const INT_TO_STRING: [u8; 6] = [
    b'$', b'A', b'C', 
    b'G', b'N', b'T'
];

/// for complementing in the integer space; note that $ and N go to themselves
pub const COMPLEMENT_INT: [u8; 6] = [0, 5, 3, 2, 4, 1]; //$ACGNT -> $TGCNA

/// builds up the STRING_TO_INT const for us
const fn build_stoi() -> [u8; 256] {
    let mut ret: [u8; 256] = [4; 256];

    ret['$' as usize] = 0;
    ret['A' as usize] = 1;
    ret['C' as usize] = 2;
    ret['G' as usize] = 3;
    ret['N' as usize] = 4;
    ret['T' as usize] = 5;

    ret['a' as usize] = 1;
    ret['c' as usize] = 2;
    ret['g' as usize] = 3;
    ret['n' as usize] = 4;
    ret['t' as usize] = 5;
    
    ret
}

/// Helper function that reverse complements an integer vector
/// # Arguments
/// * `seq` - the sequence to reverse complement in integer format
/// # Examples
/// ```rust
/// use fmlrc::string_util::reverse_complement_i;
/// let seq: Vec<u8> = vec![0, 1, 2, 3, 4, 5]; //"$ACGNT"
/// let rev_comp = reverse_complement_i(&seq);
/// assert_eq!(rev_comp, vec![1, 4, 2, 3, 5, 0]); //"ANCGT$"
/// ```
#[inline]
pub fn reverse_complement_i(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev() //reverse
        .map(|&c| COMPLEMENT_INT[c as usize]) //complement
        .collect::<Vec<u8>>() //collect and return
}

/// Helper function that converts a string to the corresponding Vec<u8> representation
/// # Arguments
/// * `seq` - the sequence to convert to integer
/// # Examples
/// ```rust
/// use fmlrc::string_util::convert_stoi;
/// let test = "ACGTN$";
/// let converted = convert_stoi(&test);
/// assert_eq!(converted, vec![1, 2, 3, 5, 4, 0]);
/// ```
#[inline]
pub fn convert_stoi(seq: &str) -> Vec<u8> {
    seq.bytes()
        .map(|c| STRING_TO_INT[c as usize])
        .collect::<Vec<u8>>()
}

/// Helper function that converts an integer array to its corresponding String representation
/// # Arguments
/// * `iseq` - the integer sequence to convert to String
/// # Examples
/// ```rust
/// use fmlrc::string_util::convert_itos;
/// let test: Vec<u8> = vec![0, 1, 2, 3, 4, 5];
/// let converted = convert_itos(&test);
/// assert_eq!(&converted, "$ACGNT");
/// ```
#[inline]
pub fn convert_itos(iseq: &[u8]) -> String {
    let ret_vec = iseq.iter()
        .map(|&v| INT_TO_STRING[v as usize])
        .collect::<Vec<u8>>();
    unsafe {
        //this is a no-alloc way to convert the collected Vec<u8> into a String
        String::from_utf8_unchecked(ret_vec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_simple_index() {
        let test: Vec<u8> = vec![0, 1, 2, 3, 4, 5];
        let rev_comp = reverse_complement_i(&test);
        assert_eq!(rev_comp, vec![1, 4, 2, 3, 5, 0]);
    }

    #[test]
    fn test_convert_stoi() {
        let test = "ACGTN$";
        let converted = convert_stoi(&test);
        assert_eq!(converted, vec![1, 2, 3, 5, 4, 0]);
    }

    #[test]
    fn test_convert_itos() {
        let test: Vec<u8> = vec![0, 1, 2, 3, 4, 5];
        let converted = convert_itos(&test);
        assert_eq!(&converted, "$ACGNT");
    }

    #[test]
    fn test_reverse_complement_i() {
        let seq: Vec<u8> = vec![0, 1, 2, 3, 4, 5]; //"$ACGNT"
        let rev_comp = reverse_complement_i(&seq);
        assert_eq!(rev_comp, vec![1, 4, 2, 3, 5, 0]); //"ANCGT$"
    }
}