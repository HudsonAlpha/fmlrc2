
extern crate log;

use log::info;
use std::fs::{File, OpenOptions};
use std::io::prelude::*;
use std::io::{BufWriter, Read};

pub const VC_LEN: usize = 6;      //$ A C G N T
pub const LETTER_BITS: usize = 3; //defined
pub const NUMBER_BITS: usize = 5; //8-letterBits
pub const NUM_POWER: usize = 32;  //2**numberBits
pub const MASK: u8 = 0x07;        //255 >> numberBits
pub const COUNT_MASK: u8 = 0x1F;

//TODO: convert_to_vec currently pulls the whole compressed BWT into memory, which I'm largely okay with
// if we ever decide to change it, we could try to build an iterator following this guide: 
// https://burgers.io/wrapped-iterators-in-rust; seems like a bit of a pain though

/// this function will convert a stream of characters into the compressed vector representation
/// # Arguments
/// * `bwt` - the stream of characters to be compressed, allowed characters are "$ACGNT"; newline characters ('\n') are ignored
/// # Examples
/// ```rust
/// use fmlrc::bwt_converter::convert_to_vec;
/// use std::io::Cursor;
/// let seq = Cursor::new("ACGNT$");
/// let vec = convert_to_vec(seq);
/// assert_eq!(vec.len(), 6);
/// ```
pub fn convert_to_vec(bwt: impl Read) -> Vec<u8> {
    let mut translate: [u8; 256] = [255; 256];
    let valid_symbols = "$ACGNT";
    for (x, c) in valid_symbols.bytes().enumerate() {
        translate[c as usize] = x as u8;
    }
    
    let mut ret = Vec::<u8>::new();
    let mut curr: u8 = 36; //'$' - can be any valid character as long as count below is 0
    let mut count: u64 = 0;
    let mut sym_count: [u64; 6] = [0; 6];
    for c in bwt.bytes() {
        let ch = c.unwrap();
        if ch == curr {
            count += 1;
        }
        else if translate[ch as usize] == 255 {
            //10 is newline, which we can ignore
            if ch != 10 {
                panic!("Unexpected symbol in input: char \"{}\"", curr);
            }
        }
        else {
            //end of run, add the total to the vector
            //debug stuff - symCount[translator[currSym]] += currCount;
            sym_count[translate[curr as usize] as usize] += count;
            while count > 0 {
                let write_byte: u8 = translate[curr as usize] | ((count as u8 & COUNT_MASK) << LETTER_BITS) as u8;
                ret.push(write_byte);
                count >>= NUMBER_BITS;
            }
            curr = ch;
            count = 1;
        }
    }

    //do the last run now
    if translate[curr as usize] == 255 {
        panic!("Unexpected symbol in input: char \"{}\"", curr);
    }
    else {
        //end of run, add the total to the vector
        //debug stuff - symCount[translator[currSym]] += currCount;
        sym_count[translate[curr as usize] as usize] += count;
        while count > 0 {
            let write_byte: u8 = translate[curr as usize] | ((count as u8 & COUNT_MASK) << LETTER_BITS) as u8;
            ret.push(write_byte);
            count >>= NUMBER_BITS;
        }
    }
    info!("Converted BWT with symbol counts: {:?}", sym_count);
    info!("RLE-BWT byte length: {:?}", ret.len());

    ret
}

/// This will take some iterable byte array and save it to a file in numpy format. 
/// This primarily adds a numpy data type and shape fields that occupies the first 96 bytes of the file.
/// The intended use is to pass in a compressed BWT for saving them, but really any byte array can be stored this way.
/// # Arguments
/// * `bwt` - a data type implementing Read that represents the compressed BWT
/// * `filename` - the filename to save the output to
pub fn save_bwt_numpy(bwt: impl Read, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
    let npy_file: File = File::create(filename)?;
    let mut buffer = BufWriter::new(npy_file);
    buffer.write_all(&[32; 95])?;
    buffer.write_all(&[10; 1])?;
    let mut num_bytes: u64 = 0;
    for c in bwt.bytes() {
        buffer.write_all(&[c?])?;
        num_bytes += 1;
    }
    buffer.flush()?;

    let header_string = b"\x93NUMPY\x01\x00\x56\x00{\'descr\': \'|u1\', \'fortran_order\': False, \'shape\': (";
    let header_tail = b", ), }"; //added a space after ',' here, so slightly different but functionally identical
    let mut npy_file: File = OpenOptions::new().write(true).open(filename)?;

    //header format - "header_string" -> length of data -> "header_tail"
    npy_file.write_all(header_string)?;
    npy_file.write_all(num_bytes.to_string().as_bytes())?;
    npy_file.write_all(header_tail)?;
    npy_file.flush()?;
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ropebwt2_util::create_bwt_from_strings;
    use std::io::Cursor;
    use tempfile::{Builder, NamedTempFile};

    #[test]
    fn test_convert_to_vec() {
        let seq = "ACGNT$";
        let seq = Cursor::new(seq);
        let vec = convert_to_vec(seq);
        assert_eq!(vec.len(), 6);
        assert_eq!(vec[0], 8+1);
        assert_eq!(vec[1], 8+2);
        assert_eq!(vec[2], 8+3);
        assert_eq!(vec[3], 8+4);
        assert_eq!(vec[4], 8+5);
        assert_eq!(vec[5], 8+0);
    }

    #[test]
    fn test_newline() {
        //test newlines at start, end, mid-run, and between chars
        let seq = "\n$$\n$$\nAAA\n";
        let seq = Cursor::new(seq);
        let vec = convert_to_vec(seq);
        assert_eq!(vec.len(), 2);
        assert_eq!(vec[0], (4<<3)+0);
        assert_eq!(vec[1], (3<<3)+1);
    }

    #[test]
    fn test_compression() {
        let seq = "A".repeat(32+32*32*3);
        let seq = Cursor::new(seq);
        let vec = convert_to_vec(seq);
        assert_eq!(vec.len(), 3);
        assert_eq!(vec[0], 1);
        assert_eq!(vec[1], 9);
        assert_eq!(vec[2], 1+(3 << 3));

        let seq = "A".repeat(31)+&("C".repeat(31));
        let seq = Cursor::new(seq);
        let vec = convert_to_vec(seq);
        assert_eq!(vec.len(), 2);
        assert_eq!(vec[0], 249);
        assert_eq!(vec[1], 250);

        let seq = "N".repeat(32767);
        let seq = Cursor::new(seq);
        let vec = convert_to_vec(seq);
        assert_eq!(vec.len(), 3);
        assert_eq!(vec[0], 4+(0x1F << 3));
        assert_eq!(vec[1], 4+(0x1F << 3));
        assert_eq!(vec[2], 4+(0x1F << 3));
    }

    #[test]
    fn test_bwt_conversion() {
        //bwt = "GTN$$ACCC$G" -> [11, 13, 12, 16, 9, 26, 8, 11]
        let data: Vec<&str> = vec!["CCGT", "N", "ACG"];
        let bwt: String = create_bwt_from_strings(&data).unwrap();
        let seq = Cursor::new(bwt);
        let vec = convert_to_vec(seq);

        let expected_result: Vec<u8> = vec![11, 13, 12, 16, 9, 26, 8, 11];
        assert_eq!(expected_result, vec);
    }

    #[test]
    fn test_save_bwt_numpy() {
        let seq = "A".repeat(32+32*32*3);
        let header_string = b"\x93NUMPY\x01\x00\x56\x00{\'descr\': \'|u1\', \'fortran_order\': False, \'shape\': (3, ), }";
        let mut expected_result: Vec<u8> = header_string.to_vec();
        while expected_result.len() < 95 {
            expected_result.push(32);
        }
        expected_result.push(10);
        expected_result.push(1);
        expected_result.push(9);
        expected_result.push(1+(3 << 3));

        let seq = Cursor::new(seq);
        let vec = convert_to_vec(seq);

        let file: NamedTempFile = Builder::new().prefix("temp_data_").suffix(".npy").tempfile().unwrap();
        let filename: String = file.path().to_str().unwrap().to_string();
        
        save_bwt_numpy(&vec[..], &filename).unwrap();
        let read_file = File::open(&filename).unwrap();
        let mut read_result: Vec<u8> = Vec::<u8>::new();
        for c in read_file.bytes() {
            read_result.push(c.unwrap());
        }
        assert_eq!(expected_result, read_result);
    }
}