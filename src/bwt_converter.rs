
extern crate log;

use log::{info, error};
use std::io::Read;

pub const VC_LEN: usize = 6;      //$ A C G N T
pub const LETTER_BITS: usize = 3; //defined
pub const NUMBER_BITS: usize = 5; //8-letterBits
pub const NUM_POWER: usize = 32;  //2**numberBits
pub const MASK: u8 = 0x07;        //255 >> numberBits
pub const COUNT_MASK: u8 = 0x1F;

pub fn convert_to_vec(bwt: impl Read) -> Vec<u8> {
    let mut translate: [u8; 256] = [255; 256];
    let valid_symbols = "$ACGNT";
    let mut x = 0;
    for c in valid_symbols.bytes() {
        translate[c as usize] = x;
        x += 1;
    }
    
    let mut ret = Vec::<u8>::new();
    let mut curr: u8 = 36; //'$' - can be any valid character as long as count below is 0
    let mut count: u64 = 0;
    let mut sym_count: [u64; VC_LEN] = [0; VC_LEN];
    for c in bwt.bytes() {
        let ch = c.unwrap();
        if ch == curr {
            count += 1;
        }
        else {
            if translate[ch as usize] == 255 {
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
                    count = count >> NUMBER_BITS;
                }
                curr = ch;
                count = 1;
            }
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
            count = count >> NUMBER_BITS;
        }
    }
    info!("Converted BWT with symbol counts: {:?}", sym_count);

    ret
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

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
}