
pub mod bv_bwt;
pub mod bwt_converter;
pub mod indexed_bit_vec;
pub mod read_correction;
pub mod ropebwt2_util;
pub mod stats_util;
pub mod string_util;

use std::io::Cursor;
use crate::bv_bwt::BitVectorBWT;
use crate::bwt_converter::convert_to_vec;
use crate::ropebwt2_util::create_bwt_from_strings;
use crate::string_util::convert_stoi;
use crate::read_correction::*;

fn get_constant_bwt() -> BitVectorBWT {
    //build the dataset
    let const_string = "AACGGATCAAGCTTACCAGTATTTACGT";
    let rep_count = 30;
    let mut data: Vec<&str> = vec![];
    for i in {0..rep_count} {
        data.push(&const_string);
    }

    //build our BWT
    let bwt_string = create_bwt_from_strings(&data).unwrap();
    let bwt_cursor = Cursor::new(bwt_string);
    let vec = convert_to_vec(bwt_cursor);
    let mut bv_bwt: BitVectorBWT = BitVectorBWT::new();
    bv_bwt.load_vector(vec);
    bv_bwt
}

fn bridge_test() {
    //build our test and verify it's fine
    let bwt = get_constant_bwt();
    let query = convert_stoi(&"AACGGATCAAGCTTACCAGTATTTACGT");
    assert_eq!(bwt.count_kmer(&query), 30);

    //okay, now test the actual bridging
    let seed = convert_stoi(&"AACGGAT");
    let target = convert_stoi(&"TTTACGT");
    let min_count = 15;
    let branch_lim = 1;
    let max_branch_len = 40;
    let bridges = bridge_kmers(&bwt, &seed, &target, min_count, branch_lim, max_branch_len);
    assert_eq!(bridges.len(), 1);
    assert_eq!(bridges[0], query);
}

pub fn test() {
    //println!("Hello, lib!");
    //indexed_bit_vec::test_sub();
    //let test_ibv = indexed_bit_vec::IndexedBitVec::with_capacity(10);
    //test_ibv.test_sub();
    //println!("{} {}", test_ibv.index.len(), test_ibv.index.capacity());
    //println!("{} {}", test_ibv.bitvec.len(), test_ibv.bitvec.capacity());

    /*
    let seq = "TG$$CAGCCG";
    let seq = Cursor::new(seq);
    let vec = bwt_converter::convert_to_vec(seq);
    
    let mut bwt = bv_bwt::BitVectorBWT::new();
    bwt.load_vector(vec);

    println!("{:?}", bwt.count_pileup(vec![2, 3, 3], 2));
    */
    bridge_test();
}