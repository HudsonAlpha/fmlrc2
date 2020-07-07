
/*!
# FM-Index Long Read Corrector v2
This library provides access to the functionality used by FMLRC2 to perform read correction using a Burrows Wheeler Transform (BWT).
Currently, the BWT is assumed to have been generated externally (typically with a tool like ropebwt2) and stored in the same numpy format as FMLRC v1.
FMLRC load a binary representation of the BWT into memory for performing very fast queries at the cost of memory usage.
This particular implementation is accelerated over FMLRC v1 by using a cache to pre-compute common queries to the BWT.

## Example
```rust
use fmlrc::bv_bwt::BitVectorBWT;
use fmlrc::bwt_converter::convert_to_vec;
use fmlrc::ropebwt2_util::create_bwt_from_strings;
use fmlrc::string_util::convert_stoi;
use std::io::Cursor;

//example with in-memory BWT
let data: Vec<&str> = vec!["ACGT", "CCGG"];
let seq = create_bwt_from_strings(&data).unwrap();
let cursor_seq = Cursor::new(seq);
let vec_form = convert_to_vec(cursor_seq);
let mut bwt = BitVectorBWT::new();
bwt.load_vector(vec_form);
//bwt.load_numpy_file(filename); <- if in a numpy file

//do a count
let kmer: Vec<u8> = convert_stoi(&"ACGT");
let kmer_count = bwt.count_kmer(&kmer); //ACGT
assert_eq!(kmer_count, 1);
```
*/

/// Contains the bit vector implementation of the BWT
pub mod bv_bwt;
/// Contains the function for reformating a BWT string into the expected run-length format or numpy file
pub mod bwt_converter;
/// Contains bit vector with basic rank support; other crates exist with this, but they tended to be slow for some reason
pub mod indexed_bit_vec;
/// Contains a wrapper around the rust-bio FASTA writer, but forces an ordering on the reads
pub mod ordered_fasta_writer;
/// Contains the logic for performing the read correction
pub mod read_correction;
/// Contains wrapper functions for `ropebwt2`, most will fail if `ropebwt2` is not on the PATH
pub mod ropebwt2_util;
/// Contains special statistics functions, mainly an ignored median score
pub mod stats_util;
/// Contains inline functions for converting between strings and integer formats
pub mod string_util;
