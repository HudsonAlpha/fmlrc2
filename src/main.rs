
extern crate argparse;
extern crate env_logger;
extern crate exitcode;
extern crate log;

use log::{info, error};
use std::fs::{File, OpenOptions};
use std::path::Path;

fn main() {
    //initialize logging for our benefit later
    env_logger::from_env(env_logger::Env::default().default_filter_or("info")).init();

    //non-cli parameters
    const JOB_SLOTS: usize = 10000;

    //this is the CLI block, params that get populated appear before
    let mut bwt_fn: String = String::new();
    let mut long_read_fn: String = String::new();
    let mut corrected_read_fn: String = String::new();
    let mut kmer_sizes: Vec<usize> = vec![21, 59];
    let mut threads: usize = 1;
    let mut begin_id: u64 = 0;
    let mut end_id: u64 = 0xFFFFFFFFFFFFFFFF;
    let mut min_count: u64 = 5;
    let mut min_frac: f64 = 0.1;
    let mut branch_factor: f64 = 4.0;
    let mut print_version: bool = false;
    let mut verbose_mode: bool = false;
    {
        let mut ap = argparse::ArgumentParser::new();
        ap.set_description("FM-index Long Read Corrector - Rust implementation");
        //optional parameters
        ap.refer(&mut print_version).add_option(&["-v", "--version"], argparse::StoreTrue, "print version number and exit");
        ap.refer(&mut verbose_mode).add_option(&["-V", "--verbose"], argparse::StoreTrue, "enable verbose output");
        ap.refer(&mut kmer_sizes).add_option(&["-k", "-K", "--kmer_size"], argparse::Collect, "k-mer sizes for correction, can be specified multiple times (default: \"-k 21 -K 59\")");
        ap.refer(&mut threads).add_option(&["-t", "--threads"], argparse::Store, "number of correction threads (default: 1)");
        ap.refer(&mut begin_id).add_option(&["-b", "--begin_index"], argparse::Store, "index of read to start with (default: 0)");
        ap.refer(&mut end_id).add_option(&["-e", "--end_index"], argparse::Store, "index of read to end with (default: end of file)");
        ap.refer(&mut min_count).add_option(&["-m", "--min_count"], argparse::Store, "absolute minimum k-mer count to consisder a path (default: 5)");
        ap.refer(&mut min_frac).add_option(&["-f", "--min_dynamic_count"], argparse::Store, "dynamic minimum k-mer count fraction of median to consider a path (default: 0.1)");
        ap.refer(&mut branch_factor).add_option(&["-B", "--branch_factor"], argparse::Store, "branching factor for correction, scaled by k (default: 4.0)");
        
        //main required parameters
        ap.refer(&mut bwt_fn).add_argument("comp_msbwt.npy", argparse::Store, "The compressed BWT file with high accuracy reads").required();
        ap.refer(&mut long_read_fn).add_argument("long_reads.fa", argparse::Store, "The FASTX file with uncorrected reads").required();
        ap.refer(&mut corrected_read_fn).add_argument("corrected_reads.fa", argparse::Store, "The FASTA file to write corrected reads to").required();
        ap.parse_args_or_exit();
    }

    info!("Input parameters (required):");
    info!("\tBWT: \"{}\"", bwt_fn);
    let bwt_file: File = match File::open(&bwt_fn) {
        Ok(file) => file,
        Err(e) => {
            error!("Failed to open BWT file: {:?}", e);
            std::process::exit(exitcode::NOINPUT);
        }
    };
    
    info!("\tInput reads: \"{}\"", long_read_fn);
    let read_file: File = match File::open(&long_read_fn) {
        Ok(file) => file,
        Err(e) => {
            error!("Failed to open input reads file: {:?}", e);
            std::process::exit(exitcode::NOINPUT);
        }
    };

    info!("\tOutput corrected reads: \"{}\"", corrected_read_fn);
    let write_file: File = match File::create(&corrected_read_fn) {
        Ok(file) => file,
        Err(e) => {
            error!("Failed to create output corrected reads file: {:?}", e);
            std::process::exit(exitcode::NOINPUT);
        }
    };
    
    info!("Execution Parameters:");
    info!("\tverbose: {}", verbose_mode);
    info!("\tthreads: {}", threads);
    info!("Correction Parameters:");
    info!("\treads to correct: [{}, {})", begin_id, end_id);
    if begin_id > end_id {
        error!("--begin_index set to value larger than --end_index");
        std::process::exit(exitcode::DATAERR);
    }
    info!("\tk-mer sizes: {:?}", kmer_sizes);
    info!("\tabs. mininimum count: {}", min_count);
    info!("\tdyn. minimimum fraction: {}", min_frac);
    if min_frac < 0.0  || min_frac > 1.0 {
        error!("--min_dynamic_count must be within the range [0, 1]");
        std::process::exit(exitcode::DATAERR);
    }
    info!("\tbranching factor: {}", branch_factor);
    if branch_factor <= 0.0 {
        error!("--branch_factor must be greater than 0.0");
        std::process::exit(exitcode::DATAERR);
    }
    
    //fmlrc::test();
}
