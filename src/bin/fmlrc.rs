
extern crate argparse;
extern crate env_logger;
extern crate exitcode;
extern crate log;
extern crate needletail;

use log::{info, error};
use needletail::parse_sequence_path;
use std::fs::File;
use std::sync::{Arc, mpsc};
use threadpool::ThreadPool;

use fmlrc::bv_bwt::BitVectorBWT;
use fmlrc::ordered_fasta_writer::OrderedFastaWriter;
use fmlrc::read_correction::{CorrectionParameters, CorrectionResults, LongReadFA, correction_job};

const VERSION: Option<&'static str> = option_env!("CARGO_PKG_VERSION");

fn main() {
    //initialize logging for our benefit later
    env_logger::from_env(env_logger::Env::default().default_filter_or("info")).init();

    //non-cli parameters
    const JOB_SLOTS: u64 = 10000;
    const UPDATE_INTERVAL: u64 = 10000;
    let version_string: String = "fmlrc v".to_string()+&VERSION.unwrap_or("?").to_string();

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
    let mut verbose_mode: bool = false;
    {
        let mut ap = argparse::ArgumentParser::new();
        ap.set_description("FM-index Long Read Corrector - Rust implementation");
        //optional parameters
        ap.add_option(&["-v", "--version"], argparse::Print(version_string), "print version number and exit");
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
    match File::open(&bwt_fn) {
        Ok(_) => {},
        Err(e) => {
            error!("Failed to open BWT file: {:?}", e);
            std::process::exit(exitcode::NOINPUT);
        }
    };
    
    info!("\tInput reads: \"{}\"", long_read_fn);
    match File::open(&long_read_fn) {
        Ok(_) => {},
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
    let mut fasta_writer = OrderedFastaWriter::new(&write_file);
    
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
    
    //TODO make some of these hard-coded into params?
    let my_params: CorrectionParameters = CorrectionParameters {
        kmer_sizes,
        min_count,
        max_branch_attempt_length: 10000,
        branch_limit_factor: branch_factor,
        branch_buffer_factor: 1.3,
        tail_buffer_factor: 1.05,
        frac: min_frac,
        verbose: verbose_mode
    };
    let arc_params: Arc<CorrectionParameters> = Arc::new(my_params);

    //first load the BWT into memory
    let mut bwt: BitVectorBWT = BitVectorBWT::new();
    match bwt.load_numpy_file(&bwt_fn) {
        Ok(_) => {},
        Err(e) => {
            error!("Failed to load BWT file: {:?}", e);
            std::process::exit(exitcode::IOERR);
        }
    };
    let arc_bwt: Arc<BitVectorBWT> = Arc::new(bwt);

    //we need to set up the multiprocessing components now
    let pool = ThreadPool::new(threads);
    let (tx, rx) = mpsc::channel();

    //now needletail open the reads to correct
    let mut read_index: u64 = 0;
    let mut jobs_queued: u64 = 0;
    let mut results_received: u64 = 0;

    info!("Starting read correction processes...");
    let parsing_result = parse_sequence_path(
        long_read_fn,
        |_| {},
        |record| {
            if read_index >= begin_id && read_index < end_id {
                //if we've filled our queue, then we should wait until we get some results back
                if jobs_queued - results_received >= JOB_SLOTS {
                    let rx_value: CorrectionResults = rx.recv().unwrap();
                    if verbose_mode {
                        info!("Job #{:?}: {:.2} -> {:.2}", rx_value.read_index, rx_value.avg_before, rx_value.avg_after);
                    }
                    match fasta_writer.write_correction(rx_value) {
                        Ok(()) => {},
                        Err(e) => {
                            error!("Failed while writing read correction: {:?}", e);
                            std::process::exit(exitcode::IOERR);
                        }
                    };
                    results_received += 1;
                    if results_received % UPDATE_INTERVAL == 0 {
                        info!("Processed {} reads...", results_received);
                    }
                }

                //clone the transmit channel and submit the pool job
                let tx = tx.clone();
                let arc_bwt = arc_bwt.clone();
                let arc_params = arc_params.clone();
                let read_data: LongReadFA = LongReadFA {
                    read_index: jobs_queued,
                    label: String::from_utf8((*record.id).to_vec()).unwrap(),
                    seq: String::from_utf8((*record.seq).to_vec()).unwrap()
                };
                //println!("Submitting {:?}", jobs_queued);
                pool.execute(move|| {
                    let correction_results: CorrectionResults = correction_job(arc_bwt, read_data, arc_params);
                    tx.send(correction_results).expect("channel will be there waiting for the pool");
                });
                jobs_queued += 1;
            }
            read_index += 1;
        }
    );

    match parsing_result {
        Ok(_) => {},
        Err(e) => {
            error!("Failed to parse long read file: {:?}", e);
            std::process::exit(exitcode::IOERR);
        }
    };

    while results_received < jobs_queued {
        let rx_value: CorrectionResults = rx.recv().unwrap();
        if verbose_mode {
            info!("Job #{:?}: {:.2} -> {:.2}", rx_value.read_index, rx_value.avg_before, rx_value.avg_after);
        }
        match fasta_writer.write_correction(rx_value) {
            Ok(()) => {},
            Err(e) => {
                error!("Failed while writing read correction: {:?}", e);
                std::process::exit(exitcode::IOERR);
            }
        };
        results_received += 1;
        if results_received % UPDATE_INTERVAL == 0 {
            info!("Processed {} reads...", results_received);
        }
    }
    info!("Finished processing {} total reads in range [{}, {})", results_received, begin_id, end_id);
}
