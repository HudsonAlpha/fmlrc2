
extern crate clap;
extern crate env_logger;
extern crate exitcode;
extern crate log;
extern crate needletail;

use clap::{Arg, App, value_t, values_t};
use log::{info, error};
use needletail::parse_fastx_file;
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
    
    //this is the CLI block, params that get populated appear before
    let bwt_fn: String;
    let long_read_fn: String;
    let corrected_read_fn: String;
    let mut kmer_sizes: Vec<usize> = vec![21, 59];
    let mut threads: usize = 1;
    let mut begin_id: u64 = 0;
    let mut end_id: u64 = 0xFFFFFFFFFFFFFFFF;
    let mut min_count: u64 = 5;
    let mut min_frac: f64 = 0.1;
    let mut branch_factor: f64 = 4.0;
    let mut cache_size: usize = 8;
    let verbose_mode: bool;
    
    let matches = App::new("FMLRC2")
        .version(VERSION.unwrap_or("?"))
        .author("J. Matthew Holt <jholt@hudsonalpha.org>")
        .about("FM-index Long Read Corrector - Rust implementation")
        .arg(Arg::with_name("verbose_mode")
            .short("v")
            .long("verbose")
            .help("enable verbose output"))
        .arg(Arg::with_name("kmer_sizes")
            .short("k")
            .long("K")
            .multiple(true)
            .takes_value(true)
            .help("k-mer sizes for correction, can be specified multiple times (default: \"-k 21 59\")"))
        .arg(Arg::with_name("threads")
            .short("t")
            .long("threads")
            .takes_value(true)
            .help("number of correction threads (default: 1)"))
        .arg(Arg::with_name("begin_id")
            .short("b")
            .long("begin_index")
            .takes_value(true)
            .help("index of read to start with (default: 0)"))
        .arg(Arg::with_name("end_id")
            .short("e")
            .long("end_index")
            .takes_value(true)
            .help("index of read to end with (default: end of file)"))
        .arg(Arg::with_name("min_count")
            .short("m")
            .long("min_count")
            .takes_value(true)
            .help("absolute minimum k-mer count to consisder a path (default: 5)"))
        .arg(Arg::with_name("min_frac")
            .short("f")
            .long("min_dynamic_count")
            .takes_value(true)
            .help("dynamic minimum k-mer count fraction of median to consider a path (default: 0.1)"))
        .arg(Arg::with_name("branch_factor")
            .short("B")
            .long("branch_factor")
            .takes_value(true)
            .help("branching factor for correction, scaled by k (default: 4.0)"))
        .arg(Arg::with_name("cache_size")
            .short("C")
            .long("cache_size")
            .takes_value(true)
            .help("the length of k-mer to precompute in cache (default: 8)"))
        .arg(Arg::with_name("COMP_MSBWT.NPY")
            .help("The compressed BWT file with high accuracy reads")
            .required(true)
            .index(1))
        .arg(Arg::with_name("LONG_READS.FA")
            .help("The FASTX file with uncorrected reads")
            .required(true)
            .index(2))
        .arg(Arg::with_name("CORRECTED_READS.FA")
            .help("The FASTA file to write corrected reads to")
            .required(true)
            .index(3))
        .get_matches();
    
    //pull out required values
    bwt_fn = matches.value_of("COMP_MSBWT.NPY").unwrap().to_string();
    long_read_fn = matches.value_of("LONG_READS.FA").unwrap().to_string();
    corrected_read_fn = matches.value_of("CORRECTED_READS.FA").unwrap().to_string();

    //now check options
    verbose_mode = matches.is_present("verbose_mode");
    kmer_sizes = values_t!(matches.values_of("kmer_sizes"), usize).unwrap_or(kmer_sizes);
    threads = value_t!(matches.value_of("threads"), usize).unwrap_or(threads);
    begin_id = value_t!(matches.value_of("begin_id"), u64).unwrap_or(begin_id);
    end_id = value_t!(matches.value_of("end_id"), u64).unwrap_or(end_id);
    min_count = value_t!(matches.value_of("min_count"), u64).unwrap_or(min_count);
    min_frac = value_t!(matches.value_of("min_frac"), f64).unwrap_or(min_frac);
    branch_factor = value_t!(matches.value_of("branch_factor"), f64).unwrap_or(branch_factor);
    cache_size = value_t!(matches.value_of("cache_size"), usize).unwrap_or(cache_size);

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
    info!("\tcache size: {}", cache_size);
    info!("Correction Parameters:");
    info!("\treads to correct: [{}, {})", begin_id, end_id);
    if begin_id > end_id {
        error!("--begin_index set to value larger than --end_index");
        std::process::exit(exitcode::DATAERR);
    }
    kmer_sizes.sort_unstable();
    info!("\tk-mer sizes: {:?}", kmer_sizes);
    info!("\tabs. mininimum count: {}", min_count);
    info!("\tdyn. minimimum fraction: {}", min_frac);
    if !(0.0..=1.0).contains(&min_frac) {
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
        //TODO: make this 0.4 a CLI param? we did a prelim test with 0.3 and prec/recall went up but total indels/subs went up *shrug*
        midpoint_ed_factor: 0.4,
        tail_buffer_factor: 1.05,
        frac: min_frac,
        verbose: verbose_mode
    };
    let arc_params: Arc<CorrectionParameters> = Arc::new(my_params);

    //first load the BWT into memory
    let mut bwt: BitVectorBWT = BitVectorBWT::with_cache_size(cache_size);
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
    match parse_fastx_file(&long_read_fn) {
        Ok(mut fastx_reader) => {
            while let Some(raw_record) = fastx_reader.next() {
                let record = match raw_record {
                    Ok(record) => { record },
                    Err(e) => {
                        error!("Invalid record while parsing long read file: {:?}", e);
                        std::process::exit(exitcode::IOERR);
                    }
                };
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
                        label: String::from_utf8(record.id().to_vec()).unwrap(),
                        seq: String::from_utf8(record.seq().to_vec()).unwrap()
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
        },
        Err(e) => {
            error!("Failed to open long read file: {:?}", e);
            std::process::exit(exitcode::IOERR);
        }
    }

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
