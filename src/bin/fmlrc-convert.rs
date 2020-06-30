
extern crate argparse;
extern crate env_logger;
extern crate exitcode;

use log::{info, error};
use std::fs::File;
use std::io;

use fmlrc::bwt_converter::{convert_to_vec, save_bwt_numpy};

const VERSION: Option<&'static str> = option_env!("CARGO_PKG_VERSION");

fn main() {
    //initialize logging for our benefit later
    env_logger::from_env(env_logger::Env::default().default_filter_or("info")).init();

    //non-cli parameters
    let version_string: String = "fmlrc-convert v".to_string()+&VERSION.unwrap_or("?").to_string();

    //this is the CLI block, params that get populated appear before
    let mut in_fn: String = "stdin".to_string();
    let mut bwt_fn: String = String::new();
    {
        let mut ap = argparse::ArgumentParser::new();
        ap.set_description("FM-index Long Read Corrector - Rust implementation");
        //optional parameters
        ap.add_option(&["-v", "--version"], argparse::Print(version_string), "print version number and exit");
        ap.refer(&mut in_fn).add_option(&["-i", "--input"], argparse::Store, "The raw BWT (default: stdin)");
        
        //main required parameters
        ap.refer(&mut bwt_fn).add_argument("comp_msbwt.npy", argparse::Store, "The location to store the compressed BWT").required();
        ap.parse_args_or_exit();
    }

    info!("Input parameters (required):");
    info!("\tInput BWT: \"{}\"", in_fn);
    let input_reader: Box<dyn io::Read>;
    if &in_fn == "stdin" {
        input_reader = Box::new(io::stdin());
    } else {
        match File::open(&bwt_fn) {
            Ok(fp) => {
                input_reader = Box::new(fp);
            },
            Err(e) => {
                error!("Failed to open BWT file: {:?}", e);
                std::process::exit(exitcode::NOINPUT);
            }
        };
    }

    info!("\tOutput BWT: \"{}\"", bwt_fn);
    match File::create(&bwt_fn) {
        Ok(file) => file,
        Err(e) => {
            error!("Failed to create output BWT file: {:?}", e);
            std::process::exit(exitcode::NOINPUT);
        }
    };

    //this is where the work happens
    let comp_bwt = convert_to_vec(input_reader);
    save_bwt_numpy(&comp_bwt[..], &bwt_fn).unwrap();

    info!("RLE-BWT conversion complete.");
}