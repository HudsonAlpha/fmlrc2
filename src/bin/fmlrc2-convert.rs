
extern crate clap;
extern crate env_logger;
extern crate exitcode;

use clap::{Arg, App, value_t};
use log::{info, error};
use std::fs::File;
use std::io;

use fmlrc::bwt_converter::{convert_to_vec, save_bwt_numpy};

const VERSION: Option<&'static str> = option_env!("CARGO_PKG_VERSION");

fn main() {
    //initialize logging for our benefit later
    env_logger::from_env(env_logger::Env::default().default_filter_or("info")).init();

    //this is the CLI block, params that get populated appear before
    let mut in_fn: String = "stdin".to_string();
    let bwt_fn: String;

    let matches = App::new("FMLRC2 BWT Converter")
        .version(VERSION.unwrap_or("?"))
        .author("J. Matthew Holt <jholt@hudsonalpha.org>")
        .about("FMLCR BWT Converter - Rust implementation")
        .arg(Arg::with_name("in_fn")
            .short("i")
            .long("--input")
            .takes_value(true)
            .help("The raw BWT (default: stdin)"))
        .arg(Arg::with_name("COMP_MSBWT.NPY")
            .help("The location to store the compressed BWT")
            .required(true)
            .index(1))
        .get_matches();
    
    //pull out required values
    bwt_fn = matches.value_of("COMP_MSBWT.NPY").unwrap().to_string();
    
    //optional values
    in_fn = value_t!(matches.value_of("in_fn"), String).unwrap_or(in_fn);

    info!("Input parameters (required):");
    info!("\tInput BWT: \"{}\"", in_fn);
    let input_reader: Box<dyn io::Read> = if &in_fn == "stdin" {
        Box::new(io::stdin())
    } else {
        Box::new(match File::open(&in_fn) {
            Ok(fp) => {
                fp
            },
            Err(e) => {
                error!("Failed to open BWT file: {:?}", e);
                std::process::exit(exitcode::NOINPUT);
            }
        })
    };

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