
extern crate env_logger;
extern crate log;

use log::{info, error};

fn main() {
    //initialize logging for our benefit later
    env_logger::from_env(env_logger::Env::default().default_filter_or("info")).init();

    fmlrc::test();
}
