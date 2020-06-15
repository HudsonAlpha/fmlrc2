
extern crate subprocess;

use std::io::Read;
use subprocess::{Exec, Redirection};

/// Wrapper script for using ropebwt2 with collection of strings.
/// This is primarily for performing easy tests within the Rust environment.
/// For production, we recommend running the `ropebwt2` command separately.
/// # Arguments
/// * `data` - a string instance where individual strings are separated by newline characters ('\n')
/// # Examples
/// ```rust
/// use fmlrc::ropebwt2_util::create_bwt_from_strings;
/// let data: Vec<&str> = vec!["CCGT", "ACG"];
/// assert_eq!(create_bwt_from_strings(&data).unwrap(), "GT$$ACCCG\n".to_string());
/// ```
pub fn create_bwt_from_strings(data: &Vec<&str>) -> Result<String, Box<dyn std::error::Error>> {
    let join_data = data.join("\n");
    let out = (
        Exec::cmd("sort") |
        Exec::cmd("tr").arg("NT").arg("TN") |
        Exec::cmd("ropebwt2").arg("-LR") |
        Exec::cmd("tr").arg("NT").arg("TN")
    ).stdin(join_data.as_str()).stdout(Redirection::Pipe).capture()?.stdout_str();
    Ok(out)
}
/*
TODO: seems like we need a function that will take a list of compressed FASTQ files and write the BWT out
    will also need a way to test it; looks like a combination of tempfile and flate2 to create gzipped files, then
    rust-bio to write fastq records easily
    Follow up links:
    https://docs.rs/tempfile/3.1.0/tempfile/struct.Builder.html
    https://docs.rs/bio/0.31.0/bio/io/fastq/index.html
    https://docs.rs/flate2/1.0.14/flate2/struct.GzBuilder.html
    https://rust-cli.github.io/book/tutorial/testing.html
*/

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_from_strings() {
        let data: Vec<&str> = vec!["CCGT", "ACG", "N"];
        assert_eq!(create_bwt_from_strings(&data).unwrap(), "GTN$$ACCC$G\n".to_string());
    }
}