
extern crate subprocess;
extern crate log;

use log::{info, error};
use std::io::Read;
use subprocess::{Exec, Redirection};

/// Wrapper function for using ropebwt2 with collection of strings.
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
    info!("Concatenating string vector...");
    let join_data = data.join("\n");
    info!("Running BWT construction pipeline from strings...");
    let out = (
        Exec::cmd("sort") |
        Exec::cmd("tr").arg("NT").arg("TN") |
        Exec::cmd("ropebwt2").arg("-LR") |
        Exec::cmd("tr").arg("NT").arg("TN")
    ).stdin(join_data.as_str()).stdout(Redirection::Pipe).capture()?.stdout_str();
    info!("BWT construction from strings complete.");
    Ok(out)
}

/// Wrapper function for using ropebwt2 with a collection of gzipped FASTQ files.
/// Returns a streamable result from the ropebwt2 command.
/// # Argument
/// * `fastqs` - a vector of fastq filenames
pub fn stream_bwt_from_fastqs(fastqs: &Vec<&str>) -> Result<Box<dyn Read>, Box<dyn std::error::Error>> {
    let join_filenames = fastqs.join(" ");
    info!("Streaming BWT construction pipeline from gzipped FASTQ files:\n{:?}", fastqs);
    let mut initial_command = Exec::cmd("gunzip").arg("-c");
    for fq in fastqs {
        initial_command = initial_command.arg(fq);
    }
    let out = (
        initial_command |
        Exec::cmd("awk").arg("NR % 4 == 2") |
        Exec::cmd("sort") |
        Exec::cmd("tr").arg("NT").arg("TN") |
        Exec::cmd("ropebwt2").arg("-LR") |
        Exec::cmd("tr").arg("NT").arg("TN")
    ).stdout(Redirection::Pipe).stream_stdout()?;
    Ok(Box::new(out))
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
    use flate2::{Compression, GzBuilder};
    use std::fs::File;
    use std::io::{self, Write};
    use tempfile::{Builder, NamedTempFile};
    
    #[test]
    fn test_from_strings() {
        let data: Vec<&str> = vec!["CCGT", "ACG", "N"];
        assert_eq!(create_bwt_from_strings(&data).unwrap(), "GTN$$ACCC$G\n".to_string());
    }

    fn write_strings_to_fqgz(data: Vec<&str>) -> NamedTempFile {
        let mut file: NamedTempFile = Builder::new().prefix("temp_data_").suffix(".fq.gz").tempfile().unwrap();
        let mut gz = GzBuilder::new().write(file, Compression::default());
        let mut i: usize = 0;
        for s in data {
            writeln!(gz, "@seq_{}\n{}\n{}\n{}", i, s, "+", "F".repeat(s.len())).unwrap();
            i += 1;
        }

        //have to keep the file handle or everything blows up
        gz.finish().unwrap()
    }

    #[test]
    fn test_from_fastqs() {
        //create a temporary file and get the filename out
        let data: Vec<&str> = vec!["CCGT", "N"];
        let file = write_strings_to_fqgz(data);
        let data2: Vec<&str> = vec!["ACG"];
        let file2 = write_strings_to_fqgz(data2);

        let fastq_filenames: Vec<&str> = vec![
            &file.path().to_str().unwrap(),
            &file2.path().to_str().unwrap()
        ];
        
        //now read the string in and verify correctness
        let mut bwt_stream = stream_bwt_from_fastqs(&fastq_filenames).unwrap();
        let mut buffer: String = String::new();
        let num_bytes = bwt_stream.read_to_string(&mut buffer).unwrap();
        assert_eq!(num_bytes, 12);
        assert_eq!(buffer, "GTN$$ACCC$G\n".to_string());
    }
}