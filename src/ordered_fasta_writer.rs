
use bio::io::fasta;
use std::collections::HashMap;
use std::fs;
use std::io;

use crate::read_correction::CorrectionResults;

/// This is a wrapper for the bio::io::fasta::Writer that forces reads to be written in a specified order.
/// Reads that are not ready to be written are stored until the reads before it are written.
/// # Examples
/// ```rust
/// use std::fs;
/// use tempfile::{Builder, NamedTempFile};
/// use fmlrc::ordered_fasta_writer::OrderedFastaWriter;
/// use fmlrc::read_correction::CorrectionResults;
/// 
/// // File setup
/// let file: NamedTempFile = Builder::new().prefix("out_fasta_").suffix(".fa").tempfile().unwrap();
/// let temp_filename: String = file.path().to_str().unwrap().to_string(); 
/// {
///     let mut fasta_writer = OrderedFastaWriter::new(&file);
///
///     //write some corrections out of order
///     let corr_result_a = CorrectionResults {
///         read_index: 0,
///         label: "a".to_string(),
///         original_seq: "GCTA".to_string(),
///         corrected_seq: "CCTA".to_string(),
///         avg_before: 0.0,
///         avg_after: 0.0
///     };
///     let corr_result_b = CorrectionResults {
///         read_index: 1,
///         label: "b".to_string(),
///         original_seq: "ACGT".to_string(),
///         corrected_seq: "ACCT".to_string(),
///         avg_before: 0.0,
///         avg_after: 0.0
///     };
///     fasta_writer.write_correction(corr_result_b).unwrap();
///     fasta_writer.write_correction(corr_result_a).unwrap();
///     //call flush manually if the file writer stays in scope
///     //fasta_writer.flush().unwrap();
/// }
/// //check that the results are in order
/// let file_string = fs::read_to_string(temp_filename).unwrap();
/// let expected = ">a\nCCTA\n>b\nACCT\n";
/// assert_eq!(expected, file_string);
/// ```
pub struct OrderedFastaWriter<W: io::Write> {
    /// the actual fasta writer
    writer: fasta::Writer<W>,
    /// contains results we aren't ready to write yet
    map_store: HashMap<u64, CorrectionResults>,
    /// the index for the next read to write
    current_index: u64
}

impl <W: io::Write> OrderedFastaWriter<W> {
    /// Creates an `OrderedFastaWriter` that is wrapping a file buffer.
    /// # Arguments
    /// `writer` - a buffer implementing `std::io::Write`
    pub fn new(writer: W) -> Self {
        OrderedFastaWriter {
            writer: fasta::Writer::new(writer),
            map_store: HashMap::<u64, CorrectionResults>::new(),
            current_index: 0
        }
    }

    /// Writes a correction to the file or buffers it if not ready to write.
    /// # Arguments
    /// `correction` - a read correction from fmlrc containing the read index, label, and corrected sequence
    pub fn write_correction(&mut self, correction: CorrectionResults) -> io::Result<()> {
        if correction.read_index < self.current_index {
            return Err(io::Error::new(io::ErrorKind::Other, "Read index is smaller than next expected index"));
        }
        match self.map_store.insert(correction.read_index, correction) {
            None => {},
            Some(_) => {
                return Err(io::Error::new(io::ErrorKind::Other, "Read index was already present in the map_store"));
            }
        };
        self.drain_map_store()
    }

    fn drain_map_store(&mut self) -> io::Result<()> {
        while !self.map_store.is_empty() {
            match self.map_store.remove(&self.current_index) {
                Some(correction) => {
                    let record: fasta::Record = fasta::Record::with_attrs(&correction.label, None, &correction.corrected_seq.as_bytes());
                    match self.writer.write_record(&record) {
                        Ok(()) => {},
                        Err(e) => return Err(e)
                    };
                    self.current_index += 1;
                },
                None => {
                    break;
                }
            };
        }
        Ok(())
    }

    /// Flushes the buffer, call before trying to read anything.
    /// # Arguments
    /// `correction` - a read correction from fmlrc containing the read index, label, and corrected sequence
    pub fn flush(&mut self) -> io::Result<()> {
        match self.drain_map_store() {
            Ok(()) => {},
            Err(e) => return Err(e)
        };
        self.writer.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::{Builder, NamedTempFile};

    #[test]
    fn test_fasta_writing() {
        let file: NamedTempFile = Builder::new().prefix("out_fasta_").suffix(".fa").tempfile().unwrap();
        let temp_filename: String = file.path().to_str().unwrap().to_string();
        {
            let mut ofw = OrderedFastaWriter::new(&file);
            
            let correction_a = CorrectionResults {
                read_index: 0,
                label: "a".to_string(),
                original_seq: "GCTA".to_string(),
                corrected_seq: "CCTA".to_string(),
                avg_before: 0.0,
                avg_after: 0.0
            };

            let correction_b = CorrectionResults {
                read_index: 1,
                label: "b".to_string(),
                original_seq: "ACGT".to_string(),
                corrected_seq: "ACCT".to_string(),
                avg_before: 0.0,
                avg_after: 0.0
            };

            let correction_c = CorrectionResults {
                read_index: 2,
                label: "c".to_string(),
                original_seq: "ACGT".to_string(),
                corrected_seq: "GGTT".to_string(),
                avg_before: 0.0,
                avg_after: 0.0
            };
            
            //submit them out of order, b -> c -> a
            ofw.write_correction(correction_b).unwrap();
            ofw.write_correction(correction_c).unwrap();
            ofw.write_correction(correction_a).unwrap();
        }
        
        let file_string = fs::read_to_string(temp_filename).unwrap();
        let expected = ">a\nCCTA\n>b\nACCT\n>c\nGGTT\n";

        assert_eq!(expected, file_string);
    }
}