[![Crates.io](https://img.shields.io/crates/d/fmlrc.svg)](https://crates.io/crates/fmlrc)
[![Crates.io](https://img.shields.io/crates/v/fmlrc.svg)](https://crates.io/crates/fmlrc)
[![Crates.io](https://img.shields.io/crates/l/fmlrc.svg)](https://crates.io/crates/fmlrc)
[![Build status](https://github.com/HudsonAlpha/rust-fmlrc/actions/workflows/quickstart-ci.yml/badge.svg)](https://github.com/HudsonAlpha/rust-fmlrc/actions)

# rust-fmlrc

This repo contains the source code for FMLRC v2, based on the same methodology used by the original [FMLRC](https://github.com/holtjma/fmlrc).  In benchmarks, the results between FMLRC v1 and v2 are nearly identical, but tests have shown that v2 uses approximately 50% of the run and CPU time compared to v1.

## Installation
All installation options assume you have installed [Rust](https://www.rust-lang.org) along with the `cargo` crate manager for Rust.

### From Cargo
```bash
cargo install fmlrc
fmlrc2 -h
fmlrc2-convert -h
```

### From GitHub
```bash 
git clone https://github.com/HudsonAlpha/rust-fmlrc.git
cd rust-fmlrc
#testing optional, some tests will fail if ropebwt2 is not installed or cannot be found on PATH
cargo test --release
cargo build --release
./target/release/fmlrc2 -h
./target/release/fmlrc2-convert -h
```

## Usage
### BWT Building
#### msbwt2 Construction Approach
For most users, it is recommended to use the [msbwt2](https://github.com/HudsonAlpha/rust-msbwt) crate to [build the BWT](https://github.com/HudsonAlpha/rust-msbwt#msbwt-building).
This approach is generally simpler (requiring only one command) and more flexible (accepting both FASTQ and FASTA files at once).
While it is generally competitive with the `ropebwt2` construction approach (see below) for memory and CPU usage, it is not parallelized and typically runs slower by wall-clock time.

#### ropebwt2 Construction Approach
If you are familiar with more complicated shell commands, then `ropebwt2` can also be used to build the BWT.  
For _most_ short-read datasets, this approach is faster than `msbwt2-build` but also more complicated (multiple piped commands) and less flexible (fixed to FASTQ in the below example).
Given one or more FASTQ files of accurate reads (`reads.fq.gz` with extras labeled as `[reads2.fq.gz ...]`), you can use the following command from this crate to create a BWT at `comp_msbwt.npy`.
Note that this command requires the [ropebwt2](https://github.com/lh3/ropebwt2) executable to be installed:
```
gunzip -c reads.fq.gz [reads2.fq.gz ...] | \
    awk 'NR % 4 == 2' | \
    sort | \
    tr NT TN | \
    ropebwt2 -LR | \
    tr NT TN | \
    fmlrc2-convert comp_msbwt.npy
```

Note: If you are **only** using the BWT for correction, then the `sort` can be removed from the above command. This will reduce construction time significantly, but loses the read recovery property of the BWT.

### Correction
Assuming the accurate-read BWT is built (`comp_msbwt.npy`) and uncorrected reads are available (fasta/fastq, gzip optional, `uncorrected.fq.gz`), invoking FMLRC v2 is fairly simple:
```
fmlrc2 [OPTIONS] <comp_msbwt.npy> <uncorrected.fq.gz> <corrected_reads.fa>
```

Currently, only uncompressed FASTA is supported for output reads.

#### Options to consider
1. `-h` - see full list of options and exit
2. `-k`, `--K` - sets the k-mer sizes to use, default is `[21, 59]`; all values are sorted from lowest to highest prior to correction
3. `-t`, `--threads` - number of correction threads to use (default: 1)
4. `-C`, `--cache_size` - the length of sequences to pre-compute (i.e. `C`-mers); will reduce CPU-time of queries by `O(C)` but *increases* cache memory usage by `O(6^C)`; default of `8` uses ~25MB; if memory is not an issue, consider using `10` with ~1GB cache footprint (or larger if memory _really_ isn't an issue)

## FMLRC v2 core differences
1. Implemented in [Rust](https://www.rust-lang.org) instead of C++ - this comes will all the benefits of Rust including `cargo`, such as easy installation of the binary and supporting structs/functions along with documentation
2. Unlimited `k`/`K` parameters - FMLRC v1 allowed 1 or 2 sizes for `k` only; FMLRC v2 can have the option set as many times as desired at increased CPU time (for example, a 3-pass correction with `k=[21, 59, 79]`) 
3. Call caching - FMLRC v2 pre-computes all _k_-mers of a given size. This reduces the run-time significantly by cutting reducing calls to the FM-index.
4. Input handling - thanks to [needletail](https://crates.io/crates/needletail), the uncorrected reads can be in FASTA/FASTQ and may or may not be gzip compressed.
5. SIMD accelerated alignment - thanks to [triple_accel](https://crates.io/crates/triple_accel), the correction alignment step can be accelerated with SIMD instructions when available
6. Unit testing - FMLRC v2 has unit testing through the standard Rust testing framework (i.e. `cargo test`)

## Benchmarks
Thus far, all benchmarks have focused on a relatively small _E. coli_ dataset for verifying correctness.
The files for this dataset can be found in the original [fmlrc example](https://github.com/holtjma/fmlrc/blob/master/example/run_example.sh).
The exact same BWT and uncorrected long read files were used for both fmlrc v1 and fmlrc v2.
[ELECTOR](https://github.com/kamimrcht/ELECTOR) was used to evaluate the results.
All fmlrc executions were run on a Macbook Pro with 2.2 GHz Intel Core i7 processor (8 cores) with 16 GB of RAM.
Run times were gathered using Mac OSX `time`.
All parameters were set to defaults except for `-C 10` in FMLRC v2.

The following table summarizes the results. 
The actual corrections are _nearly_ identical (there are slight differences not reflected in summary metrics).
However, FMLRC v2 runs in less than half the time from both real time and CPU time perspectives. 
While not explicitly measured, FMLRC v2 does use ~1GB of extra memory due to the 10-mer cache (`-C 10`).

| Metric | FMLRC v1.0.0 | FMLRC2 v0.1.5 (`-C 10`) |
| - | - | - |
| Recall | 0.9825 | 0.9826 |
| Precision | 0.9815 | 0.9816 |
| Real time | 7m29.908s | **2m38.214s** |
| CPU time | 51m34.704s | **17m54.289s** |

## Reference
FMLRC v2 does not currently have a pre-print or paper. If you use FMLRC v2, please cite the FMLRC v1 paper:

[Wang, Jeremy R. and Holt, James and McMillan, Leonard and Jones, Corbin D. FMLRC: Hybrid long read error correction using an FM-index. BMC Bioinformatics, 2018. 19 (1) 50.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2051-3)

## License
Licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license
   ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Contribution
Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.