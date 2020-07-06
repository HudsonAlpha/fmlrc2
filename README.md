# rust-fmlrc

This repo contains the source code for FMLRC v2, based on the same methodology used by the original [FMLRC](https://github.com/holtjma/fmlrc).  In benchmarks, the results between FMLRC v1 and v2 are nearly identical, but tests have shown that v2 uses approximately 50% of the run and CPU time compared to v1.

## Installation
All installation options assume you have installed [Rust](https://www.rust-lang.org) along with the `cargo` crate manager for Rust.

### From Cargo
Coming soon.

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
### MSBWT Building
The Multi-String Burrows Wheeler Transform (MSBWT or BWT) must be build prior to read correction. FMLRC v2 uses the same data structure as v1, so the [original instructions](https://github.com/holtjma/fmlrc#building-the-short-read-bwt) can be used.

Given a FASTQ file of accurate reads (`reads.fq.gz`), you can also use the following command from this crate to create a BWT at `comp_msbwt.npy`.  Note that this command requires the [ropebwt2](https://github.com/lh3/ropebwt2) executable to be installed:
```
gunzip -c reads.fq.gz | \
    awk 'NR % 4 == 2' | \
    sort | \
    tr NT TN | \
    ropebwt2 -LR | \
    tr NT TN | \
    fmlrc2-convert comp_msbwt.npy
```

#### Optional Construction Speedup
If you are **only** using the BWT for correction, then the `sort` can be removed from the above command. This will reduce construction time significantly, but loses the read recovery property of the BWT.

### Correction
Assuming the accurate-read BWT is built (`comp_msbwt.npy`) and uncorrected reads are available (fastq/fastq, gzip optional, `uncorrected.fq.gz`), invoking FMLRC v2 is fairly simple:
```
fmlrc2 [OPTIONS] <comp_msbwt.npy> <uncorrected.fq.gz> <corrected_reads.fa>
```

Currently, only uncompressed FASTA is supported for output reads.

#### Options to consider
1. `-h` - see full list of options and exit
2. `-k`, `-K` - sets the k-mer sizes to use, default is `[21, 59]`; all values are sorted from lowest to highest prior to correction
3. `-t`, `--threads` - number of correction threads to use (default: 1)
4. `-C`, `--cache_size` - the length of sequences to pre-compute (i.e. `C`-mers); will reduce CPU-time of queries by `O(C)` but *increases* cache memory usage by `O(6^C)`; default of `8` uses ~25MB; if memory is not an issue, consider using `10` with ~1GB cache footprint (or larger if memory _really_ isn't an issue)

## FMLRC v2 core differences
1. Implemented in [Rust](https://www.rust-lang.org) instead of C++ - this comes will all the benefits of Rust including `cargo`, such as easy installation of the binary and supporting structs/functions along with documentation
2. Unlimited `k`/`K` parameters - FMLRC v1 allowed 1 or 2 sizes for `k` only; FMLRC v2 can have the option set as many times as desired at increased CPU time (for example, a 3-pass correction with `k=[21, 59, 79]`) 
3. Call caching - FMLRC v2 pre-computes all _k_-mers of a given size. This reduces the run-time significantly by cutting reducing calls to the FM-index.
4. Input handling - thanks to [needletail](https://crates.io/crates/needletail), the uncorrected reads can be in FASTA/FASTQ and may or may not be gzip compressed.
5. SIMD accelerated alignment - thanks to [triple_accel](https://crates.io/crates/triple_accel), the correction alignment step can be accelerated with SIMD instructions when available
6. Unit testing - FMLRC v2 has unit testing through the standard Rust testing framework (i.e. `cargo test`)

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