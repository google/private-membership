# InsPIRe: Communication-Efficient PIR with Silent Preprocessing

This is an implementation of the InsPIRe protocol for Private Information Retreival (PIR).
This scheme was introduced in the paper ["InsPIRe: Communication-Efficient PIR with Silent Preprocessing"](https://eprint.iacr.org/2025/1352)

## Reproducing Experimental Evaluation
To reproduce the experimental results shown in the paper, refer to the `evaluation` directory.

## Running

To build and run this code:
1. Ensure you are running on Ubuntu, and that AVX-512 is available on the CPU (you can run `lscpu` and look for the `avx512f` flag).
<!-- Our benchmarks were collected using `TODO: MACHINE`, which has all necessary CPU features. -->
2. Run `sudo apt-get update && sudo apt-get install -y build-essential`.
3. [Install Rust using rustup](https://www.rust-lang.org/tools/install) using `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`.
  - Select `1) Proceed with installation (default)` when prompted
  - After installation, configure the current shell as instructed by running `source "$HOME/.cargo/env"`
4. Clone the repository and run `cd InsPIRe`.
5. Run `cargo run --release --bin inspire -- --num-items 4096 --item-size-bits 8192 --dim0 1024` to run the code on a random 4 MB database to build the project.
The first time you run this command, Cargo will download and install the necessary libraries to build the code;
later calls will not take as long. Stability warnings can be safely ignored. 

### Options
To pass arguments, make sure to run `cargo run --release -- <ARGS>` (the ` -- ` is important).
Passing `--verbose` or setting the environment variable `RUST_LOG=debug`
will enable detailed logging. PIR results are checked for correctness, unless the `--online-only` flag is used.
The full command-line parameters are as follows:

```
Usage: cargo run --release -- [OPTIONS] <NUM_ITEMS> [ITEM_SIZE_BITS] [TRIALS] [OUT_REPORT_JSON]

Arguments:
  <NUM_ITEMS>        Number of items in the database
  [ITEM_SIZE_BITS]   Size of each item in bits (optional, default 1)
  [PROTOCOL_TYPE]    Protocol to use, (default DoublePIR w/ InspiRING)
  [SECOND_LEVEL_PACKING_MASK] Packing type in the first layer (SimplePIR & InsPIRe) and second layer over the mask (DoublePIR & InsPIRE)
  [SECOND_LEVEL_PACKING_BODY] Packing type for body in second layer (DoublePIR & InsPIRE)
  [TRIALS]           Number of trials (optional, default 5) to run the YPIR scheme 
                     and average performance measurements over (with one additional warmup trial excluded)
  [OUT_REPORT_JSON]  Output report file (optional) where results will be written in JSON

Options:
  -v, --verbose  Verbose mode (optional) if set, the program will print debug logs to stderr
  -h, --help     Print help
  -V, --version  Print version
```

### Setting the arguements
Below are template examples for all possible ways of running the protocol

```
# SimplePIR + InpiRING packing
cargo run --release --bin run -- --num-items $NUM_ITEMS --item-size-bits $ITEM_SIZE_IN_BITS --protocol-type SimplePIR --second-level-packing-mask InspiRING --second-level-packing-body NoPacking --gammas $gamma_0 

# InsPIRe_0
cargo run --release --bin run -- --num-items $NUM_ITEMS --item-size-bits $ITEM_SIZE_IN_BITS --protocol-type DoublePIR --second-level-packing-mask InspiRING --second-level-packing-body InspiRING --gammas $gamma_0,$gamma_1,$gamma_2

# InsPIRe^(2)
cargo run --release --bin run -- --num-items $NUM_ITEMS --item-size-bits $ITEM_SIZE_IN_BITS --protocol-type InsPIRe --second-level-packing-mask InspiRING --second-level-packing-body InspiRING --gammas $gamma_0,$gamma_1,$gamma_2

# InsPIRe 
cargo run --release --bin inspire -- --num-items $NUM_ITEMS --item-size-bits $ITEM_SIZE_IN_BITS --dim0 $dim_0

```