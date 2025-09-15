# Experimental Evalutation of InsPIRe
The following document outlines the steps to reproduce all results from the paper.
This document assume the project has been correctly built using the instructions in the parent directory.
Run `cargo run --release --bin inspire -- --num-items 4096 --item-size-bits 8192 --dim0 1024` to ensure correct build the project before proceeding.

For reproducing tables from the paper you will need to install the necessary python requirements.
```
    virtualenv .venv
    . .venv/bin/activate
    pip install -r requirements.txt
```

## PIR as a function of DB size

To generate Table 3 and Table 4 from the paper, run the following:
```
run-pir-experiment-small.sh
run-pir-experiment-medium.sh
run-pir-experiment-large.sh
run-inspire-squared.sh
run-inspire.sh
```

Then, to visualize all results, run the following scripts:

To produce (an extended version of) Table 2 from the conference version or Tables 3,4, and 5 from the full version, run the following:
```
python visualize_all.py
```

To generate the appendix tables (from the full version), run the following:
```
python visualize_pir.py
python visualize_inspire_squared.py
python visualize_inspire.py
```

## Packing time as a function of number of packed items, for various values of gamma

Run the following to generate the requires logs
```
cargo run --release --bin pack
```
then run the following to generate plots
```
python visualize_packing.py
```


## Related Work

### Hintless PIR

```
git clone git@github.com:RasoulAM/hintless_pir-clone.git
cd hintless_pir-clone
sudo apt install bazel-bootstrap
sudo curl -L https://github.com/bazelbuild/bazelisk/releases/latest/download/bazelisk-linux-amd64 -o /usr/local/bin/bazel
sudo chmod +x /usr/local/bin/bazel
export USE_BAZEL_VERSION=7.1.2
bazel test //...
bazel build //hintless_simplepir:hintless_simplepir_benchmarks
cd ..
./run-hintlesspir.sh
python parse_hintlesspir.py
```