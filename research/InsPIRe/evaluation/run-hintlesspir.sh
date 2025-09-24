#!/bin/bash
RESULTS_DIR="results/hintlesspir"
mkdir -p $RESULTS_DIR

./hintless_pir-clone/bazel-bin/hintless_simplepir/hintless_simplepir_benchmarks --benchmark_filter=BM_HintlessPir --benchmark_time_unit=ms --num_rows=$((2**10)) --num_cols=$((2**10)) > $RESULTS_DIR/log-10-10.txt
./hintless_pir-clone/bazel-bin/hintless_simplepir/hintless_simplepir_benchmarks --benchmark_filter=BM_HintlessPir --benchmark_time_unit=ms --num_rows=$((2**15)) --num_cols=$((2**15)) > $RESULTS_DIR/log-15-15.txt

# # Uncomment this line to perform full experiments
# ./hintless_pir-clone/bazel-bin/hintless_simplepir/hintless_simplepir_benchmarks --benchmark_filter=BM_HintlessPir --benchmark_time_unit=ms --num_rows=$((2**16)) --num_cols=$((2**17)) > $RESULTS_DIR/log-16-17.txt
# ./hintless_pir-clone/bazel-bin/hintless_simplepir/hintless_simplepir_benchmarks --benchmark_filter=BM_HintlessPir --benchmark_time_unit=ms --num_rows=$((2**17)) --num_cols=$((2**18)) > $RESULTS_DIR/log-17-18.txt
