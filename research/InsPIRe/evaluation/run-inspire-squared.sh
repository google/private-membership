#!/bin/bash

RESULTS_DIR="results/inspire-squared"
mkdir -p "$RESULTS_DIR"

mkdir -p saved_targets
cp ../target/release/run saved_targets/run
TARGET=saved_targets/run

echo "Starting"

# Uncomment this line to perform full experiments
# for log_db_bits in 31 33 36 38; do
for log_db_bits in 33; do
  db_bits=$((2**log_db_bits))
  item_size_bits=1
  for gamma_0 in 16 32 64; do
    gamma_1=1024
    gamma_2=$gamma_0
    for factor in 4 2 1; do
        output_file="$RESULTS_DIR/${db_bits}-${item_size_bits}-InsPIRe-gamma=${gamma_0}-${gamma_1}-${gamma_2}-factor=${factor}.json"
        $TARGET --protocol-type InsPIRe \
          --num-items $db_bits \
          --item-size-bits 1 \
          --second-level-packing-mask InspiRING \
          --second-level-packing-body InspiRING \
          --performance-factor $factor \
          --gammas $gamma_0,$gamma_1,$gamma_2 \
          --out-report-json $output_file \
          --trials 5 >> $RESULTS_DIR/log.txt
    done
  done
done