#!/bin/bash

RESULTS_DIR="results/inspire_0"
mkdir -p $RESULTS_DIR

mkdir -p saved_targets
cp ../target/release/run saved_targets/
TARGET=saved_targets/run

# Uncomment this line to perform full experiments
# for db_mb in 1024 8192 32768; do
for db_mb in 1024; do
  echo "DB Size=${db_mb}"
  num_items=$((2**20 * 8 * $db_mb)) # Compute 2^log_num_items
  item_size_bits=1 

  # doublepir variant
  for pair in "InspiRING NoPacking 1024" "InspiRING NoPacking 2048" "CDKS CDKS 2048"; do
    read -r second_level_packing_mask second_level_packing_body gamma <<< "$pair"
    output_file="$RESULTS_DIR/db=${db_mb}MB-DoublePIR-${second_level_packing_mask}-${second_level_packing_body}-gamma=${gamma}.json"
    $TARGET --protocol-type DoublePIR \
      --num-items $num_items \
      --item-size-bits $item_size_bits \
      --second-level-packing-mask $second_level_packing_mask \
      --second-level-packing-body $second_level_packing_body \
      --gammas $gamma \
      --out-report-json $output_file \
      --trials 5
  done
done
