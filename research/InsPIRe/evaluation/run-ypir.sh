#!/bin/bash

TARGET=ypir/target/release/run

RESULTS_DIR="results/simpleypir"
mkdir -p $RESULTS_DIR

# Define configurations
configs=(
  "12 $((32 * 1024 * 8))" # log_db_num_items = 12, item_size_bits = 32 KB

  # ## Uncomment the following lines to run full experiments
  # "15 $((32 * 1024 * 8))" # log_db_num_items = 15, item_size_bits = 32 KB
  # "18 $((32 * 1024 * 8))" # log_db_num_items = 18, item_size_bits = 32 KB
)

# Iterate over configurations
for config in "${configs[@]}"; do
  read log_db_num_items item_size_bits <<< "$config"
  num_items=$((2**log_db_num_items))

  output_file="$RESULTS_DIR/${log_db_num_items}.txt"

  $TARGET $num_items $item_size_bits 1 5 $output_file --is-simplepir 

  num_items=$((2**log_db_num_items))
  db_mb=$((num_items*item_size_bits / 8 / 1024 / 1024))

  sed -i "2i\
    \"specs\": {\
      \"inputDatabaseSizeMb\": $db_mb\
    }," $output_file

done

RESULTS_DIR="results/ypir"
mkdir -p $RESULTS_DIR

# ## Uncomment the following lines to run full experiments
# for db_mb in 256 1024 8192 32768; do
for db_mb in 256 1024; do
  echo "DB Size=${db_mb}"
  num_items=$((2**20 * 8 * $db_mb)) # Compute 2^log_num_items
  item_size_bits=1 

  # doublepir variant
  output_file="$RESULTS_DIR/db=${db_mb}MB-DoublePIR.txt"
  $TARGET $num_items $item_size_bits 1 5 $output_file 

  sed -i "2i\
    \"specs\": {\
      \"inputDatabaseSizeMb\": $db_mb\
    }," $output_file

done
