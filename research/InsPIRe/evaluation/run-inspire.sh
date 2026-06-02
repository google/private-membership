#!/bin/bash

RESULTS_DIR="results/inspire"
mkdir -p $RESULTS_DIR

mkdir -p saved_targets
cp ../target/release/inspire saved_targets/inspire
TARGET=saved_targets/inspire


echo "PIR over a database consisting of 1-bit entries for three database sizes, 1 GB, 8 GB, and 32 GB."

# Define combinations where the first element is log2_num_items
# and the second is a space-separated list of log_dim0 values
declare -A combinations_map
combinations_map[33]="17 16 15 14 13 12 11 10"

# # Uncomment this line to perform full experiments
# combinations_map[36]="17 16 15 14 13 12"
# combinations_map[38]="20 19 18 17 16 15"

item_size_bits=1
for log2_num_items in "${!combinations_map[@]}"; do
  log_dim0_list="${combinations_map[$log2_num_items]}"

  for log_dim0 in $log_dim0_list; do
    echo $log_dim0,$log2_num_items
    num_items=$((2**log2_num_items))
    dim0=$((2**log_dim0))

    output_file="$RESULTS_DIR/${num_items}-${item_size_bits}-dim0=${dim0}.json"
    $TARGET --num-items $num_items \
      --item-size-bits $item_size_bits \
      --dim0 $dim0 \
      --out-report-json $output_file \
      --trials 5 >> $RESULTS_DIR/log.txt \
      --label "one-bit-payload"

  done
done


echo "PIR over a database consisting of 64B entries for three database sizes, 1 GB, 8 GB, and 32 GB."

# Define combinations where the first element is log2_num_items
# and the second is a space-separated list of log_dim0 values
declare -A combinations_map
combinations_map[33]="17 16 15 14 13 12"

# # Uncomment this line to perform full experiments
# combinations_map[36]="19 18 17 16 15 14 13"
# combinations_map[38]="20 19 18 17 16 15"

item_size_bits=$(((2**9)))
for log2_num_bits in "${!combinations_map[@]}"; do
  log_dim0_list="${combinations_map[$log2_num_bits]}"
  log2_num_items=$((log2_num_bits-9))

  for log_dim0 in $log_dim0_list; do

    echo $log_dim0,$log2_num_items
    num_items=$((2**log2_num_items))
    dim0=$((2**log_dim0))
    
    output_file="$RESULTS_DIR/${num_items}-${item_size_bits}-dim0=${dim0}.json"
    $TARGET --num-items $num_items \
      --item-size-bits $item_size_bits \
      --dim0 $dim0 \
      --out-report-json $output_file \
      --trials 5 >> $RESULTS_DIR/log.txt \
      --label "512-bit-payload"

  done
done

echo "PIR over a database consisting of 32 KB entries for three database sizes, 1 GB, 8 GB, and 32 GB."

# Define combinations where the first element is log2_num_items
# and the second is a space-separated list of log_dim0 values
declare -A combinations_map
combinations_map[33]="15 14 13 12 11 10 9 8"

# # Uncomment this line to perform full experiments
# combinations_map[36]="18 17 16 15 14 13 12"
# combinations_map[38]="20 19 18 17 16 15 14"

item_size_bits=$((32*1024*8))
for log2_num_bits in "${!combinations_map[@]}"; do
  log_dim0_list="${combinations_map[$log2_num_bits]}"
  log2_num_items=$((log2_num_bits-18))

  for log_dim0 in $log_dim0_list; do

    echo $log_dim0,$log2_num_items
    num_items=$((2**log2_num_items))
    dim0=$((2**log_dim0))
    
    output_file="$RESULTS_DIR/${num_items}-${item_size_bits}-dim0=${dim0}.json"
    $TARGET --num-items $num_items \
      --item-size-bits $item_size_bits \
      --dim0 $dim0 \
      --out-report-json $output_file \
      --trials 5 >> $RESULTS_DIR/log.txt \
      --label "32KB-payload"

  done
done