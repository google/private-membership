RESULTS_DIR="results/pir-experiments/large"
mkdir -p $RESULTS_DIR

mkdir -p saved_targets
cp ../target/release/run saved_targets
TARGET=saved_targets/run

# Define configurations
configs=(
  # "3 $((32 * 1024 * 8))"  # log_db_num_items = 3, item_size_bits = 32 KB
  "12 $((32 * 1024 * 8))" # log_db_num_items = 12, item_size_bits = 32 KB
  "15 $((32 * 1024 * 8))" # log_db_num_items = 15, item_size_bits = 32 KB
  "18 $((32 * 1024 * 8))" # log_db_num_items = 18, item_size_bits = 32 KB
  "19 $((64 * 1024 * 8))" # log_db_num_items = 19, item_size_bits = 64 KB
)

# Iterate over configurations
for config in "${configs[@]}"; do
  read log_db_num_items item_size_bits <<< "$config"
  num_items=$((2**log_db_num_items))

  gamma=1024
  output_file="$RESULTS_DIR/${log_db_num_items}-SimplePIR-InspiRING-gamma=${gamma}.json"

  $TARGET --protocol-type SimplePIR \
    --num-items $num_items \
    --item-size-bits $item_size_bits \
    --second-level-packing-mask InspiRING \
    --second-level-packing-body NoPacking \
    --out-report-json $output_file \
    --gammas $gamma \
    --trials 5

  gamma=2048
  output_file="$RESULTS_DIR/${log_db_num_items}-SimplePIR-InspiRING-gamma=${gamma}.json"
  $TARGET --protocol-type SimplePIR \
    --num-items $num_items \
    --item-size-bits $item_size_bits \
    --second-level-packing-mask InspiRING \
    --second-level-packing-body NoPacking \
    --out-report-json $output_file \
    --gammas $gamma \
    --trials 5

  output_file="$RESULTS_DIR/${log_db_num_items}-SimplePIR-CDKS.json"

  gamma=2048
  $TARGET --protocol-type SimplePIR \
    --num-items $num_items \
    --item-size-bits $item_size_bits \
    --second-level-packing-mask CDKS \
    --second-level-packing-body NoPacking \
    --gammas $gamma \
    --out-report-json $output_file \
    --trials 5

done
