#!/bin/bash

OUR_DIR="results/application-experiments"

mkdir $OUR_DIR &&

############## Peer routing
cargo run --release --bin run -- --num-items 256 --item-size-bits $((3 * 1024 * 4)) --protocol-type SimplePIR --second-level-packing-mask InspiRING --second-level-packing-body NoPacking --small-params --gammas 1024 $ONLINE_ONLY --trials 10 --out-report-json $OUR_DIR/peer-routing.json &&

############## Provider Advertisements
bins_list=(4096 8192 16384)
for bins in "${bins_list[@]}"; do
    item_size_bits=$((200000 * 8 * 110 / $bins))
    cargo run --release --bin run -- --num-items $bins --item-size-bits "${item_size_bits}" --protocol-type SimplePIR --second-level-packing-mask InspiRING --second-level-packing-body NoPacking --gammas 2048 $ONLINE_ONLY --trials 5 --out-report-json $OUR_DIR/prov-adv-simplepir-bins=${bins}.json
    cargo run --release --bin inspire -- --num-items $bins --item-size-bits "${item_size_bits}" --trials 5 --dim0 4096 --out-report-json $OUR_DIR/prov-adv-inspire-bins=${bins}.json
done

############## Bitswap
# Define the dimensions to iterate through
dims=(32768 16384 8192 4096)
NUM_ITEMS=$((2**14))
ITEM_SIZE_BITS=$((256*1024*8))

for dim0 in "${dims[@]}"; do
    cargo run --release --bin inspire -- --num-items ${NUM_ITEMS} --item-size-bits ${ITEM_SIZE_BITS} --dim0 ${dim0} --trials 5 --out-report-json $OUR_DIR/bitswap-${dim0}.json 
done


############## Device Enrollment
NUM_ITEMS_list=($((2**24)) $((2**25)) $((2**26)) $((2**27)))
ITEM_SIZE_BITS=$((64*8))
dim0=$((2**14))
for NUM_ITEMS in "${NUM_ITEMS_list[@]}"; do
    cargo run --release --bin inspire -- --num-items ${NUM_ITEMS} --item-size-bits ${ITEM_SIZE_BITS} --dim0 ${dim0} --trials 5 --out-report-json $OUR_DIR/device-enrollment-num=${NUM_ITEMS}-dim=${dim0}.json 
    cargo run --release --bin run -- --num-items ${NUM_ITEMS} --item-size-bits ${ITEM_SIZE_BITS} --protocol-type SimplePIR --second-level-packing-mask CDKS --second-level-packing-body NoPacking --gammas 2048 --trials 5 --out-report-json $OUR_DIR/device-enrollment-num=${NUM_ITEMS}-simpleypir.json 
done

echo Done