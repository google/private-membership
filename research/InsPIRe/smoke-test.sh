OUT_FILE="smoke-test.txt"


echo Start > $OUT_FILE &&

# # InsPIRe
cargo run --release --bin inspire -- --num-items $((2**12)) --item-size-bits $((1024*8)) --dim0 1024 >> $OUT_FILE &&
cargo run --release --bin inspire -- --num-items $((2**15)) --item-size-bits $((32*1024*8)) --dim0 8192 >> $OUT_FILE &&


# SimplePIR + InspiRING
cargo run --release --bin run -- --num-items 4096 --item-size-bits $((6 * 1024 * 8)) --protocol-type SimplePIR --second-level-packing-mask InspiRING --second-level-packing-body NoPacking --gammas 2048 >> $OUT_FILE &&
cargo run --release --bin run -- --num-items 2049 --item-size-bits 8192 --protocol-type SimplePIR --second-level-packing-mask InspiRING --second-level-packing-body NoPacking --gammas 1024 >> $OUT_FILE &&
cargo run --release --bin run -- --num-items 2048 --item-size-bits $((32*1024*8)) --protocol-type SimplePIR --second-level-packing-mask InspiRING --second-level-packing-body NoPacking --gammas 2048 >> $OUT_FILE &&
cargo run --release --bin run -- --num-items 2049 --item-size-bits $((32*1024*8)) --protocol-type SimplePIR --second-level-packing-mask InspiRING --second-level-packing-body NoPacking --gammas 2048 >> $OUT_FILE &&
cargo run --release --bin run -- --num-items 256 --item-size-bits $((3*1024*4)) --protocol-type SimplePIR --second-level-packing-mask InspiRING --second-level-packing-body NoPacking --gammas 1024 >> $OUT_FILE &&

# InsPIRe_0
cargo run --release --bin run -- --num-items 2048 --item-size-bits 8 --protocol-type DoublePIR --second-level-packing-mask InspiRING --second-level-packing-body InspiRING --gammas 2048 >> $OUT_FILE &&
cargo run --release --bin run -- --num-items 2048 --item-size-bits 8 --protocol-type DoublePIR --second-level-packing-mask InspiRING --second-level-packing-body InspiRING --gammas 1024 >> $OUT_FILE &&
cargo run --release --bin run -- --num-items 2048 --item-size-bits 8 --protocol-type DoublePIR --second-level-packing-mask InspiRING --second-level-packing-body InspiRING --gammas 512 >> $OUT_FILE &&

# InsPIRe_0 without packing part of the second hint
cargo run --release --bin run -- --num-items 2048 --item-size-bits 8 --protocol-type DoublePIR --second-level-packing-mask InspiRING --second-level-packing-body NoPacking --gammas 2048 >> $OUT_FILE &&
cargo run --release --bin run -- --num-items 2048 --item-size-bits 8 --protocol-type DoublePIR --second-level-packing-mask InspiRING --second-level-packing-body NoPacking --gammas 1024 >> $OUT_FILE &&
cargo run --release --bin run -- --num-items 2048 --item-size-bits 8 --protocol-type DoublePIR --second-level-packing-mask InspiRING --second-level-packing-body NoPacking --gammas 512 >> $OUT_FILE &&

# InsPIRe^(2)
cargo run --release --bin run -- --num-items 2048 --item-size-bits 2048 --protocol-type InsPIRe --second-level-packing-mask InspiRING --second-level-packing-body InspiRING --gammas 128,1024,128 >> $OUT_FILE &&
cargo run --release --bin run -- --num-items 2048 --item-size-bits 2048 --protocol-type InsPIRe --second-level-packing-mask InspiRING --second-level-packing-body InspiRING --gammas 256,1024,128 >> $OUT_FILE &&
cargo run --release --bin run -- --num-items $((2**26)) --item-size-bits 1 --protocol-type InsPIRe --second-level-packing-mask InspiRING --second-level-packing-body InspiRING --gammas 16,1024,1024 >> $OUT_FILE &&

# DoublePIR
cargo run --release --bin run -- --num-items 2048 --item-size-bits 8 --protocol-type DoublePIR --second-level-packing-mask CDKS --second-level-packing-body CDKS --gammas 2048 >> $OUT_FILE &&

# Other
cargo run --release --bin run -- --num-items 2048 --item-size-bits 2048 --protocol-type InsPIRe --second-level-packing-mask InspiRING --second-level-packing-body NoPacking --gammas 16,1024,16 >> $OUT_FILE &&

#####

# TODO: not working at the moment
# cargo run --release --bin run -- --num-items 2048 --item-size-bits $((32*1024*14)) --protocol-type SimplePIR --second-level-packing-mask CDKS --second-level-packing-body NoPacking --gammas 2048 >> $OUT_FILE &&
# cargo run --release --bin run -- --num-items 10000 --item-size-bits 8 --protocol-type DoublePIR --second-level-packing-mask CDKS --second-level-packing-body NoPacking >> $OUT_FILE &&

echo Done