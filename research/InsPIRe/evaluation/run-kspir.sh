#!/bin/bash

RESULTS_DIR="results/kspir"
mkdir -p $RESULTS_DIR

TARGET=kspir-clone/build/tests/test-pir

echo Start &&
$TARGET 256 > $RESULTS_DIR/256.txt &&
echo =============================== &&
$TARGET 512 > $RESULTS_DIR/512.txt &&
echo =============================== &&
$TARGET 1024 > $RESULTS_DIR/1024.txt &&
echo =============================== &&

### Uncomment the following lines to run full experiments
# $TARGET 2048 > $RESULTS_DIR/2048.txt &&
# echo =============================== &&
# $TARGET 4096 > $RESULTS_DIR/4096.txt &&
# echo =============================== &&
# $TARGET 8192 > $RESULTS_DIR/8192.txt &&
# echo =============================== &&
# $TARGET 16384 > $RESULTS_DIR/16384.txt &&
# echo =============================== &&
# $TARGET 32768 > $RESULTS_DIR/32768.txt &&
echo Done