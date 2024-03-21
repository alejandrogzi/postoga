#!/usr/bin/env bash 

LOG="./test_out/postoga.log"
MODE="base"
TESTDIR="./supply/test"
OUTDIR="./test_out"


if [[ -d $TESTDIR ]]; then
    ./postoga.py $MODE --td $TESTDIR -o $OUTDIR --to gff -th 0.5
else
    echo "Directory $TESTDIR not found, clone the repository again or contact the developer"
fi

head -n 9 $LOG