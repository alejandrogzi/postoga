#!/usr/bin/env bash 

LOG="./supply/test/postoga.log"
MODE="base"
DIR="./supply/test"


if [[ -d $DIR ]]; then
    ./postoga.py $MODE --path $DIR --to gff -th 0.5
else
    echo "Directory $DIR not found, clone the repository again or contact the developer"
fi

head -n 9 $LOG