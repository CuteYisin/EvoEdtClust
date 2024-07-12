#!/bin/bash

sim=0.9

INPUT_FASTA=./testdata/testdata50-10.fa
OUTPUT_DIR=./output/testdata50-10

if [ ! -d ${OUTPUT_DIR} ]; then
    mkdir -p ${OUTPUT_DIR}
fi

make clean
make -j

/usr/bin/time --format="%e\t%M" -a -o ./_TimeMemory.txt ./evoEdtClust.exe -t ${sim} -p 1 ${INPUT_FASTA} ${OUTPUT_DIR} > log.txt
