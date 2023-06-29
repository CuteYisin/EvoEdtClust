#!/bin/bash

BATCH=50-10

#INPUT_FASTA=./testdata/testdata${BATCH}.fa
INPUT_FASTA=./testdata/mod_testdata${BATCH}.fa
#INPUT_FASTA=./testdata/metaFamily_8339.fasta
#INPUT_FASTA=/data/cabins/yxxiang/real/MERC_10000000.fasta
OUTPUT_DIR=./output

make clean
make -j

rm -rf ${OUTPUT_DIR}
mkdir ${OUTPUT_DIR}

./hierarchical_lsh.exe -m 1 -t 0.9 ${INPUT_FASTA} ${OUTPUT_DIR} 
#/usr/bin/time --format="%e\t%M" -a -o ./_TimeMemory.txt ./hierarchical_lsh.exe -m 1 -t 0.9 ${INPUT_FASTA} ${OUTPUT_DIR} > ${OUTPUT_DIR}/log.txt
