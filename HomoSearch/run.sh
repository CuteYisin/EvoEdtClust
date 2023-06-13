#!/bin/bash

BATCH=50-1

#INPUT_FASTA=./testdata/testdata${BATCH}.fa
INPUT_FASTA=./testdata/metaFamily_8339.fasta
OUTPUT_DIR=./output

make clean
make -j


#rm -rf ${OUTPUT_DIR}
#mkdir ${OUTPUT_DIR}

./hierarchical_lsh -m $1 -t $2 ${INPUT_FASTA} ${OUTPUT_DIR}
