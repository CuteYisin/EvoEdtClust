#!/bin/bash

BATCH=50-1
#BATCH=50-50

#INPUT_FASTA=/data/lobby/lsh_simul/testdata${BATCH}.fa
INPUT_FASTA=/data/cabins/yxxiang/stable_lsh/testdata/testdata58.fa
OUTPUT_DIR=./output

make clean
make -j


#rm -rf ${OUTPUT_DIR}
mkdir ${OUTPUT_DIR}/$1-$2

./distance_space_embedding -w $1 -k $2 ${INPUT_FASTA} ${OUTPUT_DIR}
