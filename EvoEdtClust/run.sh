#!/bin/bash

sim=0.8

#INPUT_FASTA=./testdata/mod_testdata50-10.fa
#OUTPUT_DIR=./output/mod_testdata

#INPUT_FASTA=/data/cabins/yxxiang/simulate/testdata50-10.fa
#OUTPUT_DIR=./output/testdata50-10
#groundtruth_file=/data/cabins/yxxiang/simulate/result/testdata50-10/groundTruth/groundTruth_0.8.txt

#INPUT_FASTA=/data/cabins/yxxiang/simulate/testdata100-100.fa
#OUTPUT_DIR=./output/testdata100-100

INPUT_FASTA=/data/cabins/yxxiang/real/metaFamily_8339.fasta
OUTPUT_DIR=./output/metaFamily_8339
GDDIR=/data/cabins/yxxiang/real/result/metaFamily_8339

#INPUT_FASTA=/data/cabins/yxxiang/real/Homo_sapiens_human_37744.fasta
#OUTPUT_DIR=./output/Homo_sapiens_human_37744
#GDDIR=/data/cabins/yxxiang/real/result/Homo_sapiens_human_37744

#INPUT_FASTA=/data/cabins/yxxiang/real/MERC_10000000.fasta
#OUTPUT_DIR=./output/MERC_10000000

make clean
make -j

rm -rf ${OUTPUT_DIR}
mkdir ${OUTPUT_DIR}

./evoEdtClust.exe -t ${sim} ${INPUT_FASTA} ${OUTPUT_DIR}
python ./evaluation/onlytestEvo.py \
    ${GDDIR}/index.lookup \
    ${GDDIR}/groundTruth/groundTruth_${sim}.txt \
    ${OUTPUT_DIR}/EvoEdtClust_${sim}.out \
    ${GDDIR}/Linclust/LinSparseMatrix_${sim}0.txt \
    ${GDDIR}/ALFATClust/ALFATSparseMatrix_${sim}0.txt

#/usr/bin/time --format="%e\t%M" -a -o ./_TimeMemory.txt ./evoEdtClust.exe -t 0.9 ${INPUT_FASTA} ${OUTPUT_DIR} > ${OUTPUT_DIR}/log.txt
