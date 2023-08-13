#!/bin/bash

sim=0.9

INPUT_FASTA=/data/cabins/yxxiang/simulate/testdata50-10.fa
OUTPUT_DIR=./output/testdata50-10

#INPUT_FASTA=/data/cabins/yxxiang/simulate/testdata50-1.fa
#OUTPUT_DIR=./output/testdata50-1
#groundtruth_file=/data/cabins/yxxiang/simulate/result/testdata50-10/groundTruth/groundTruth_0.8.txt

#INPUT_FASTA=/data/cabins/yxxiang/simulate/testdata100-100.fa
#OUTPUT_DIR=./output/testdata100-100

#INPUT_FASTA=/data/cabins/yxxiang/real/metaFamily_8339.fasta
#OUTPUT_DIR=./output/metaFamily_8339

#INPUT_FASTA=/data/cabins/yxxiang/real/Homo_sapiens_human_37744.fasta
#OUTPUT_DIR=./output/Homo_sapiens_human_37744
#GDDIR=/data/cabins/yxxiang/real/result/Homo_sapiens_human_37744

#INPUT_FASTA=/data/cabins/yxxiang/real/MERC_1000000.fasta
#OUTPUT_DIR=./output/MERC_1000000

#INPUT_FASTA=/data/cabins/yxxiang/MERC/MERC_300000.fasta
#OUTPUT_DIR=./output/MERC_300000

if [ ! -d ${OUTPUT_DIR} ]; then
    mkdir -p ${OUTPUT_DIR}
fi

make clean
make -j

/usr/bin/time --format="%e\t%M" -a -o ./_TimeMemory.txt ./evoEdtClust.exe -t ${sim} -p 1 ${INPUT_FASTA} ${OUTPUT_DIR} > log.txt
#echo "== if there's no diff, there's no content between here"
#diff ${OUTPUT_DIR}/reserved_EvoEdtClust_0.9.out ${OUTPUT_DIR}/EvoEdtClust_0.9.out
#echo "== and here"
