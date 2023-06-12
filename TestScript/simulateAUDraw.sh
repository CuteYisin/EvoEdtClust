#!/bin/bash
WORK_DIR=/data/cabins/yxxiang/simulate
PROJECT_DIR=/data/cabins/yxxiang/HomoPic/project

number_alg=2
label=Linclust,ALFATClust
sim_threshold=0.9

if [ ! -d ${WORK_DIR}/result/pic ];then
        mkdir ${WORK_DIR}/result/pic
fi

for folder in ${WORK_DIR}/result/*; do
	if [ -d "$folder" ] && [ -d "${folder}/groundTruth" ]; then
		echo "=== Folder $folder contains groundTruth, plot ROC and PR curve."
		
		fasta_file=${WORK_DIR}/$(basename ${folder}).fasta
		if [ ! -f ${fasta_file} ]; then
			fasta_file=${WORK_DIR}/$(basename ${folder}).fa
		fi
		
		number_seq=$(grep -c ">" ${fasta_file})

		groundtruth_file=${folder}/groundTruth/groundTruth_${sim_threshold}.txt
		input_file=${folder}/Linclust/LinSparseMatrix_,${folder}/ALFATClust/ALFATSparseMatrix_
		output_pic=${WORK_DIR}/result/pic/$(basename ${folder})
		python AUDraw.py ${number_alg} ${number_seq} ${groundtruth_file} ${label} ${input_file} ${output_pic} #${LOG_FILE}
	fi
done
