#!/bin/bash
WORK_DIR=/data/cabins/yxxiang/real                  # The folder of fasta file
PROJECT_DIR=/data/cabins/yxxiang/HomoPic/project        # The direction of your project saved
sim_threshold=0.9

# Compiling groundTruth
make -f ${PROJECT_DIR}/groundTruth/makefile clean >> ${PROJECT_DIR}/groundTruth/Compilation.log
make -f ${PROJECT_DIR}/groundTruth/makefile -j >> ${PROJECT_DIR}/groundTruth/Compilation.log

# Generating ground truth only for sequence sets with 500,000 sequences and less
find ${WORK_DIR} -maxdepth 1 \( -name "*.fa" -o -name "*.fasta" \) | while read file; do
	if [ "$(wc -l < "$file")" -le 100000 ]; then
		filename=$(basename $file)
		echo "=== Ground truth for ${filename}  is being generated..."
		groundtruth_dir=${WORK_DIR}/result/${filename%.*}/groundTruth
		mkdir -p ${groundtruth_dir}
		echo -n -e ${sim_threshold}"\t" >> ${groundtruth_dir}/_TimeMemory.txt
		/usr/bin/time --format="%e\t%M" -a -o ${groundtruth_dir}/_TimeMemory.txt ${PROJECT_DIR}/groundTruth/groundTruth -i ${file} -o ${groundtruth_dir}/groundTruth_${sim_threshold}.txt -t ${sim_threshold} > ${groundtruth_dir}/log
	fi
done
