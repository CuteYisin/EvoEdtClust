#!/bin/bash
WORK_DIR=/data/cabins/yxxiang/simulate
PROJECT_DIR=/data/cabins/yxxiang/HomoPic/project

number_alg=9
sim_threshold_Linclust=0.90
sim_threshold_ALFATClust=0.80
label=Linclust_50,Linclust_100,Linclust_200,ALFATClust_50,ALFATClust_100,ALFATClust_200,BruteForce_50,BruteForce_100,BruteForce_200
picture_file=${WORK_DIR}/result/pic
time_memory_file=""

for l in $(echo $label | tr ',' ' '); do
  time_memory_file="${time_memory_file},${picture_file}/timeMemory${l}.txt"
done
time_memory_file=${time_memory_file#","}

for seq_length in 50 100 200;do
	> "${picture_file}/timeMemoryBruteForce_${seq_length}.txt"
	> "${picture_file}/timeMemoryALFATClust_${seq_length}.txt"
	> "${picture_file}/timeMemoryLinclust_${seq_length}.txt"
	for folder in ${WORK_DIR}/result/testdata${seq_length}*; do
    	if [ -d "${folder}/groundTruth" ] && [ -f "${folder}/groundTruth/_TimeMemory.txt" ]; then
        	folder_name=$(basename "${folder}")
        	fasta_file="${WORK_DIR}/${folder_name}.fa"
        	if [ -f "${fasta_file}" ]; then
            		num_seq=$(grep -c ">" "${fasta_file}")
            	while IFS=$'\t' read -r col1 col2 col3; do
                	echo -e "${num_seq}\t${col2}\t${col3}" >> "${picture_file}/timeMemoryBruteForce_${seq_length}.txt"
            	done < "${folder}/groundTruth/_TimeMemory.txt"
        	fi
    	fi

	if [ -d "${folder}/ALFATClust" ] && [ -f "${folder}/ALFATClust/_TimeMemory.txt" ]; then 
                folder_name=$(basename "${folder}")
                fasta_file="${WORK_DIR}/${folder_name}.fa"
                if [ -f "${fasta_file}" ]; then
                        num_seq=$(grep -c ">" "${fasta_file}")
                while IFS=$'\t' read -r col1 col2 col3; do
			if [ "${col1}" == "${sim_threshold_ALFATClust}" ]; then
                        	echo -e "${num_seq}\t${col2}\t${col3}" >> "${picture_file}/timeMemoryALFATClust_${seq_length}.txt"
                	fi
		done < "${folder}/ALFATClust/_TimeMemory.txt"
                fi
        fi

	if [ -d "${folder}/Linclust" ] && [ -f "${folder}/Linclust/_TimeMemory.txt" ]; then
                folder_name=$(basename "${folder}")
                fasta_file="${WORK_DIR}/${folder_name}.fa"
                if [ -f "${fasta_file}" ]; then
                        num_seq=$(grep -c ">" "${fasta_file}")
                while IFS=$'\t' read -r col1 col2 col3; do
                        if [ "${col1}" == "${sim_threshold_Linclust}" ]; then
                                echo -e "${num_seq}\t${col2}\t${col3}" >> "${picture_file}/timeMemoryLinclust_${seq_length}.txt"
                        fi
                done < "${folder}/Linclust/_TimeMemory.txt"
                fi
        fi
	done
done

python ${PROJECT_DIR}/TimeMemoryDraw.py ${number_alg} ${label} ${time_memory_file} ${picture_file}

