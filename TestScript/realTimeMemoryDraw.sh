#!/bin/bash
WORK_DIR=/data/cabins/yxxiang/real
PROJECT_DIR=/data/cabins/yxxiang/HomoPic/project

number_alg=12
label=Linclust_50,Linclust_60,Linclust_70,Linclust_80,Linclust_90,ALFATClust_50,ALFATClust_60,ALFATClust_70,ALFATClust_80,ALFATClust_90,BruteForce_70,BruteForce_90
picture_file=${WORK_DIR}/result/pic
time_memory_file=""

for l in $(echo $label | tr ',' ' '); do
  time_memory_file="${time_memory_file},${picture_file}/timeMemory${l}.txt"
done
time_memory_file=${time_memory_file#","}

for sim_threshold in 50 60 70 80 90; do
	if [ "${sim_threshold}" == "70" ] || [ "${sim_threshold}" == "90" ]; then
		> "${picture_file}/timeMemoryBruteForce_${sim_threshold}.txt"
	fi
	> "${picture_file}/timeMemoryALFATClust_${sim_threshold}.txt"
	> "${picture_file}/timeMemoryLinclust_${sim_threshold}.txt"
done

for folder in ${WORK_DIR}/result/*; do
	if [ -d "${folder}/groundTruth" ] && [ -f "${folder}/groundTruth/_TimeMemory.txt" ]; then
        	folder_name=$(basename "${folder}")
        	fasta_file="${WORK_DIR}/${folder_name}.fasta"
        	if [ -f "${fasta_file}" ]; then
            		num_seq=$(grep -c ">" "${fasta_file}")
            		while IFS=$'\t' read -r col1 col2 col3; do
				if [ "${col1}" == "0.7" ]; then
                			echo -e "${num_seq}\t${col2}\t${col3}" >> "${picture_file}/timeMemoryBruteForce_70.txt"
            			fi
				if [ "${col1}" == "0.9" ]; then
                                        echo -e "${num_seq}\t${col2}\t${col3}" >> "${picture_file}/timeMemoryBruteForce_90.txt"
                                fi
			done < "${folder}/groundTruth/_TimeMemory.txt"
        	fi
    	fi

	if [ -d "${folder}/ALFATClust" ] && [ -f "${folder}/ALFATClust/ALFATSparseMatrix_0.50.txt" ]; then 
                folder_name=$(basename "${folder}")
                fasta_file="${WORK_DIR}/${folder_name}.fasta"
                if [ -f "${fasta_file}" ]; then
                        num_seq=$(grep -c ">" "${fasta_file}")
                	while IFS=$'\t' read -r col1 col2 col3; do
				if [ "${col1}" == "0.50" ]; then
                        		echo -e "${num_seq}\t${col2}\t${col3}" >> "${picture_file}/timeMemoryALFATClust_50.txt"
                		fi
				if [ "${col1}" == "0.60" ]; then
                                        echo -e "${num_seq}\t${col2}\t${col3}" >> "${picture_file}/timeMemoryALFATClust_60.txt"
                                fi
				if [ "${col1}" == "0.70" ]; then
                                        echo -e "${num_seq}\t${col2}\t${col3}" >> "${picture_file}/timeMemoryALFATClust_70.txt"
                                fi
				if [ "${col1}" == "0.80" ]; then
                                        echo -e "${num_seq}\t${col2}\t${col3}" >> "${picture_file}/timeMemoryALFATClust_80.txt"
                                fi
				if [ "${col1}" == "0.90" ]; then
                                        echo -e "${num_seq}\t${col2}\t${col3}" >> "${picture_file}/timeMemoryALFATClust_90.txt"
                                fi
			done < "${folder}/ALFATClust/_TimeMemory.txt"
                fi
        fi

	if [ -d "${folder}/Linclust" ] && [ -f "${folder}/Linclust/_TimeMemory.txt" ]; then
                folder_name=$(basename "${folder}")
                fasta_file="${WORK_DIR}/${folder_name}.fasta"
                if [ -f "${fasta_file}" ]; then
                        num_seq=$(grep -c ">" "${fasta_file}")
                	while IFS=$'\t' read -r col1 col2 col3; do
                        	if [ "${col1}" == "0.50" ]; then
                                	echo -e "${num_seq}\t${col2}\t${col3}" >> "${picture_file}/timeMemoryLinclust_50.txt"
                        	fi
				if [ "${col1}" == "0.60" ]; then
                                        echo -e "${num_seq}\t${col2}\t${col3}" >> "${picture_file}/timeMemoryLinclust_60.txt"
                                fi
				if [ "${col1}" == "0.70" ]; then
                                        echo -e "${num_seq}\t${col2}\t${col3}" >> "${picture_file}/timeMemoryLinclust_70.txt"
                                fi
				if [ "${col1}" == "0.80" ]; then
                                        echo -e "${num_seq}\t${col2}\t${col3}" >> "${picture_file}/timeMemoryLinclust_80.txt"
                                fi
				if [ "${col1}" == "0.90" ]; then
                                        echo -e "${num_seq}\t${col2}\t${col3}" >> "${picture_file}/timeMemoryLinclust_90.txt"
                                fi
                	done < "${folder}/Linclust/_TimeMemory.txt"
                fi
        fi
	done

python ${PROJECT_DIR}/TimeMemoryDraw.py ${number_alg} ${label} ${time_memory_file} ${picture_file}

