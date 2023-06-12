#!/bin/bash
WORK_DIR=/data/cabins/yxxiang/real                  # The folder of fasta file
PROJECT_DIR=/data/cabins/yxxiang/HomoPic/project        # The direction of your project saved
sim_list="0.50 0.60 0.70 0.80 0.90"                     # The list of similarity thresholds for Linclust and ALFATClust

# Looking up .fasta file in the WORK_DIR
mkdir -p ${WORK_DIR}/result
find ${WORK_DIR} -maxdepth 1 \( -name "*.fa" -o -name "*.fasta" \) | while read file; do
	## Creating folders with the same name for each data set
	filename=$(basename ${file})
	echo ${filename}
	folder_dir=${WORK_DIR}/result/${filename%.*}
    	mkdir -p ${folder_dir}
	
	## Making lookup table
	ID=0
	while read -r line
	do
		if [ "${line:0:1}" = ">" ];then
			echo -e $line"\t"$ID >> ${folder_dir}/index.lookup
			let ID++
		fi
	done < ${file}
	
	## Compiling clusterMon
	make -f ${PROJECT_DIR}/clusterMon/makefile clean >> ${PROJECT_DIR}/clusterMon/Compilation.log
	make -f ${PROJECT_DIR}/clusterMon/makefile -j >> ${PROJECT_DIR}/clusterMon/Compilation.log

	## Calling Linclust
	mkdir -p ${folder_dir}/Linclust
	lin_dir=${folder_dir}/Linclust
	
	source ${PROJECT_DIR}/conda.sh
	conda activate mmseq2
	for T in ${sim_list}
	do
		echo "=== Linclust ${filename} with a similarity threshold of $T is running..."
		echo -n -e $T"\t" >> ${lin_dir}/_TimeMemory.txt
		/usr/bin/time --format="%e\t%M" -a -o ${lin_dir}/_TimeMemory.txt mmseqs easy-cluster ${file} ${lin_dir}/Linclu_${T} ${lin_dir}/tmp --min-seq-id $T --threads 1 > ${lin_dir}/log
		### Generating sparse matrices for Linclust
		${PROJECT_DIR}/clusterMon/clusterMon -m Linclust -i ${lin_dir}/Linclu_${T}_all_seqs.fasta -o ${lin_dir}/LinSparseMatrix_$T.txt -l ${folder_dir}/index.lookup
	done
	conda deactivate
	echo ""


	## Calling ALFATClust
	mkdir -p ${folder_dir}/ALFATClust
	alfat_dir=${folder_dir}/ALFATClust
	
	update_file="${alfat_dir}/${filename%.*}_updata.fasta"
	awk '{
  		if (substr($0,1,1)==">") {
    			gsub(/[^[:alnum:]>]/,"")
    			print
  		} else {
    			print
  		}
	}' ${file} > ${update_file}
	### Making lookup table for ALFAT
        ID=0
        while read -r line
        do
                if [ "${line:0:1}" = ">" ];then
                        echo -e $line"\t"$ID >> ${alfat_dir}/index_ALFAT.lookup
                        let ID++
                fi
        done < ${update_file}

	source ${PROJECT_DIR}/conda.sh
	conda activate ALFATClust
	for T in ${sim_list}
        do
                echo "=== ALFATClust ${filename} with a similarity threshold of $T is running..."
                echo -n -e $T"\t" >> ${alfat_dir}/_TimeMemory.txt
                /usr/bin/time --format="%e\t%M" -a -o ${alfat_dir}/_TimeMemory.txt ./alfatclust.py -i ${update_file} -o ${alfat_dir}/ALFATclu_$T.fasta -p -l $T -t 1 > ${alfat_dir}/log
		### Generating sparse matrices for ALFATClust
                ${PROJECT_DIR}/clusterMon/clusterMon -m ALFATClust -i ${alfat_dir}/ALFATclu_$T.fasta -o ${alfat_dir}/ALFATSparseMatrix_$T.txt -l ${alfat_dir}/index_ALFAT.lookup
        done
	conda deactivate
	echo ""

	## Calling LSH
	### T.B.A.

done
