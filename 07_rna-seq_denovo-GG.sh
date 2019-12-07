#!/bin/bash

# INPUT: Local dos arquivos BAM final
input=$1

if [ ! ${input} ]
then   
        echo "Missing input directory"
        exit
else   
        if [ ! -d ${input} ]
        then   
                echo "Wrong input directory ${input}"
                exit
        fi
fi


# OUTPUT: diretório para armazenar o resultado do processo de montagem
output=$2

if [ ! ${output} ]
then   
        echo "Missing output directory"
        exit
else   
        if [ ! -d ${output} ]
        then   
                echo "Wrong output directory ${output}"
                exit
        fi
fi

num_threads="10"
mem_gb="10G"

##################################################################################################################################################

# Arquivos e diretórios de saída (output) 

basedir_out="${output}/12_denovo-trinity_GG"

aligned_out="${basedir_out}/trinity_GG_input"

mkdir -p ${aligned_out}

##################################################################################################################################################

if [ ! -e "${aligned_out}/All.sorted.bam" ]; then
	echo -e "Collecting alignments ..."
	
	bamfiles=()
	
	bamfiles=( $( find ${input} -name 'Aligned.out.sorted.bam' ) )
	
	samtools merge -f ${aligned_out}/All.sorted.bam ${bamfiles[*]}
	
	samtools sort --threads ${num_threads} ${aligned_out}/All.sorted.bam > ${aligned_out}/All.csorted.bam
	
	rm -f ${aligned_out}/All.sorted.bam
fi

##################################################################################################################################################
	
echo -e "Assembling step (Trinity) ..."
	
Trinity --KMER_SIZE 27 \
	--output ${basedir_out} \
	--seqType fq \
	--max_memory ${mem_gb} \
	--CPU ${num_threads} \
	--min_per_id_same_path 95 \
	--max_diffs_same_path  5 \
	--path_reinforcement_distance 5 \
	--group_pairs_distance 500 \
	--min_glue 5 \
	--min_contig_length 600 \
	--min_kmer_cov 3 \
	--genome_guided_bam ${aligned_out}/All.csorted.bam \
	--genome_guided_max_intron 1 \
	 > ${basedir_out}/Trinity.log.out.txt \
	2> ${basedir_out}/Trinity.log.err.txt
