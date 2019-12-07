#!/bin/bash


# Criando pastas para armazenamento dos resultados
if [ ! -d ./output/00_processed/01_fastqc_pre/ ]; then
	echo "Criando diretorios FASTQC" 
	mkdir -p ./output/00_processed/01_fastqc_pre;
fi;

if [ ! -d ./output/00_processed/04_fastqc_pos/ ]; then
    	mkdir -p ./output/00_processed/04_fastqc_pos;
fi;

if [ ! -d ./output/00_processed/02_atropos/ ]; then
	echo "Criando diretorio Atropos"
	mkdir -p ./output/00_processed/02_atropos;
fi;

if [ ! -d ./output/00_processed/03_prinseq/ ]; then
        echo "Criando diretorio PRINSEQ"
        mkdir -p ./output/00_processed/03_prinseq;
fi;


# Declarando a variÃ¡vel de arquivos de entrada
indir="$1"

for f in `ls ${indir}/*_R1.fastq`; do
        bn=`basename ${f} _R1.fastq`
	
	# FASTQC - Obtendo as informaÃ§Ãµes de qualidade antes do processamento	
	echo "[FASTQC] Obtendo as estatisticas da biblioteca ${bn} antes do processamento..."
	fastqc -t 2 \
		${indir}/${bn}_R1.fastq \
   		-o ./output/00_processed/01_fastqc_pre/

	fastqc -t 2 \
   		${indir}/${bn}_R2.fastq \
   		-o ./output/00_processed/01_fastqc_pre/
	
	# ATROPOS - Identicando e Removendo sequencias de adaptadores e leituras de baixa qualidade/curtas demais. ######################
	# Parametros relevantes: --aligner: tipo de algoritmo utilizado para detecÃ§Ã£o e remoÃ§Ã£o de adaptador; -e: taxa de erro permitida;
	# -n: nÂº de adaptadores a serem removidos; -m: tamanho minimo para remoÃ§Ã£o da read; --op-order: ordem das operaÃ§Ãµes; ############
	# -O: sobreposiÃ§Ã£o minima entre adaptador e leitura; -q: qualidade minima para corte na extremidade 5'; #########################
	# -A/a: sequencia de adaptador a ser removida; -o/p: local/nome dos outputs #####################################################

	eval "$(pyenv init -)"
	pyenv activate atropos # Iniciando ambiente para rodar ATROPOS

	
	echo "[Atropos] Realizando a remocao de adaptadores com Atopros no modo INSERT"
	atropos trim --aligner insert \
             -e 0.1 \
             -n 2 \
             -m 1 \
             --op-order GAWCQ \
             --match-read-wildcards \
             -O 20 \
             -q 25 \
             -T 2 \
             --correct-mismatches conservative \
             --pair-filter any \
             -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
             -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT  \
             -o ./output/00_processed/02_atropos/${bn}_R1.atropos_insert.fastq \
             -p ./output/00_processed/02_atropos/${bn}_R2.atropos_insert.fastq \
             -pe1 ${indir}/${bn}_R1.fastq \
             -pe2 ${indir}/${bn}_R2.fastq \
             --untrimmed-output ./output/00_processed/02_atropos/${bn}_R1.atropos_untrimmed.fastq \
             --untrimmed-paired-output ./output/00_processed/02_atropos/${bn}_R2.atropos_untrimmed.fastq \
             > ./output/00_processed/02_atropos/${bn}.atropos.log.out.txt \
             2> ./output/00_processed/02_atropos/${bn}.atropos.log.err.txt
	

        echo "[Atropos] Realizando a remocao de adaptadores com Atopros no modo ADAPTER"
	atropos trim    --aligner adapter \
             -e 0.1 \
             -n 2 \
             -m 1 \
             --match-read-wildcards \
             -O 3 \
             -q 20 \
             --pair-filter both \
             -pe1 ./output/00_processed/02_atropos/${bn}_R1.atropos_untrimmed.fastq \
             -pe2 ./output/00_processed/02_atropos/${bn}_R2.atropos_untrimmed.fastq \
             -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
             -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT  \
             -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
             -G CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \
             -T 2 \
             -o ./output/00_processed/02_atropos/${bn}_R1.atropos_adapter.fastq  \
             -p ./output/00_processed/02_atropos/${bn}_R2.atropos_adapter.fastq \
             > ./output/00_processed/02_atropos/${bn}.atropos_adapter.log.out.txt \
             2>  ./output/00_processed/02_atropos/${bn}.atropos_adapter.log.err.txt

	pyenv deactivate # Encerrando ambiente ATROPOS
	
	# Concatenando as reads processadas nos dois modos do ATROPOS em um Ãºnico arquivo
	cat     ./output/00_processed/02_atropos/${bn}_R1.atropos_insert.fastq \
        ./output/00_processed/02_atropos/${bn}_R1.atropos_adapter.fastq \
	> ./output/00_processed/02_atropos/${bn}_R1.atropos_final.fastq

	cat     ./output/00_processed/02_atropos/${bn}_R2.atropos_insert.fastq \
        ./output/00_processed/02_atropos/${bn}_R2.atropos_adapter.fastq \
   	> ./output/00_processed/02_atropos/${bn}_R2.atropos_final.fastq
	
	# Removendo arquivos anteriores Ã  fusÃ£o
	rm -f ./output/00_processed/02_atropos/${bn}_R1.atropos_insert.fastq \
      	./output/00_processed/02_atropos/${bn}_R1.atropos_adapter.fastq \
      	./output/00_processed/02_atropos/${bn}_R2.atropos_insert.fastq \
      	./output/00_processed/02_atropos/${bn}_R2.atropos_adapter.fastq
	
	# PRINSEQ - Processamento de regiÃµes e reads de baixa qualidade e/ou tamanho curto ############################################################
	# Parametros relevantes: -out_format: Define o tipo de saida, 3 se refere a FASTQ; ############################################################
	# -trim_qual_***: Parametros referentes ao processamento por qualidade; -trim_tail: Define tamanho de corte em extremidades (right e left); ###
	# -out_good: Local/nome de saÃ­da das reads processadas; -out_bad: Local de saÃ­da das reads de baixa qualidade que foram removidas no processo #
	# (Null para nÃ£o obter esse arquivo); -lc_method: MetÃ³do de detecÃ§Ã£o e remoÃ§Ã£o de regiÃµes de baixa complexidade; ############################## 
	echo "[PRINSEQ] Realizando o processamento da sequencia ${bn}.fastq"
        prinseq-lite.pl -fastq ./output/00_processed/02_atropos/${bn}_R1.atropos_final.fastq \
                -fastq2 ./output/00_processed/02_atropos/${bn}_R2.atropos_final.fastq \
                -out_format 3 \
                -trim_qual_window 3 \
                -trim_qual_step 1 \
                -trim_qual_right 30 \
                -trim_qual_type mean \
                -trim_qual_rule lt \
                -out_good ./output/00_processed/03_prinseq/${bn}.prinseq \
                -out_bad null \
                -lc_method dust \
                -lc_threshold 30 \
                -min_len 20 \
                -trim_tail_right 5 \
                -trim_tail_left 5\
                -ns_max_p 80 \
                -noniupac \
                > ./output/00_processed/03_prinseq/${bn}.prinseq.out.log \
                2> ./output/00_processed/03_prinseq/${bn}.prinseq.err.log

	# FASTQC - Obtendo as informaÃ§Ãµes de qualidade apÃ³s o processamento 
	echo "[FASTQC] Obtendo as estatisticas da biblioteca ${bn} apos processamento..."
	fastqc -t 2 \
		./output/00_processed/03_prinseq/${bn}.prinseq_1.fastq \
		-o ./output/00_processed/04_fastqc_pos/
	
	fastqc -t 2 \
		./output/00_processed/03_prinseq/${bn}.prinseq_2.fastq \
		-o ./output/00_processed/04_fastqc_pos/

	# SE EXISTIR <SAMPLE_NAME>.prinseq_1_singletons.fastq
	if [ -e "./output/00_processed/03_prinseq/${bn}.prinseq_1_singletons.fastq" ]; then
	fastqc -t 2 \
		./output/00_processed/03_prinseq/${bn}.prinseq_1_singletons.fastq \
		-o ./output/00_processed/04_fastqc_pos/
	fi

	# SE EXISTIR <SAMPLE_NAME>.prinseq_2_singletons.fastq
	if [ -e "./output/00_processed/03_prinseq/${bn}.prinseq_2_singletons.fastq" ]; then
	fastqc -t 2 \
		./output/00_processed/03_prinseq/${bn}.prinseq_2_singletons.fastq \
		-o ./output/00_processed/04_fastqc_pos/	
	fi
done
