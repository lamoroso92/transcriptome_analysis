#!/bin/bash


# Criando pastas para armazenamento dos resultados
if [ ! -d ./output/processed/fastqc/pre/ ]; then
        echo "Criando diretorios FASTQC"
        mkdir -p ./output/processed/fastqc/pre;
fi;

if [ ! -d ./output/processed/fastqc/pos/ ]; then
        mkdir -p ./output/processed/fastqc/pos;
fi;

if [ ! -d ./output/processed/atropos/ ]; then
        echo "Criando diretorio Atropos"
        mkdir -p ./output/processed/atropos;
fi;

if [ ! -d ./output/processed/prinseq/ ]; then
        echo "Criando diretorio PRINSEQ"
        mkdir -p ./output/processed/prinseq;
fi;


# Declarando a variável de arquivos de entrada
indir="$1"

for f in `ls ${indir}/*_R1.fastq`; do
        bn=`basename ${f} _R1.fastq`

        # FASTQC - Obtendo as informações de qualidade antes do processamento
        echo "[FASTQC] Obtendo as estatisticas da biblioteca ${bn} antes do processamento..."
        fastqc -t 2 \
                ${indir}/${bn}_R1.fastq \
                -o ./output/processed/fastqc/pre/

        fastqc -t 2 \
                ${indir}/${bn}_R2.fastq \
                -o ./output/processed/fastqc/pre/

        # ATROPOS - Identicando e Removendo sequencias de adaptadores e leituras de baixa qualidade/curtas demais. ######################
        # Parametros relevantes: --aligner: tipo de algoritmo utilizado para detecção e remoção de adaptador; -e: taxa de erro permitida;
        # -n: nº de adaptadores a serem removidos; -m: tamanho minimo para remoção da read; --op-order: ordem das operações; ############
        # -O: sobreposição minima entre adaptador e leitura; -q: qualidade minima para corte na extremidade 5'; #########################
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
             -o ./output/processed/atropos/${bn}_R1.atropos_insert.fastq \
             -p ./output/processed/atropos/${bn}_R2.atropos_insert.fastq \
             -pe1 ${indir}/${bn}_R1.fastq \
             -pe2 ${indir}/${bn}_R2.fastq \
             --untrimmed-output ./output/processed/atropos/${bn}_R1.atropos_untrimmed.fastq \
             --untrimmed-paired-output ./output/processed/atropos/${bn}_R2.atropos_untrimmed.fastq \
             > ./output/processed/atropos/${bn}.atropos.log.out.txt \
             2> ./output/processed/atropos/${bn}.atropos.log.err.txt


        echo "[Atropos] Realizando a remocao de adaptadores com Atopros no modo ADAPTER"
        atropos trim    --aligner adapter \
             -e 0.1 \
             -n 2 \
             -m 1 \
             --match-read-wildcards \
             -O 3 \
             -q 20 \
             --pair-filter both \
             -pe1 ./output/processed/atropos/${bn}_R1.atropos_untrimmed.fastq \
             -pe2 ./output/processed/atropos/${bn}_R2.atropos_untrimmed.fastq \
             -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
             -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT  \
             -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
             -G CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \
             -T 2 \
             -o ./output/processed/atropos/${bn}_R1.atropos_adapter.fastq  \
             -p ./output/processed/atropos/${bn}_R2.atropos_adapter.fastq \
             > ./output/processed/atropos/${bn}.atropos_adapter.log.out.txt \
             2>  ./output/processed/atropos/${bn}.atropos_adapter.log.err.txt

        pyenv deactivate # Encerrando ambiente ATROPOS

        # Concatenando as reads processadas nos dois modos do ATROPOS em um único arquivo
        cat     ./output/processed/atropos/${bn}_R1.atropos_insert.fastq \
        ./output/processed/atropos/${bn}_R1.atropos_adapter.fastq \
        > ./output/processed/atropos/${bn}_R1.atropos_final.fastq

        cat     ./output/processed/atropos/${bn}_R2.atropos_insert.fastq \
        ./output/processed/atropos/${bn}_R2.atropos_adapter.fastq \
        > ./output/processed/atropos/${bn}_R2.atropos_final.fastq

        # Removendo arquivos anteriores à fusão
        rm -f ./output/processed/atropos/${bn}_R1.atropos_insert.fastq \
        ./output/processed/atropos/${bn}_R1.atropos_adapter.fastq \
        ./output/processed/atropos/${bn}_R2.atropos_insert.fastq \
        ./output/processed/atropos/${bn}_R2.atropos_adapter.fastq

        # PRINSEQ - Processamento de regiões e reads de baixa qualidade e/ou tamanho curto ############################################################
        # Parametros relevantes: -out_format: Define o tipo de saida, 3 se refere a FASTQ; ############################################################
        # -trim_qual_***: Parametros referentes ao processamento por qualidade; -trim_tail: Define tamanho de corte em extremidades (right e left); ###
        # -out_good: Local/nome de saída das reads processadas; -out_bad: Local de saída das reads de baixa qualidade que foram removidas no processo #
        # (Null para não obter esse arquivo); -lc_method: Metódo de detecção e remoção de regiões de baixa complexidade; ##############################
        echo "[PRINSEQ] Realizando o processamento da sequencia ${bn}.fastq"
        prinseq-lite.pl -fastq ./output/processed/atropos/${bn}_R1.atropos_final.fastq \
                -fastq2 ./output/processed/atropos/${bn}_R2.atropos_final.fastq \
                -out_format 3 \
                -trim_qual_window 3 \
                -trim_qual_step 1 \
                -trim_qual_right 30 \
                -trim_qual_type mean \
                -trim_qual_rule lt \
                -out_good ./output/processed/prinseq/${bn}.prinseq \
                -out_bad null \
                -lc_method dust \
                -lc_threshold 30 \
                -min_len 20 \
                -trim_tail_right 5 \
                -trim_tail_left 5\
                -ns_max_p 80 \
                -noniupac \
                > ./output/processed/prinseq/${bn}.prinseq.out.log \
                2> ./output/processed/prinseq/${bn}.prinseq.err.log

        # FASTQC - Obtendo as informações de qualidade após o processamento
        echo "[FASTQC] Obtendo as estatisticas da biblioteca ${bn} apos processamento..."
        fastqc -t 2 \
                ./output/processed/prinseq/${bn}.prinseq_1.fastq \
                -o ./output/processed/fastqc/pos/

        fastqc -t 2 \
                ./output/processed/prinseq/${bn}.prinseq_2.fastq \
                -o ./output/processed/fastqc/pos/

        # SE EXISTIR <SAMPLE_NAME>.prinseq_1_singletons.fastq
        if [ -e "./output/processed/prinseq/${bn}.prinseq_1_singletons.fastq" ]; then
        fastqc -t 2 \
                ./output/processed/prinseq/${bn}.prinseq_1_singletons.fastq \
                -o ./output/processed/fastqc/pos/
        fi

        # SE EXISTIR <SAMPLE_NAME>.prinseq_2_singletons.fastq
        if [ -e "./output/processed/prinseq/${bn}.prinseq_2_singletons.fastq" ]; then
        fastqc -t 2 \
                ./output/processed/prinseq/${bn}.prinseq_2_singletons.fastq \
                -o ./output/processed/fastqc/pos/
        fi
done
