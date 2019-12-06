#!/bin/bash

for biogroup in A B; do
	for rep in 1 2; do
		echo "Gerando reads para amostra ${biogroup} réplica ${rep} ..."
		
		# Gerando as sequencias:
		# -r: local do arquivo contendo os 'fasta' dos genes de interesse
		# -a: local dos arquivos de abundancia
		# -t: quantidade de sequencias geradas
		# -i: tamanho das leituras em pb
		# -s: desvio padrão das sequencias considerando a abundancia fornecida

		generate_fragments.py -r ../ref/transcriptoma.fa \
		   -a ../ref/abundance_${biogroup}.txt \
		   -o ../tmp.frags_${biogroup}_${rep} \
		   -t 25000 \
		   -i 300 \
		   -s 30
		
		# Renomeando as sequencias e adicionando identificador numerico para cada uma
		cat ../tmp.frags_${biogroup}_${rep}.1.fasta | renameSeqs.pl \
		   -if FASTA \
		   -of FASTA \
		   -p SAMPLE${biogroup}${rep} \
		   -w 1000 | \
		   sed 's/^>\(\S\+\).*/>\1/' \
		   > ../frags_${biogroup}${rep}.fa

		# Transformando as sequencias geradas em "paired-end" 2x151 com adaptadores para tornar mais realistico a simulacao
		cat ../frags_${biogroup}${rep}.fa | simNGS -a \
		AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
		-p paired \
		/usr/local/bioinfo/simNGS/data/s_4_0099.runfile \
		-n 151 > ../SAMPLE${biogroup}${rep}.fastq 2> ../SAMPLE${biogroup}${rep}.err.txt
		
		# Criando local de saida
		mkdir -p ../raw
		
		# Separando os pares em R1 e R2
		deinterleave_pairs ../SAMPLE${biogroup}${rep}.fastq \
		   -o ../raw/SAMPLE${biogroup}${rep}_R1.fastq \
		      ../raw/SAMPLE${biogroup}${rep}_R2.fastq

		# Removendo arquivos intermediarios
		rm -f ../tmp.frags_${biogroup}_${rep}.1.fasta ../frags_${biogroup}${rep}.fa ../SAMPLE${biogroup}${rep}.fastq ../SAMPLE${biogroup}${rep}.err.txt

		echo "Número de reads ${biogroup}${rep} R1:" $(echo "$(cat ../raw/SAMPLE${biogroup}${rep}_R1.fastq | wc -l)/4" | bc)
		echo "Número de reads ${biogroup}${rep} R2:" $(echo "$(cat ../raw/SAMPLE${biogroup}${rep}_R2.fastq | wc -l)/4" | bc)

	done
done
