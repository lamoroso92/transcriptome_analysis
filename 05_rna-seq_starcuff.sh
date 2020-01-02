#!/bin/bash

## Definindo e checando as variáveis

## Número de threads
num_threads=10

## Diretório de saída dos resultados 
outdir="$1"

if [ ! ${outdir} ]; then
	echo "Missing output directory."
	exit
fi

if [ ! -d ${outdir} ]; then
	echo "Wrong output directory (${outdir})."
	exit
fi

## Local e nome do arquivo GTF
refgtf="$2"

if [ ! ${refgtf} ]; then
	echo "Missing GTF file."
	exit
fi

if [ ! -e "${refgtf}" ]; then
	echo "Not found GTF file (${refgtf})."
	exit
fi

## Local e nome do arquivo do Genoma
refseq="$3"

if [ ! ${refseq} ]; then
	echo "Missing GENOME fasta file."
	exit
fi

if [ ! -e "${refseq}" ]; then
	echo "Not found GENOME fasta file (${refseq})."
	exit
fi

#########################################################################################################################################

## Criando pastas de outputs
mkdir -p ${outdir}/01_star_index
mkdir -p ${outdir}/02_star_out_pe
mkdir -p ${outdir}/03_star_out_se
mkdir -p ${outdir}/04_star_out_final
mkdir -p ${outdir}/05_cufflinks
mkdir -p ${outdir}/06_cuffmerge
mkdir -p ${outdir}/07_cuffcompare
mkdir -p ${outdir}/08_cuffquant
mkdir -p ${outdir}/09_cuffnorm
mkdir -p ${outdir}/10_cuffdiff

#########################################################################################################################################

## Obtendo nomes dos arquivos de Reads (Nome das anmostas)
## Criando arquivos singletons para aqueles que não tiveram, assim mantendo a compatibilidade com os programas
for r1 in `ls ${outdir}/00_processed/03_prinseq/*.prinseq_1.fastq`; do
	r1_singletons=`echo ${r1} | sed 's/prinseq_1.fastq/prinseq_1_singletons.fastq/'`
	if [ ! -e "${r1_singletons}" ]; then
		touch ${r1_singletons}
	fi
	
	r2=`echo ${r1} | sed 's/prinseq_1.fastq/prinseq_2.fastq/'`
	
	if [ ! -e "${r2}" ]; then
		echo "Read 2 (${r2}) paired with Read 1 ($r1) not found."
		exit
	fi
	
	r2_singletons=`echo ${r2} | sed 's/prinseq_2.fastq/prinseq_2_singletons.fastq/'`
	if [ ! -e "${r2_singletons}" ]; then
		touch ${r2_singletons}
	fi

	name=`basename ${r1} | sed 's/.prinseq_1.fastq//'`
	
	## Indexando Genoma
	if [ ! -e "${outdir}/star_index/SAindex" ]; then
		echo "Indexing genome (${refseq}) ..."
		STAR 	--runThreadN        ${num_threads} \
   			--runMode           genomeGenerate \
     			--genomeFastaFiles  ${refseq} \
     			--genomeDir         ${outdir}/01_star_index \
    	 		--sjdbGTFfile       ${refgtf} \
			--genomeSAindexNbases 10 \
     			--sjdbOverhang      149 \
			--alignSJDBoverhangMin 999 \
		 > ${outdir}/01_star_index/STAR.index.log.out.txt \
		2> ${outdir}/01_star_index/STAR.index.log.err.txt

	fi

#########################################################################################################################################

## Alinhando reads paired-end no genoma
echo "STAR alignment PE with sample ${name}: ${r1} & ${r2} ..."
	
mkdir -p ${outdir}/02_star_out_pe/${name}
	
STAR 	--runThreadN  	    ${num_threads} \
       	--genomeDir	    ${outdir}/01_star_index \
       	--readFilesIn       ${r1} ${r2} \
       	--sjdbGTFfile	    ${refgtf} \
        --alignIntronMax    1 \
        --alignMatesGapMax  0 \
       	--outFileNamePrefix ${outdir}/02_star_out_pe/${name}/ \
	--outSAMtype        BAM Unsorted \
	--outSAMstrandField intronMotif \
	> ${outdir}/02_star_out_pe/${name}/STAR.alignment_pe.log.out.txt \
	2> ${outdir}/02_star_out_pe/${name}/STAR.alignment_pe.log.err.txt

#########################################################################################################################################

## Alinhando singletons no genoma
echo "STAR alignment SE with sample ${name}: ${r1_singletons} & ${r2_singletons} ..."
	
mkdir -p ${outdir}/03_star_out_se/${name}
	
STAR 	--runThreadN  	    ${num_threads} \
       	--genomeDir	    ${outdir}/01_star_index \
       	--readFilesIn       ${r1_singletons},${r2_singletons} \
       	--sjdbGTFfile	    ${refgtf} \
	--alignIntronMax    1 \
	--alignMatesGapMax  0 \
	--outSAMtype        BAM Unsorted \
	--outSAMstrandField intronMotif \
        --outFileNamePrefix ./$outdir/03_star_out_se/${name}/ \
	> ./${outdir}/03_star_out_se/${name}/STAR.alignment_se.log.out.txt \
	2> ./${outdir}/03_star_out_se/${name}/STAR.alignment_se.log.err.txt

#########################################################################################################################################

## Fundindo arquivos .bam das reads paired-end e singletons de cada amostra
echo "Merging STAR alignment PE & SE ..."
	
mkdir -p ${outdir}/04_star_out_final/${name}

## Combinar resultados do alinhamento com reads paired-end e alinhamento com reads single-end (singletons)       
samtools merge -@ ${num_threads} -f -n  ${outdir}/04_star_out_final/${name}/Aligned.out.bam \
	${outdir}/02_star_out_pe/${name}/Aligned.out.bam \
	${outdir}/03_star_out_se/${name}/Aligned.out.bam

#########################################################################################################################################

## Ordenando o resultado do alinhamento por coordenadas genomicas
## Exigencia para executar o cufflinks
echo "Sorting STAR alignment final ..."

samtools sort -@ ${num_threads} -o ${outdir}/04_star_out_final/${name}/Aligned.out.sorted.bam \
	${outdir}/04_star_out_final/${name}/Aligned.out.bam

#########################################################################################################################################

## Obtendo as estatísticas dos alinhamentos
echo "Collecting alignment statistics ..."
	
SAM_nameSorted_to_uniq_count_stats.pl ${outdir}/04_star_out_final/${name}/Aligned.out.bam > ${outdir}/04_star_out_final/${name}/Aligned.stats.txt

#########################################################################################################################################	

## Cufflinks: Montagem dos transcritos
echo "Running Cufflinks ..."

mkdir -p ${outdir}/05_cufflinks/${name}

cufflinks --output-dir ${outdir}/05_cufflinks/${name} \
	--num-threads ${num_threads} \
	--GTF-guide ${refgtf} \
	--frag-bias-correct ${refseq} \
	--multi-read-correct \
	--library-type fr-unstranded \
	--frag-len-mean 300 \
	--frag-len-std-dev 50 \
	--total-hits-norm \
	--max-frag-multihits 20 \
	--min-isoform-fraction 0.20 \
	--max-intron-length 1 \
	--min-intron-length 0 \
	--overhang-tolerance 4 \
	--max-bundle-frags 999999 \
	--max-multiread-fraction 0.45 \
	--overlap-radius 10 \
	--3-overhang-tolerance 300 \
	${outdir}/04_star_out_final/${name}/Aligned.out.sorted.bam

done
#########################################################################################################################################

## Cuffmerge: Fusão das montagens dos transcritos
echo "Running cuffmerge ..."

find ${outdir}/05_cufflinks/ -name 'transcripts.gtf' > ${outdir}/06_cuffmerge/assembly_GTF_list.txt

cuffmerge -o ${outdir}/06_cuffmerge/ \
        --ref-gtf ${refgtf} \
        --ref-sequence ${refseq} \
        --min-isoform-fraction 0.20 \
        --num-threads ${num_threads} \
         ${outdir}/06_cuffmerge/assembly_GTF_list.txt \
         > ${outdir}/06_cuffmerge/cuffmerge.log.out.txt \
        2> ${outdir}/06_cuffmerge/cuffmerge.log.err.txt

## Obtendo a contagem de códigos de classes do Cufflinks
perl -F"\t" -lane 'my ($transcript_id)=$F[8]=~/transcript_id "([^"]+)"/; my ($class_code)=$F[8]=~/class_code "([^"]+)"/; print join("\t", $transcript_id, $class_code);' \
   ./output/06_cuffmerge/merged.gtf \
   | nsort -u \
   | cut -f 2 \
   | nsort \
   | uniq -c \
   | sed 's/^ *//' \
   | awk 'BEGIN{OFS="\t";}{ print $2,$1;} ' \
   | sort -t$'\t' -k 2nr,2nr > ${outdir}/06_cuffmerge/class_code.count.txt

#########################################################################################################################################

## Cuffcompare: Compara e atribui estatísticas dos transcritos. Também pode identificar transcritos variantes
echo "Running cuffcompare ..."
cuffcompare     -r ${refgtf} \
                -s ${refseq} \
                -o ${outdir}/07_cuffcompare/cuffmerge \
                ${outdir}/06_cuffmerge/merged.gtf \
                 > ${outdir}/06_cuffmerge/cuffcompare.log.out.txt \
                2> ${outdir}/06_cuffmerge/cuffcompare.log.err.txt

biogroup_label=()
for bamfile in `ls ${outdir}/04_star_out_final/*/Aligned.out.sorted.bam`; do
        name=`basename $(dirname ${bamfile})`
        echo "Running cuffquant using sample ${name} with ${outdir}/06_cuffmerge/merged.gtf as reference ..."
        mkdir -p ${outdir}/08_cuffquant/${name}
	
	## Cuffquant: Gera os arquivos de quantificação das expressões que serão utilizados no cuffnorm/cuffdiff
        cuffquant       --output-dir ${outdir}/08_cuffquant/${name} \
                        --frag-bias-correct ${refseq} \
                        --multi-read-correct \
                        --num-threads ${num_threads} \
                        --library-type fr-unstranded \
                        --frag-len-mean 300 \
                        --frag-len-std-dev 50 \
                        --max-bundle-frags 9999999 \
                        --max-frag-multihits 20 \
                        ${outdir}/06_cuffmerge/merged.gtf \
                        ${bamfile} \
                 > ${outdir}/08_cuffquant/${name}/cuffquant.log.out.txt \
                2> ${outdir}/08_cuffquant/${name}/cuffquant.log.err.txt

        groupname=`echo ${name} | sed 's/[0-9]\+$//'`
        biogroup_label=($(printf "%s\n" ${biogroup_label[@]} ${groupname} | sort -u ))

done
biogroup_files=()

#########################################################################################################################################
echo "Running Differential Expression Analysis ..."
for label in ${biogroup_label[@]}; do
        echo -e "\tCollecting .cxb files for ${label} ..."
        group=()
        for cxbfile in `ls ${outdir}/08_cuffquant/${label}*/abundances.cxb`; do
                echo -e "\t\tFound ${cxbfile}"
                group=(${group[@]} "${cxbfile}")
        done
        biogroup_files=(${biogroup_files[@]} $(IFS=, ; echo "${group[*]}") )
done

## Cuffnorm: Normaliza as contagens absolutas
echo -e "\tRunning cuffnorm & cuffdiff ..."
echo -e "\t\tLabels.: " $(IFS=, ; echo "${biogroup_label[*]}")
echo -e "\t\tFiles..: " ${biogroup_files[*]}

echo -e "\t\t\tGenerating abundance matrices (cuffnorm) ..."

cuffnorm        --output-dir ${outdir}/09_cuffnorm \
                --labels $(IFS=, ; echo "${biogroup_label[*]}") \
                --num-threads ${num_threads} \
                --library-type fr-unstranded \
                --library-norm-method geometric \
                --output-format simple-table \
                ${outdir}/06_cuffmerge/merged.gtf \
                ${biogroup_files[*]}

## de-normalize (cuffnorm): obtenção das contagens brutas dos transcritos
echo -e "\t\t\tGenerating raw abundance matrix (De-normalize cuffnorm) ..."
de-normalize-cuffnorm.R \
   --in=${outdir}/09_cuffnorm/genes.count_table \
   --st=${outdir}/09_cuffnorm/samples.table \
   --out=${outdir}/09_cuffnorm/genes.raw.count_table

## Cuffdiff: Identifica os transcritos diferencialmente expressos
echo -e "\t\t\tAnalysing differential expression (cuffdiff) ..."

cuffdiff        --output-dir ${outdir}/10_cuffdiff \
                --labels $(IFS=, ; echo "${biogroup_label[*]}") \
                --frag-bias-correct ${refseq} \
                --multi-read-correct \
                --num-threads ${num_threads} \
                --library-type fr-unstranded \
                --frag-len-mean 300 \
                --frag-len-std-dev 50 \
                --max-bundle-frags 9999999 \
                --max-frag-multihits 20 \
                --total-hits-norm \
                --min-reps-for-js-test 2 \
                --library-norm-method geometric \
                --dispersion-method per-condition \
                --min-alignment-count 10 \
                ${outdir}/06_cuffmerge/merged.gtf \
                ${biogroup_files[*]}

mv ./Log.out ./output/04_star_out_final/Log.out

echo -e ""
echo -e ""
echo -e "\t\t\tAnalysis finished!"
echo -e ""
echo -e ""
