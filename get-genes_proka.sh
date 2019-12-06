rm -f ../ref/transcriptoma.fa

IFS=$'\n'
for accline in $(cat ../ref/ACCs.txt); do
        acc=`echo ${accline} | cut -f 1`
        seqref=`echo ${accline} | cut -f 2`
        chr_start=`echo ${accline} | cut -f 3`
        chr_stop=`echo ${accline} | cut -f 4`
        strand=`echo ${accline} | cut -f 5`

        echo "Pegando FASTA para ${acc}  [${seqref}:${chr_start}-${chr_stop}(${strand})] ..."

        efetch -db nucleotide -id ${seqref}  -format fasta \
        -chr_start ${chr_start} \
        -chr_stop ${chr_stop} \
        -strand ${strand} | \
        sed "s/^>.*/>${acc}/" \
        >> ../ref/transcriptoma.fa
done
