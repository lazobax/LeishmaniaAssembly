#!/bin/bash

set -e

# Define variables
read_path_str="data"
outdir="data/filtered"
DBpath="DB"
threads=8
adapterlength=44
pctmatch=97
x="onePerDownSample" # Replace with your sample identifier
mkdir -p $outdir

# Check if the input FASTQ file exists
if [ -s ${read_path_str}/${x}.fastq ]; then

    # Decompress FASTQ file (commented out for this example)
    # pigz -cd -p ${threads} ${read_path_str}/${x}.fastq.gz > ${read_path_str}/${x}.fastq

    echo "Convert FASTQ to FASTA"
    cat ${read_path_str}/${x}.fastq | sed -n '1~4s/^@/>/p;2~4p' > ${read_path_str}/${x}.fasta

    echo "Run BLAST to identify contamination"
    blastn -db $DBpath/pacbio_vectors_db -query ${read_path_str}/${x}.fasta \
        -num_threads ${threads} -task blastn -reward 1 -penalty -5 -gapopen 3 \
        -gapextend 3 -dust no -soft_masking true -evalue 700 -searchsp 1750000000000 \
        -outfmt 6 > ${outdir}/${x}.contaminant.blastout

    echo "Create blocklist of contaminated reads"
    cat ${outdir}/${x}.contaminant.blastout | grep 'NGB0097' \
        | awk -v OFS='\t' -v var1="${adapterlength}" -v var2="${pctmatch}" \
          '{if (($2 ~ /NGB00972/ && $3 >= var2 && $4 >= var1) || ($2 ~ /NGB00973/ && $3 >= 97 && $4 >= 34)) print $1}' \
        | sort -u > ${outdir}/${x}.blocklist

    echo "Remove contaminated reads from FASTQ"
    cat ${read_path_str}/${x}.fastq | paste - - - - | grep -v -f ${outdir}/${x}.blocklist -F \
        | tr "\t" "\n" | pigz -p ${threads} --fast > ${outdir}/${x}.filt.fastq.gz

    # Calculate and write statistics
    f=$(cat ${outdir}/${x}.blocklist | wc -l) # Number of adapter contaminated reads
    r1=$(cat ${read_path_str}/${x}.fastq | wc -l) 
    r2=$(awk -v r1=$r1 'BEGIN{ans=r1/4; print ans}') # Number of reads
    p1=$(awk -v n1=$r2 -v n2=$f 'BEGIN{ans=n2/n1*100; print ans}') # Proportion of adapter contaminated reads
    r3=$(awk -v r2=$r2 -v f=$f 'BEGIN{ans=r2-f; print ans}') # Number of reads retained
    p2=$(awk -v p1=$p1 'BEGIN{ans=100-p1; print ans}') # Proportion of reads retained

    echo "For the ${x} dataset:" >> ${outdir}/${x}.stats
    echo "Removing reads containing adapters a minimum of ${adapterlength} bp in length and ${pctmatch}% match." >> ${outdir}/${x}.stats
    echo "" >> ${outdir}/${x}.stats
    echo "Number of reads:" $r2 >> ${outdir}/${x}.stats
    echo "Number of adapter contaminated reads:" $f "($p1% of total)" >> ${outdir}/${x}.stats
    echo "Number of reads retained:" $r3 "($p2% of total)" >> ${outdir}/${x}.stats
    echo "" >> ${outdir}/${x}.stats
    echo "Finished on $(date)" >> ${outdir}/${x}.stats

    # Clean up intermediate files
    # rm ${read_path_str}/${x}.fastq
    # rm ${read_path_str}/${x}.fasta
fi
