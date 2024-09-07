params.kmerLen = 31
params.readLen = 17000


workflow {

    originalReads = channel.fromPath("Lei012024_HiFi.fastq.gz")
    originalReads.view()

    downSample(originalReads)
    database = channel.fromPath("DB/*").collect()
    adapterFilt(downSample.out, database)

    fastqcReads(adapterFilt.out[0])

    countKmers(adapterFilt.out[0])

    genomeScopeScript = channel.fromPath("genomescope.R")
    genomeScope(countKmers.out, genomeScopeScript)

}

process downSample{

    conda 'bioconda::seqkit=2.8.2'
    cpus 8

    input:
        path originalReads
    output:
        path 'downSampledReads.fastq'
    script:
    """
    seqkit -j ${task.cpus} sample -p 0.02 -s 100 $originalReads > downSampledReads.fastq
    """

}

process adapterFilt{

    conda 'blast=2.10.1 pigz'
    publishDir "results/qc/adapterFilt"

    input:
        path reads //fastq file
        path "DB/*"
    output:
        path "reads.filt.fastq"
        path "stats.txt"
    script:
    def adapterlength=44
    def pctmatch=97
    """
    echo "Convert FASTQ to FASTA"
    cat ${reads} | sed -n '1~4s/^@/>/p;2~4p' > reads.fasta

    echo "Run BLAST to identify contamination"
    blastn -db DB/pacbio_vectors_db -query reads.fasta -num_threads ${task.cpus} -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust no -soft_masking true -evalue 700 -searchsp 1750000000000 -outfmt 6 > reads.contaminant.blastout

    echo "Create blocklist of contaminated reads"
    cat reads.contaminant.blastout | grep 'NGB0097' | awk -v OFS='\t' -v var1="${adapterlength}" -v var2="${pctmatch}" '{if ((\$2 ~ /NGB00972/ && \$3 >= var2 && \$4 >= var1) || (\$2 ~ /NGB00973/ && \$3 >= 97 && \$4 >= 34)) print \$1}' | sort -u > reads.blocklist

    echo "Remove contaminated reads from FASTQ"
    cat ${reads} | paste - - - - | grep -v -f reads.blocklist -F | tr "\t" "\n" > reads.filt.fastq

    # Calculate and write statistics
    f=\$(wc -l < reads.blocklist) # Number of adapter contaminated reads
    r1=\$(cat ${reads} | wc -l) 
    r2=\$(awk -v r1=\$r1 'BEGIN{ans=r1/4; print ans}') # Number of reads
    p1=\$(awk -v n1=\$r2 -v n2=\$f 'BEGIN{ans=n2/n1*100; print ans}') # Proportion of adapter contaminated reads
    r3=\$(awk -v r2=\$r2 -v f=\$f 'BEGIN{ans=r2-f; print ans}') # Number of reads retained
    p2=\$(awk -v p1=\$p1 'BEGIN{ans=100-p1; print ans}') # Proportion of reads retained

    echo "For the ${reads} dataset:" > stats.txt
    echo "Removing reads containing adapters a minimum of ${adapterlength} bp in length and ${pctmatch}% match." >stats.txt
    echo "Number of reads:" \$r2 > stats.txt
    echo "Number of adapter contaminated reads:" \$f "(\$p1% of total)" > stats.txt
    echo "Number of reads retained:" \$r3 "(\$p2% of total)"  > stats.txt
    """

}

process fastqcReads{

    conda 'fastqc'
    publishDir 'results/qc/fastqc'

    input:
        path reads
    output:
        path '*.html'
    script:
    """
    fastqc ${reads}
    """

}

process countKmers{

    conda 'bioconda::kmer-jellyfish=2.3.1'

    input:
        path reads
    output:
        path 'reads.histo'
    script:
    """
    jellyfish count -C -m ${params.kmerLen} -s 1000000000 -t ${task.cpus} ${reads} -o reads.jf
    jellyfish histo -t ${task.cpus} reads.jf > reads.histo
    """
}

process genomeScope{

    conda 'r-base r-essentials r-ggplot2 r-data.table'
    publishDir 'results/qc/genomeScope'

    input:
        path histogram
        path genomeScopeScript
    output:
        path 'results/*'
    script:
    """
    Rscript ${genomeScopeScript} ${histogram} ${params.kmerLen} ${params.readLen} results
    """
}