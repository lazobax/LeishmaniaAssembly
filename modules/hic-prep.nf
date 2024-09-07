
workflow pre_hic{
    take:
        hic_1 = Channel.fromPath("s3://sra-pub-src-17/SRR18440378/ltarentolae_hic_R1.fastq.gz.1")
        hic_2 = Channel.fromPath("s3://sra-pub-src-17/SRR18440378/ltarentolae_hic_R2.fastq.gz.1")
    

}

process get_hic_data{

    input:
        path "reads1"
        path "reads2"
    output:
        path "hic_1.fastq"
        path "hic_2.fastq"
    script:
    '''
    wget reads1
    wget reads2
    '''
}

process trim_hic_reads{

    conda 'bioconda::fastp'

    input:
        path "hic_1.fastq"
        path "hic_2.fastq"
    output:
        path "hic_1_trimmed.fastq"
        path "hic_2_trimmed.fastq"
    script:
    '''
    fastp -i hic_1.fastq -I hic_2.fastq -o hic_1_trimmed.fastq -O hic_2_trimmed.fastq
    '''

}

process qc_hic{
    
    conda 'bioconda::fastqc'
    
    input:
        path "hic_1_trimmed.fastq"
        path "hic_2_trimmed.fastq"
    output:
        path "qc_reports/*"
    script:
    '''
    fastqc -o qc_reports hic_1_trimmed.fastq hic_2_trimmed.fastq
    '''

}

process map_hic_assembly{
    conda 'bioconda::bwa bioconda::samtools'

    input:
        path "hic_1_trimmed.fastq"
        path "hic_2_trimmed.fastq"
        path "assembly.fasta"
    output:
        path "hic_map.bam"
    script:
    '''
    bwa mem -t 22 assembly hic_1_trimmed.fastq hic_2_trimmed.fastq | samtools sort -@ 22 -o hic_map.bam
    samtools index -@ 22 hic_map.bam
    '''
}

process filter_hic{
    conda 'bioconda::samtools'

    input:
        path "hic_map.bam"
    output:
        path "hic_map_filtered.bam"
    script:
    """
    samtools view -@ 22 -q 30 -b hic_map.bam > hic_map_filtered.bam
    samtools index -@ 22 hic_map_filtered.bam
    '''
}






