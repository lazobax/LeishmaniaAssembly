
params.shortReadAccession = "SRR12368184"
params.s3bucket = "s3://ltropica-data/Lei012024_HiFi.fastq.gz"

workflow {
    short_reads = download_short_reads(params.shortReadAccession)
    // // long_reads = download_long_reads(params.s3bucket)
    
    // flye_assembly = Flye_long_read_assembly(long_reads)

    // map_short_reads_minimap2(short_reads, flye_assembly)
    // index_assembly_samtools(flye_assembly)
    // view_short_reads_samtools(map_short_reads_minimap2.out, index_assembly_samtools.out)
    // sort_short_reads_samtools(view_short_reads_samtools.out)
    // index_short_reads_samtools(sort_short_reads_samtools.out)

    // mpileUP_assembly_bcftools(sort_short_reads_samtools.out, flye_assembly)
    // norm_assembly_bcftools(mpileUP_assembly_bcftools.out, flye_assembly)
    // filter_assembly_bcftools(norm_assembly_bcftools.out)
    // index_assembly_bcftools(filter_assembly_bcftools.out)
    // consensus_assembly_bcftools(index_assembly_bcftools.out, filter_assembly_bcftools.out, flye_assembly)

    // polish_assembly_pilon(index_short_reads_samtools.out, sort_short_reads_samtools.out, flye_assembly)
    // index_polished_assembly_samtools(polish_assembly_pilon.out)
    // map_polished_assembly_minimap2(short_reads, polish_assembly_pilon.out)
    // view_polished_assembly_samtools(map_polished_assembly_minimap2.out, polish_assembly_pilon.out, index_polished_assembly_samtools.out)
    // sort_polished_assembly_samtools(view_polished_assembly_samtools.out)
    // mpileup_polished_assembly_bcftools(sort_polished_assembly_samtools.out, polish_assembly_pilon.out)
    // norm_polished_assembly_bcftools(mpileup_polished_assembly_bcftools.out, polish_assembly_pilon.out)
    // filter_polished_assembly_bcftools(norm_polished_assembly_bcftools.out)
    // index_polished_assembly_bcftools(filter_polished_assembly_bcftools.out)
    // consensus_polished_assembly_bcftools(index_polished_assembly_bcftools.out, filter_polished_assembly_bcftools.out, polish_assembly_pilon.out)
    
    // aggregate_reads(short_reads)

    // sort_polished_assembly_funannotate(consensus_polished_assembly_bcftools.out)

    // download_reference_genome_ragoo()

    // arrange_assembly_ragoo(aggregate_reads.out, sort_polished_assembly_funannotate.out, download_reference_genome_ragoo.out)
    // clean_polished_assembly_funannotate(arrange_assembly_ragoo.out)
    // final_assembly(clean_polished_assembly_funannotate.out)

    // qc_assemblies_quast(flye_assembly, consensus_assembly_bcftools.out, polish_assembly_pilon.out, consensus_polished_assembly_bcftools.out, arrange_assembly_ragoo.out, clean_polished_assembly_funannotate.out, final_assembly.out)
}

process download_short_reads {
    input:
        val accession
    output:
        path("${accession}_{1,2}.fastq")
    conda 'sra-tools=2.10.9'
    script:
    """
    prefetch ${accession}
    fasterq-dump --split-files ${accession} 
    """
}

process download_long_reads{
    input:
        "s3"
    output:
        "longRead"
    script:
    '''
    mv s3 longRead
    '''
}


process Flye_long_read_assembly{

    input:
        path "longRead.fastq.gz"
    output:
        path "assembly.fasta"
    conda "flye=2.8.2"
    script:
    '''
    flye --pacbio-hifi ${longRead} --genome-size 35m --threads 22
    '''


}

process index_assembly_samtools{

    input:
        path "assembly.fasta"
    output:
        path "assembly.fasta.fai"
    conda 'samtools=1.11'
    script:
    '''
    samtools faidx assembly.fasta > assembly.fasta.fai
    '''

}

process map_short_reads_minimap2{

    input:
        path 'read_'
        path 'assembly'
    output:
        path "PE.sam"
    conda 'minimap2=2.17'
    script:
    '''
    minimap2 -t 22 -ax sr assembly read_1 read_2 > PE.sam
    '''

}

process view_short_reads_samtools{

    input:
        path 'PE.sam'
        path 'index'
    output:
        path "PE.bam"
    conda 'samtools=1.11'
    script:
    '''
    samtools view -@ 22 index PE.sam > PE.bam
    '''

}

process sort_short_reads_samtools{
    input:
        path "PE"
    output:
        path "PE_sort"
    conda 'samtools=1.11'
    script:
    '''
    samtools sort -@ 22 PE -o PE_sort
    '''
}

process index_short_reads_samtools{
    input:
        path "PE_sort"
    output:
        path "PE_index"
    conda 'samtools=1.11'
    script:
    '''
    samtools index -@ 22 PE_sort > PE_index
    '''
}

process mpileUP_assembly_bcftools{
    input:
        path "PE_sort"
        path "assembly"
    output:
        path "PE.vcf.gz"
    conda 'bcftools=1.11'
    script:
    '''
    bcftools mpileup --threads 22 -Ou -d 1000000 -f assembly PE_sort | bcftools call -mv -Oz -o PE.vcf.gz"
    '''
}

process norm_assembly_bcftools{
    input:
        path "PE.vcf.gz"
        path "assembly"
    output:
        path "PE.norm.bcf"
    conda 'bcftools=1.11'
    script:
    '''
    bcftools norm --threads 22 -f assembly PE.vcf.gz -Ob -o PE.norm.bcf
    '''

}

process filter_assembly_bcftools{

    input:
        path "PE.norm"
    output:
        path "PE.norm.flt-indels.bcf"
    conda 'bcftools=1.11'
    script:
    '''
    bcftools filter --threads 22 --IndelGap 5 PE.norm -Ob -o PE.norm.flt-indels.bcf
    '''
    
}

process index_assembly_bcftools{

    input:
        path "PE.norm.indels.bcf"
    output:
        path "PE.norm.indels.csi"
    conda 'bcftools=1.11'
    script:
    '''
    bcftools index --threads 22 PE.norm.indels.bcf -o PE.norm.indels.csi
    '''

}

process consensus_assembly_bcftools{

    input:
        path "PE.csi"
        path "PE.bcf"
        path "assembly"
    output:
        path "Flye_consensus.fasta"
    conda 'bcftools=1.11'
    script:
    '''
    bcftools onsensus -f assembly PE.bcf > Flye_consensus.fasta
    '''
}

process polish_assembly_pilon{
    input:
        path "PEindex.bam"
        path "PEsort.bam"
        path "assembly"
    output:
        path "pilon.fasta"
    conda 'pilon=1.23'
    script:
    '''
    pilon --threads 22 -Xmx300000M --genome assembly --bam PEsort.bam 
    '''
}

process index_polished_assembly_samtools{
    input:
        path "pilon.fasta"
    output:
        path "pilon.fai"
    conda 'samtools=1.11'
    script:
    '''
    samtools faidx pilon.fasta > pilon.fai
    '''
}

process map_polished_assembly_minimap2{

    input:
        path "read_"
        path "pilon.fasta"
    output:
        path "PE.sam"
    conda 'minimap2=2.17'
    script:
    '''
    minimap2 -t 22 -ax sr pilon.fasta read_1 read_2 > PE.sam
    '''
}

process view_polished_assembly_samtools{
    input:
        path "PE.sam"
        path "polishedAssembly"
        path "polishedAssIx"
    output:
        path "PE.bam"
    conda 'samtools=1.11'
    script:
    '''
    samtools view -@ 22 -bt polishedAssIx PE.sam > PE.bam
    '''
}

process sort_polished_assembly_samtools{
    input: 
        path "PE.bam"
    output:
        path "PE_sort.bam"
    conda 'samtools=1.11'
    script:
    '''
    samtools sort -@ 22 PE.bam -o PE_sort.bam
    '''
}

process mpileup_polished_assembly_bcftools{
    input:
        path "PE_sort.bam"
        path "polishedAssembly"
    output:
        path "PE.vcf.gz"
    conda 'bcftools=1.11'
    script:
    '''
    bcftools mpileup --threads 22 -Ou -d 1000000 -f polishedAssembly PE_sort.bam | bcftools call -mv -Oz -o PE.vcf.gz
    '''
}

process norm_polished_assembly_bcftools{
    input:
        path "PE.vcf.gz"
        path "polishedAssembly"
    output:
        path "PE.norm.bcf"
    conda 'bcftools=1.11'
    script:
    '''
    bcftools norm --threads 22 -f polishedAssembly PE.vcf.gz -Ob -o PE.norm.inde.bcf
    '''
}

process filter_polished_assembly_bcftools{
    input:
        path "PE.norm.bcf"
    output:
        path "PE.norm.indels.bcf"
    conda 'bcftools=1.11'
    script:
    '''
    bcftools filter --threads 22 --IndelGap 5 PE.norm.bcf -Ob -o PE.norm.indels.bcf
    '''
}

process index_polished_assembly_bcftools{
    input:
        path "PE.norm.indels.bcf"
    output:
        path "PE.norm.indels.bcf.csi"
    conda 'bcftools=1.11'
    script:
    '''
    bcftools index --threads 22 PE.norm.indels.bcf -o PE.norm.indels.bcf.csi
    '''
}

process consensus_polished_assembly_bcftools{
    input:
        path "PE.norm.indels.bcf.csi"
        path "PE.norm.indels.bcf"
        path "polishedAssembly"
    output:
        path "pilonConsensus.fasta"
    conda 'bcftools=1.11'
    script:
    '''
    bcftools consensus -f polishedAssembly PE.norm.indels.bcf > pilonConsensus.fasta
    '''
}

process aggregate_reads{
    input:
        path "read_"
    output:
        path "short_reads"
    script:
    '''
    cat read_1 read_2 > short_reads
    '''
}

process sort_polished_assembly_funannotate{
    input:
        path "pilonConsensus"
    output:
        path "SortedPolished"
    container 'nextgenusfs/funannotate'
    script:
    '''
    funannotate sort -i pilonConsensus -o SortedPolished
    '''
}

process download_reference_genome_ragoo{

    output:
        path "reference"
    script:
    '''
    wget -O reference https://tritrypdb.org/common/downloads/release-47/LmajorFriedlin/fasta/data/TriTrypDB-47_LmajorFriedlin_Genome.fasta
    '''
}

process arrange_assembly_ragoo{
    input:
        path "shortReads"
        path "polishedAssembly"
        path "reference"
    output:
        path "ragoo_output/ragoo.fasta"
    conda 'ragoo=1.1'
    script:
    '''
    ragoo.py -t 22 -C -T sr -b -g 0 -R shortReads polishedAssembly reference
    '''
}

process clean_polished_assembly_funannotate{
    input:
        path "ragoo.fasta"
    output:
        path "funannotate_assembly.fasta"
    container 'nextgenusfs/funannotate'
    script:
    '''
    funannotate clean --exhaustive --minlen 50 --pident 10 -i ragoo.fasta -o  funannotate_assembly.fasta
    '''
}

process final_assembly{

    publishDir 'Final_Assembly/'

    input:
        path "funannotate"
    output:
        path "genome.fasta"
    script:
    '''
    cp funannotate genome.fasta; sed -i 's/_RaGOO//g' genome.fasta; sed -i 's/LmjF/LSHT/g' genome.fasta"
    '''
}

process qc_assemblies_quast{

    publishDir 'qc/'

    input:
        path "flye"
        path "flyeConsensus"
        path "pilon"
        path "pilonConsesus"
        path "ragoo.fasta"
        path "funannotate"
        path "genome"
    output:
        path "report.pdf"
    container 'quay.io/biocontainers/quast:5.0.1--py27pl526ha92aebf_0'
    script:
    '''
    quast.py --eukaryote --circos --threads 22 --output-dir . flye flyeConsensus pilon pilonConsensus ragoo.fasta funannotate genome
    '''
}


println "hi"
