nextflow.enable.dsl=2

workflow{

    long_reads = Channel.fromPath('*.gz')
    flye_output = Channel.fromPath('~/LeishmaniaAssembly/work/a6/d91f443c061c6bcf91435469274b06/assembly.fasta')

    map_reads(long_reads, flye_output)
    hapdup(map_reads.out, flye_output)

}

process map_reads{

    input:
        path "reads.fastq.gz"
        path "assembly.fasta"
    output:
        path "lr_mapping.bam"
        path "lr_mapping.bam.bai"
    conda 'minimap2 samtools'
    script:
    '''
    minimap2 -ax map-ont -t 22 assembly.fasta reads.fastq.gz | samtools sort -@ 4 -m 4G > lr_mapping.bam
    samtools index -@ 22 assembly_lr_mapping.bam
    '''


}

process hapdup{
    input:
        path "lr_mapping.bam"
        path "lr_mapping.bam.bai"
        path "assembly.fasta"
    output:
        path "hapdup_dual_{1,2}.fasta" into dual_assemblies
        path "phased_blocks_hp{1,2}.bed" into phased_blocks
        path "hapdup_phased_{1,2}.fasta" into haplotype_resolved_assemblies
    docker 'mkolmogo/hapdup:0.12'
    script:
    '''
     hapdup --assembly assembly.fasta --bam lr_mapping.bam --out-dir . -t 22 --rtype hifi
    '''
}
