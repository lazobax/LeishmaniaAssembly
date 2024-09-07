workflow{
    reads = channel.fromPath("results/qc/adapterFilt/reads.filt.fastq")
    hifiasm(reads)



}

process hifiasm{
    
    conda 'hifiasm'
    publishDir "results/assembly"
    input:
        path reads
    output:
        path "LSHT.asm.p_ctg.gfa"
        path "LSHT.asm.a_ctg.gfa"
    script:
    """
    hifiasm -o LSHT.asm --primary -t ${task.cpus} ${reads}
    """
}

