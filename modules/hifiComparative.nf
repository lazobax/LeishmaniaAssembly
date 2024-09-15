workflow {

    reads = channel.fromPath("")


}

process hifiasmPrimary{
    
    conda 'hifiasm'
    publishDir "results/assembly"
    input:
        path reads
    output:
        path "primary.p_ctg.gfa"
        path "primary.asm.a_ctg.gfa"
    script:
    """
    hifiasm -o primary --primary -t ${task.cpus} ${reads}
    """
}

process hifiasmDefault{
    
    conda 'hifiasm'
    publishDir "results/assembly"
    input:
        path reads
    output:
        path "default.p_ctg.gfa"
        path "default.asm.a_ctg.gfa"
    script:
    """
    hifiasm -o default -t ${task.cpus} ${reads}
    """
}