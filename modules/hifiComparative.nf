workflow {

    reads = channel.fromPath("*fastq")
    hifiasmPrimary(reads)
    hifiasmDefault(reads)
    quast(hifiasmPrimary.out, hifiasmDefault.out)


}

process hifiasmPrimary{
    
    conda 'hifiasm'
    publishDir "results/assembly"
    input:
        path reads
    output:
        path "primary.p_ctg.gfa"
        path "primary.a_ctg.gfa"
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
        path "default.bp.hap1.p_ctg.gfa"
        path "default.bp.hap2.p_ctg.gfa"
        path "default.bp.p_ctg.gfa"
    script:
    """
    hifiasm -o default -t ${task.cpus} ${reads}
    """
}

process quast{

    conda 'quast'
    publishDir 'results/quast'

    input:
        path pp
        path pa
        path dh1
        path dh2
        path dp

    output:
        path "quast_results/*"
    script:
    """
    quast ${pp} ${pa} ${dh1} ${dh2} ${dp}
    """


}