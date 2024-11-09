workflow {

    reads = channel.fromPath('r*fastq')

    def values = [0.5, 0.45, 0.35, 0.25, 0.15, 0.05]




}

process hifiasmDefault{
    
    tag 's value is: ${s}'
    conda 'hifiasm'
    publishDir "results/assembly"
    input:
        path reads
        val s
    output:
        path "default.bp.hap1.p_ctg.gfa"
        path "default.bp.hap2.p_ctg.gfa"
        path "default.bp.p_ctg.gfa"
    script:
    """
    hifiasm -o default -t ${task.cpus} ${reads} -s ${s}
    """
}


