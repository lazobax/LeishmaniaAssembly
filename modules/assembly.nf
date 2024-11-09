workflow{
    

}

process hifiasmDefault{
    
    conda 'hifiasm'
    publishDir "results/assembly/hifiasm"
    input:
        path reads
    output:
        path "LSHT.bp.hap1.p_ctg.gfa"
        path "LSHT.bp.hap2.p_ctg.gfa"
        path "LSHT.bp.p_ctg.gfa"
    script:
    """
    hifiasm -o LSHT -t ${task.cpus} ${reads} --hg-size 35m
    """
}


process gfaToFasta{

    publishDir 'results/assembly/fasta'

    input:
        path hap1_assembly
        path hap2_assembly
        path p_assembly
    output:
        path "LSHT.hap1_ctg.fa"
        path "LSHT.hap2_ctg.fa"
        path "LSHT.p_ctg.fa"
    script:
    """
    awk '/^S/{print ">"\$2;print \$3}' ${p_assembly} > LSHT.p_ctg.fa
    awk '/^S/{print ">"\$2;print \$3}' ${hap1_assembly} > LSHT.hap1_ctg.fa
    awk '/^S/{print ">"\$2;print \$3}' ${hap2_assembly} > LSHT.hap2_ctg.fa
    """
}