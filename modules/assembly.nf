workflow{
    

}

process hifiasmPrimary{
    
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

process gfaToFasta{

    publishDir 'results/assemblyFasta'

    input:
        path p_assembly
        path a_assembly
    output:
        path "LSHT.p_ctg.fa"
        path "LSHT.a_ctg.fa"
    script:
    """
    awk '/^S/{print ">"$2;print $3}' ${p_assembly} > LSHT.p_ctg.fa
    awk '/^S/{print ">"$2;print $3}' ${a_assembly} > LSHT.a_ctg.fa
    """
}

process quastOnAssemblies{

    conda 'quast'

    input:
        path p_assembly
        path a_assembly
    output:
        path "quast_results/report.txt"
    script:
    """
    quast ${p_assembly} $
    """
}

process busco{
    conda 'busco=5.7.1'

    input:
        path p_assembly
        path a_assembly
    output:
        path "BUSCO_*"
    script:
    """
    busco -i ${p_assembly} -m genome -l euglenozoa_odb10 -c ${task.cpus} -o p_results --metaeuk
    busco -i a_assembly} -m genome -l euglenozoa_odb10 -c ${task.cpus} -o a_results --metaeuk
    """
}

process pre_merqury{

    conda 'meryl'
    input:
        path p_assembly
        path a_assembly
    output:
        path "merged_db*"
    script:
    """
    meryl -B -s ${p_assembly} -o p_meryl -m 31
    meryl -B -s ${a_assembly} -o a_meryl -m 31

    meryl -M add -s p_meryl -s a_meryl -o merged_db
    """

}

process merqury{
    conda 'merqury'

    input:  
        path p_assembly
        path a_assembly
        path db
    output:
        path "merqury*"
    script:
    """
    merqury.sh ${db} ${p_assembly} ${a_assembly} merqury
    """
}

