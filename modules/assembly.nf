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

process gfaToFasta{

    conda 'gfastats'

    input:
        path p_assembly
        path a_assembly
    output:
        path 
    script:
    """
    gfastats ${p_assembly} -o fa
    gfastats ${a_assembly} -o fa
    """
}

process gfaStats{

    conda 'gfastats'

    input:
        path p_assembly
        path a_assembly
    output:
        path "gfastats_on_pa_assemblies.txt"
    script:
    """
    gfastats summary -i primary_contigs.gfa -g 11747160 -o primary_stats.txt
    gfastats summary -i alternate_contigs.gfa -g 11747160 -o alternate_stats.txt

    # Join the two output files
    paste primary_stats.txt alternate_stats.txt > joined_stats.txt

    # Filter out lines containing "scaffold"
    grep -vi "scaffold" joined_stats.txt > gfastats_on_pa_assemblies.txt
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
    busco -i p_assembly -m genome -l euglenozoa_odb10 -c ${task.cpus} -o p_results --metaeuk
    busco -i a_assembly -m genome -l euglenozoa_odb10 -c ${task.cpus} -o a_results --metaeuk
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

