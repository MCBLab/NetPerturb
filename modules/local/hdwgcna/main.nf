process HDWGCNA {
    tag "${sc_obj.baseName}" 

    publishDir "${params.outdir}/hdwgcna_networks", mode: 'copy'

    container 'docker://leoshow21/hdwgcna:v2'

    queue 'amd-512'
    cpus params.n_cores
    memory '150 GB'

    input:
    path sc_obj
    val  column
    val  n_cores

    output:
    path "*_weight_hdWGCNA_*.rds", emit: rank_obj, optional: true

    script:
    """
    hdwgcna.R "${sc_obj}" "${column}" "${n_cores}" "${params.target}"
    """
}