process DOWNSAMPLE {
  """
  Downsamples and splits the Seurat object into cell types
  """

  label 'process_medium'
  label "r_scrank"

  container "${ workflow.containerEngine == 'singularity' ? 'docker://juliaapolonio/scrank:latest':
            'docker.io/juliaapolonio/scrank:latest' }"

  input:
    path obj
    path target
    val column
    val species
    val n_cells
    val min_cells
    val assay
    val n_hvg
    val seed
    val sctknk_min_pct
    val batch

  output:
    path "*.RDS", emit: scrank_obj
    path "*.png", emit: umap
    path "targets_qc.txt", emit: targets
    path "target_qc.tsv", emit: target_qc
    path "cell_counts.tsv", emit: cell_counts
    path "target_expression.tsv", emit: target_expression
    path "identity_qc.tsv", emit: identity_qc

  when:
  task.ext.when == null || task.ext.when  

  script:
    """
    #!/bin/bash

    downsample_and_split.R ${obj} ${target} ${column} ${species} ${n_cells} ${min_cells} ${assay} ${n_hvg} ${seed} ${sctknk_min_pct} ${batch ?: 'null'}
    """

  stub:
    """
    touch sensitive.RDS
    touch resistant.RDS
    touch umap.png
    grep -v '^[[:space:]]*\$' ${target} > targets_qc.txt
    printf 'target\\tgene\\treason\\taction\\n' > target_qc.tsv
    printf 'identity\\tn_input\\tn_used\\tstatus\\tgenes_total\\tgenes_gene4use\\tgenes_expressed\\tgenes_filtered\\tgenes_sctknk\\tgenes_mt_rb\\n' > cell_counts.tsv
    printf 'sensitive\\t100\\t100\\tkept\\t500\\t200\\t180\\t20\\t320\\t30\\n' >> cell_counts.tsv
    printf 'resistant\\t100\\t100\\tkept\\t500\\t200\\t190\\t10\\t340\\t30\\n' >> cell_counts.tsv
    printf 'identity\\tgene\\tavg_expression\\tpct_expressing\\tvar_standardized\\tvar_percentile\\n' > target_expression.tsv
    printf 'sensitive\\tBrd4\\t1.2\\t60\\t1.4\\t80\\n' >> target_expression.tsv
    printf 'resistant\\tBrd4\\t0.8\\t45\\t1.1\\t60\\n' >> target_expression.tsv
    printf 'identity\\tidentity_id\\tn_cells\\tmedian_counts\\tmedian_genes\\tmedian_mt_rb_frac\\tsparsity_gene4use\\tpc1_var_frac\\ttop5_var_frac\\tpc1_depth_rho\\tmax_pc_depth_rho\\tmax_pc_depth\\tn_batches\\tlargest_batch_frac\\tbatch_entropy\\tbatch_r2\\n' > identity_qc.tsv
    """
}

