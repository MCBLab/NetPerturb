/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { GENIE3 } from "./modules/local/genie3/main.nf"
include { SCTENIFOLDKNK_KO } from "./modules/local/sctenifoldknk_ko/main.nf"
include { SCTENIFOLDKNK_BUILD } from "./modules/local/sctenifoldknk_build/main.nf"
include { SCTENIFOLDKNK_NULL } from "./modules/local/sctenifoldknk_null/main.nf"
include { SCRANK } from "./modules/local/scrank/main.nf"
include { HDWGCNA } from "./modules/local/hdwgcna/main.nf"
include { DOWNSAMPLE } from "./modules/local/downsample_and_split/main.nf"
include { RANK_SCORE } from "./modules/local/rank_score/main.nf"
include { MERGE } from "./modules/local/merge/main.nf"
include { MERGE_SCTENIFOLDKNK } from "./modules/local/merge_sctenifoldknk/main.nf"
include { GSEA_SCTENIFOLDKNK } from "./modules/local/gsea_sctenifoldknk/main.nf"
include { REPORT } from "./modules/local/report/main.nf"
include { NETWORK_QC } from "./modules/local/network_qc/main.nf"
include { EXTRA_TARGET_QC as EXTRA_TARGET_QC_RANK } from "./modules/local/extra_target_qc/main.nf"
include { EXTRA_TARGET_QC as EXTRA_TARGET_QC_KNOCKOUT } from "./modules/local/extra_target_qc/main.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    obj = file(params.obj)
    column = params.column
    species = params.species
    n_cells = params.n_cells
    n_cores = params.n_cores
    target = file(params.target)
    network = params.network
    sctknk = params.sctknk

    // The two tracks are independent. --network picks at most one rank-score
    // method and runs it through RANK_SCORE/MERGE/REPORT; --sctknk adds the
    // scTenifoldKnk track beside it. Either runs on its own, both run in
    // parallel when both are asked for, but at least one has to be asked for.
    if( network != null && !(network in ['genie3', 'scrank', 'hdwgcna']) ) {
        error "Invalid --network '${network}'. Supported values: genie3, scrank or hdwgcna"
    }

    if( network == null && !sctknk ) {
        error "Nothing to run. Pass --network (genie3, scrank or hdwgcna), --sctknk, or both."
    }

    // Enrichment runs on the knockout table, so without that track there is
    // nothing for a GMT to enrich. Warned rather than fatal: every other part
    // of the run is still valid.
    if( params.gsea_gmt && !sctknk ) {
        log.warn "--gsea_gmt was given without --sctknk; there is no knockout table to enrich, so no GSEA will run."
    }

    // --extra_target adds targets to a run whose networks are already built,
    // without building any again. The file is read only by EXTRA_TARGET_QC,
    // never by DOWNSAMPLE or a network step, so on -resume everything up to
    // the networks is reused, and so is every RANK_SCORE and knockout of a
    // --target: only the extras' own tasks run, then the merges and REPORT.
    // A gene the networks were not built on cannot be added this way, so it
    // is not run, and the report says so; adding it to --target instead
    // rebuilds the run around it.
    extra_target = params.extra_target ? file(params.extra_target, checkIfExists: true) : null

    // --batch is optional, but a process input cannot be null, so an unset one
    // is passed as an empty string, which DOWNSAMPLE reads as no batch.
    DOWNSAMPLE( obj, target, column, species, n_cells, params.min_cells, params.assay, params.n_hvg, params.seed, params.sctknk_min_pct, params.batch ?: "" )

    // Everything past DOWNSAMPLE reads the targets that passed its QC check
    // (present in the object, and with counts in the retained cells) rather
    // than --target itself; the ones that failed go to the report instead.
    target_qc_file = DOWNSAMPLE.out.targets

    // one target per line, ';' joining genes perturbed together
    target_ch = target_qc_file
        .splitText()
        .map { it.trim() }
        .filter { it }

    // The knockout track takes the same targets RANK_SCORE does, line for
    // line: a ';'-joined target is one joint knockout of all its genes, not
    // one knockout per gene.
    sctknk_target_ch = target_ch.unique()

    DOWNSAMPLE.out.scrank_obj
    .flatten()
    .set { sc_obj }

    // --cell_subset narrows the run to some of the identities: a
    // comma-separated list of --column values, matched against the split
    // objects DOWNSAMPLE wrote. It names each one after its identity with every
    // character outside [A-Za-z0-9_-] turned into '_', so the requested names
    // get the same treatment and either spelling matches. Both tracks read the
    // narrowed set -- the networks --network builds, and so RANK_SCORE's
    // scores, as well as the knockout track. The filter sits here, on the
    // split objects, rather than in DOWNSAMPLE, so a run resumed with a
    // different subset reuses the DOWNSAMPLE it already has, and the networks
    // of every identity it built before.
    analysis_obj = sc_obj
    if( params.cell_subset ) {
        def raw = params.cell_subset
        def wanted = ( raw instanceof Collection ? raw.collect { it.toString() } : raw.toString().split(',') as List )
            .collect { it.trim().replaceAll(/[^A-Za-z0-9_\-]/, '_') }
            .findAll { it }
            .unique()
        if( !wanted ) {
            error "--cell_subset '${raw}' names no identity."
        }
        analysis_obj = sc_obj.filter { it.baseName in wanted }

        // A name with no split object -- misspelt, or an identity DOWNSAMPLE
        // dropped for --min_cells -- is reported rather than skipped without
        // a word, and the run stops when none matched: nothing would be built,
        // and with --network on REPORT would silently never run.
        sc_obj.map { it.baseName }.collect().subscribe { found ->
            def missing = wanted - found
            if( missing.size() == wanted.size() ) {
                error "No identity matches --cell_subset '${raw}'. DOWNSAMPLE kept: ${found.sort().join(', ')}."
            }
            if( missing ) {
                log.warn "--cell_subset names ${missing.join(', ')}, which DOWNSAMPLE did not keep (misspelt, or below --min_cells); the run covers ${(wanted - missing).join(', ')} only."
            }
        }
    }

    // scTenifoldKnk in two steps. The wild-type network does not depend on
    // the target, so SCTENIFOLDKNK_BUILD makes it once per cell type; the
    // knockout on it is target-specific, so SCTENIFOLDKNK_KO runs once per
    // (cell type, target) pair, in parallel. Its table never enters
    // RANK_SCORE -- it goes to its own merge, and from there into REPORT when
    // --network is running too.
    if( sctknk ) {
        SCTENIFOLDKNK_BUILD( analysis_obj, n_cores, params.sctknk_min_pct, params.seed, params.sctknk_td_k )

        // The null model: the same knockout run for --sctknk_null random
        // genes of each cell type's network, once per cell type, which every
        // knockout of that cell type is then tested against. With
        // --sctknk_null 0 the sentinel stands in and the knockout tables
        // carry scTenifoldKnk's own statistic only. Keyed by cell type to
        // meet the network it was built from.
        sctknk_wt_keyed = SCTENIFOLDKNK_BUILD.out.wt
            .map { wt -> [ wt.name.replace('_sctknk_wt.rds', ''), wt ] }
        if( params.sctknk_null > 0 ) {
            SCTENIFOLDKNK_NULL( SCTENIFOLDKNK_BUILD.out.wt, target_qc_file, params.sctknk_null, n_cores, params.seed, params.sctknk_ndim )
            sctknk_null_keyed = SCTENIFOLDKNK_NULL.out.null_model
                .map { nm -> [ nm.name.replace('_sctknk_null.rds', ''), nm ] }
        } else {
            sctknk_null_keyed = sctknk_wt_keyed
                .map { ct, wt -> [ ct, file("${projectDir}/assets/NO_SCTKNK_NULL") ] }
        }

        // The extras the knockout networks carry join the run's targets; a
        // ';'-joined one reduced to a --target it equals is dropped upstream,
        // and an exact repeat here.
        if( extra_target ) {
            EXTRA_TARGET_QC_KNOCKOUT( SCTENIFOLDKNK_BUILD.out.wt.collect(), extra_target, target, 'knockout', 'sctenifoldknk' )
            sctknk_target_ch = sctknk_target_ch
                .mix( EXTRA_TARGET_QC_KNOCKOUT.out.targets.splitText().map { it.trim() }.filter { it } )
                .unique()
        }

        sctknk_wt_keyed
        .join( sctknk_null_keyed )
        .combine( sctknk_target_ch )
        .map { ct, wt, nm, target -> tuple( wt, nm, target ) }
        .set { sctknk_input }

        SCTENIFOLDKNK_KO( sctknk_input, n_cores, params.sctknk_plot, params.seed, params.sctknk_ndim )

        MERGE_SCTENIFOLDKNK( SCTENIFOLDKNK_KO.out.dr_table.collect(), SCTENIFOLDKNK_KO.out.status.collect() )

        // Ranks each pair's genes by how far the knockout moved them and asks
        // which gene sets sit at the top of that ranking -- the analysis the
        // scTenifoldKnk paper runs on its own output. Gene sets come from a
        // file rather than a web service so the step runs on a node with no
        // internet, which means no --gsea_gmt is simply no enrichment.
        if( params.gsea_gmt ) {
            GSEA_SCTENIFOLDKNK(
                MERGE_SCTENIFOLDKNK.out.merged_dr_table,
                file(params.gsea_gmt),
                params.gsea_min_size,
                params.gsea_max_size,
                params.seed
            )
        }
    }

    if( network ) {

        if( network == 'genie3' ) {
           GENIE3( analysis_obj, n_cores, params.seed )

            GENIE3.out.rank_obj
            .collect()
            .set { rank_cells  }
        }
        else if( network == 'scrank' ) {
            SCRANK( analysis_obj, species, target_qc_file, column, n_cores )

            SCRANK.out.rank_obj
            .collect()
            .set { rank_cells  }
        }
        else if( network == 'hdwgcna' ) {
            HDWGCNA( analysis_obj, column, n_cores, params.cut_ratio, params.hdwgcna_min_cells, params.seed )

            HDWGCNA.out.rank_obj
            .collect()
            .set { rank_cells  }
        }

        // The same for the rank-score networks: only the extras they carry
        // are scored.
        rank_target_ch = target_ch
        if( extra_target ) {
            EXTRA_TARGET_QC_RANK( rank_cells, extra_target, target, 'rank_score', network )
            rank_target_ch = target_ch
                .mix( EXTRA_TARGET_QC_RANK.out.targets.splitText().map { it.trim() }.filter { it } )
                .unique()
        }

        RANK_SCORE( obj, rank_target_ch, species, column, params.binding, params.top_connections, params.assay, params.cut_ratio, rank_cells )

        MERGE( RANK_SCORE.out.rank_scores.collect(), RANK_SCORE.out.top_connections.collect() )

        // Structure of every network this run built, for the report's Data
        // quality section: the rank-score networks, hdWGCNA's metacell and fit
        // numbers when that was the method, and the knockout networks when
        // that track ran too. Summarised once here, since REPORT only ever
        // receives tables.
        qc_networks = rank_cells
        if( network == 'hdwgcna' ) {
            qc_networks = qc_networks.mix( HDWGCNA.out.qc.collect() )
        }
        if( sctknk ) {
            qc_networks = qc_networks.mix( SCTENIFOLDKNK_BUILD.out.wt.collect() )
        }
        // The extras that were run are placed in the networks too. Without
        // any, NETWORK_QC reads targets_qc.txt itself, as it always has.
        network_qc_targets = target_qc_file
        if( extra_target ) {
            network_qc_targets = target_qc_file
                .mix( EXTRA_TARGET_QC_RANK.out.targets )
                .mix( sctknk ? EXTRA_TARGET_QC_KNOCKOUT.out.targets : Channel.empty() )
                .collectFile( name: "targets_network_qc.txt" )
        }
        NETWORK_QC( qc_networks.collect(), network_qc_targets, network )

        // REPORT always takes a scTenifoldKnk table path; when it was not
        // requested this is a sentinel empty file report.qmd recognises and
        // renders as "not run for this session" rather than a real table.
        sctknk_table = sctknk
            ? MERGE_SCTENIFOLDKNK.out.merged_dr_table
            : file("${projectDir}/assets/NO_SCTKNK_TABLE")

        // And for the knockout status table, one line per cell type x target
        // saying whether the knockout ran and how far it moved the network.
        sctknk_status = sctknk
            ? MERGE_SCTENIFOLDKNK.out.status
            : file("${projectDir}/assets/NO_SCTKNK_STATUS")

        // Same sentinel arrangement for the enrichment table, which has two
        // ways of not existing: the knockout track was not run at all, or it
        // was run without a --gsea_gmt to enrich against.
        gsea_table = (sctknk && params.gsea_gmt)
            ? GSEA_SCTENIFOLDKNK.out.gsea_table
            : file("${projectDir}/assets/NO_GSEA_TABLE")

        // Which extra targets each track ran, and which it could not. The
        // sentinel is a header-only table, read as no --extra_target.
        extra_target_qc = extra_target
            ? EXTRA_TARGET_QC_RANK.out.qc
                .mix( sctknk ? EXTRA_TARGET_QC_KNOCKOUT.out.qc : Channel.empty() )
                .collectFile( name: "extra_target_qc.tsv", keepHeader: true )
            : file("${projectDir}/assets/NO_EXTRA_TARGET_QC")

        // What the report says about the run itself: when it started, the
        // command line, and the settings that shape the result. Written as a
        // key/value table since REPORT's container cannot see the workflow
        // object. The epoch lets the report measure elapsed time without
        // caring which time zone the container runs in.
        def start = workflow.start.toInstant()
        def run_info = [
            "run_name"        : workflow.runName,
            "session_id"      : workflow.sessionId,
            "started"         : java.time.ZonedDateTime.ofInstant(start, java.time.ZoneId.systemDefault())
                                    .format(java.time.format.DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss Z")),
            "started_epoch_ms": start.toEpochMilli(),
            "command_line"    : workflow.commandLine,
            "resume"          : workflow.resume,
            "profile"         : workflow.profile,
            "container_engine": workflow.containerEngine ?: "none",
            "nextflow_version": workflow.nextflow.version,
            "pipeline_version": workflow.manifest.version,
            "revision"        : workflow.revision ?: workflow.commitId ?: "local",
            "user"            : workflow.userName,
            "launch_dir"      : workflow.launchDir,
            "work_dir"        : workflow.workDir,
            "obj"             : params.obj,
            "target"          : params.target,
            "column"          : params.column,
            "species"         : params.species,
            "assay"           : params.assay,
            "network"         : params.network,
            "sctknk"          : params.sctknk,
            "n_cells"         : params.n_cells,
            "min_cells"       : params.min_cells,
            "n_hvg"           : params.n_hvg,
            "cut_ratio"       : params.cut_ratio,
            "n_cores"         : params.n_cores,
            "seed"            : params.seed,
            "sctknk_min_pct"  : params.sctknk_min_pct,
            "sctknk_null"     : params.sctknk_null,
            "sctknk_ndim"     : params.sctknk_ndim,
            "sctknk_td_k"     : params.sctknk_td_k,
            "cell_subset"     : params.cell_subset ?: "all",
            "extra_target"    : params.extra_target ?: "none",
            "batch"           : params.batch
        ]
        run_info_file = Channel
            .of( run_info.collect { k, v -> "${k}\t${String.valueOf(v).replaceAll(/[\t\r\n]+/, ' ')}" }.join("\n") + "\n" )
            .collectFile( name: "run_info.tsv" )

        // Per-step running times come from the execution trace, which
        // Nextflow appends to one row per task as each one finishes. It is
        // read once everything REPORT depends on has finished -- which is
        // every task of the run bar REPORT itself -- and copied into the work
        // dir, since the container cannot reach the published one. Rows are
        // written asynchronously right after a task ends, hence the short
        // wait before reading. The path is the one Nextflow is actually
        // writing to, so a -with-trace <file> or a trace.file from -c is
        // followed; the nextflow.config default is only the fallback. A
        // missing trace (trace disabled) becomes a one-line note the report
        // says it lacks -- not an empty string, which collectFile would emit
        // no file for, leaving REPORT waiting forever.
        def trace_cfg = workflow.session.config.navigate('trace.file')
        def trace_path = file( trace_cfg instanceof CharSequence && trace_cfg
            ? trace_cfg.toString()
            : "${params.tracedir}/execution_trace_${params.trace_report_suffix}.txt" )
        report_upstream = MERGE.out.merged_rank_scores.mix( NETWORK_QC.out.network_qc )
        if( sctknk ) {
            report_upstream = report_upstream.mix( MERGE_SCTENIFOLDKNK.out.merged_dr_table )
        }
        if( sctknk && params.gsea_gmt ) {
            report_upstream = report_upstream.mix( GSEA_SCTENIFOLDKNK.out.gsea_table )
        }
        trace_snapshot = report_upstream
            .collect()
            .map { sleep(3000); trace_path.exists() ? trace_path.text : "# no trace file at ${trace_path}\n" }
            .collectFile( name: "trace_snapshot.tsv" )

        REPORT(
            MERGE.out.merged_rank_scores,
            DOWNSAMPLE.out.umap,
            file("${projectDir}/bin/report.qmd"),
            network,
            params.score_quantile,
            sctknk_table,
            params.sctknk_top_genes,
            gsea_table,
            params.gsea_top_terms,
            DOWNSAMPLE.out.target_qc,
            DOWNSAMPLE.out.cell_counts,
            DOWNSAMPLE.out.target_expression,
            run_info_file,
            trace_snapshot,
            DOWNSAMPLE.out.identity_qc,
            NETWORK_QC.out.network_qc,
            NETWORK_QC.out.target_network_qc,
            sctknk_status,
            extra_target_qc
        )
    }
}
