#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process NEXTCLADE {
    tag "$meta.id"
    publishDir "${params.results_dir}/nextclade", mode: 'copy'

    input:
    tuple val(meta), path(consensus)
    path dataset

    output:
    tuple val(meta), path("${meta.id}.lineage.json"), emit: lineage
    path "${meta.id}/nextclade.tsv",                  emit: tsv
    path "versions.yml",                              emit: versions

    script:
    def pangolin_arg = params.pangolin_csv ? "--pangolin_csv ${params.pangolin_csv}" : ''
    """
    mkdir -p ${meta.id}

    nextclade run \\
        --input-dataset ${dataset} \\
        --output-all ${meta.id} \\
        --jobs ${task.cpus} \\
        ${consensus}

    python3 ${projectDir}/scripts/parse_nextclade.py \\
        --nextclade-tsv ${meta.id}/nextclade.tsv \\
        --dataset ${dataset} \\
        --sample ${meta.id} \\
        --who-table ${projectDir}/data/who_variant_classification.csv \\
        --nextclade-version "\$(nextclade --version | awk '{print \$NF}')" \\
        ${pangolin_arg} \\
        --output ${meta.id}.lineage.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(nextclade --version | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}
    touch ${meta.id}/nextclade.tsv versions.yml
    echo '{}' > ${meta.id}.lineage.json
    """
}

workflow LINEAGE {
    take:
    consensus_fasta

    main:
    if (params.skip_nextclade) {
        lineage  = Channel.empty()
        tsv      = Channel.empty()
        versions = Channel.empty()
    }
    else {
        NEXTCLADE(consensus_fasta, file(params.nextclade_dataset))
        lineage  = NEXTCLADE.out.lineage
        tsv      = NEXTCLADE.out.tsv
        versions = NEXTCLADE.out.versions
    }

    emit:
    lineage  = lineage
    tsv      = tsv
    versions = versions
}
