#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process MULTIQC {
    publishDir "${params.results_dir}/qc", mode: 'copy'

    input:
    path '*'

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data

    script:
    """
    multiqc . --force --interactive --zip-data-dir
    """

    stub:
    """
    touch multiqc_report.html
    mkdir -p multiqc_data
    """
}

process CREATE_SAMPLE_REPORT {
    tag "$meta.id"
    publishDir "${params.results_dir}/reports", mode: 'copy'

    input:
    tuple val(meta), path(annotated_vcf), path(lineage_json), path(consensus_stats)

    output:
    path "${meta.id}.summary_report.txt", emit: report

    script:
    """
    if [[ "${annotated_vcf}" == *.gz ]]; then
        gunzip -c ${annotated_vcf} > variants.vcf
    else
        cp ${annotated_vcf} variants.vcf
    fi

    python3 ${projectDir}/scripts/generate_sample_report.py \\
        --annotated_vcf variants.vcf \\
        --lineage_json ${lineage_json} \\
        --consensus_stats ${consensus_stats} \\
        --sample_id ${meta.id} \\
        --platform ${meta.platform} \\
        --protocol ${meta.protocol} \\
        --primer_scheme ${meta.scheme} \\
        --output_dir .
    """

    stub:
    "touch ${meta.id}.summary_report.txt"
}

workflow GENERATE_REPORT {
    take:
    multiqc_files

    main:
    MULTIQC(multiqc_files)

    emit:
    report = MULTIQC.out.report
}

workflow GENERATE_SAMPLE_REPORTS {
    take:
    per_sample  

    main:
    CREATE_SAMPLE_REPORT(per_sample)

    emit:
    reports = CREATE_SAMPLE_REPORT.out.report
}
