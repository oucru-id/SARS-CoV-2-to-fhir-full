#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process SNPEFF {
    tag "$meta.id"
    publishDir "${params.results_dir}/variants", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.annotated_variants.vcf.gz"), emit: vcf
    tuple val(meta), path("${meta.id}.variants.tsv"),              emit: table
    path "${meta.id}.snpeff.csv",                                  emit: csv
    path "versions.yml",                                           emit: versions

    script:
    """
    if [[ "${vcf}" == *.gz ]]; then
        gunzip -c ${vcf} > input.vcf
    else
        cp ${vcf} input.vcf
    fi

    snpEff \\
        -config ${params.snpeff_config} \\
        -dataDir ${params.snpeff_db} \\
        -csvStats ${meta.id}.snpeff.csv \\
        -noLog \\
        ${params.contig} \\
        input.vcf > ${meta.id}.annotated_variants.vcf

    bgzip -f ${meta.id}.annotated_variants.vcf
    tabix -f -p vcf ${meta.id}.annotated_variants.vcf.gz

    snpSift extractFields -s "," -e "." \\
        ${meta.id}.annotated_variants.vcf.gz \\
        CHROM POS REF ALT FILTER DP AF \\
        "ANN[0].GENE" "ANN[0].EFFECT" "ANN[0].IMPACT" \\
        "ANN[0].HGVS_C" "ANN[0].HGVS_P" \\
        > ${meta.id}.variants.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(snpEff -version 2>&1 | sed 's/SnpEff\\s*//' | cut -f1)
    END_VERSIONS
    """

    stub:
    """
    echo | gzip > ${meta.id}.annotated_variants.vcf.gz
    touch ${meta.id}.variants.tsv ${meta.id}.snpeff.csv versions.yml
    """
}

workflow ANNOTATE {
    take:
    vcfs

    main:
    SNPEFF(vcfs)

    emit:
    vcf      = SNPEFF.out.vcf
    table    = SNPEFF.out.table
    csv      = SNPEFF.out.csv
    versions = SNPEFF.out.versions
}