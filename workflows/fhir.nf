#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CREATE_FHIR {
    tag "$meta.id"
    publishDir "${params.results_dir}/fhir", mode: 'copy', pattern: '*.fhir.json'

    input:
    tuple val(meta), path(annotated_vcf), path(lineage_json), path(consensus_stats)
    path org_metadata

    output:
    path "${meta.id}.fhir.json", emit: fhir_output
    path "versions.yml",         emit: versions

    script:
    """
    if [[ "${annotated_vcf}" == *.gz ]]; then
        gunzip -c ${annotated_vcf} > ${meta.id}.vcf
    else
        cp ${annotated_vcf} ${meta.id}.vcf
    fi

    python3 ${projectDir}/scripts/annotated_to_fhir.py \\
        --input ${meta.id}.vcf \\
        --output ${meta.id}.fhir.json \\
        --sample_id ${meta.id} \\
        --lineage_json ${lineage_json} \\
        --consensus_stats ${consensus_stats} \\
        --organization_metadata ${org_metadata} \\
        --reference_accession ${params.fhir_reference_accession} \\
        --gff ${params.gff}

    cat <<-END_VERSIONS > versions.yml
    "fhir_converter":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    echo '{}' > ${meta.id}.fhir.json
    touch versions.yml
    """
}

workflow FHIR {
    take:
    per_sample      
    org_metadata_ch

    main:
    CREATE_FHIR(per_sample, org_metadata_ch)

    emit:
    fhir_output = CREATE_FHIR.out.fhir_output
    versions    = CREATE_FHIR.out.versions
}
