#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process UPLOAD_TO_FHIR {
    tag "${fhir_file.baseName}"
    publishDir "${params.results_dir}/fhir_upload", mode: 'copy'

    input:
    path(fhir_file)

    output:
    path "${fhir_file.baseName}.upload.json", emit: upload_result

    script:
    def token_file  = "${projectDir}/data/access_token.json"
    def api_key_arg = params.api_key ? "--api_key '${params.api_key}'" : ""
    """
    python3 ${projectDir}/scripts/upload_fhir.py \\
        --fhir_file ${fhir_file} \\
        --fhir_server_url '${params.fhir_server_url}' \\
        --token_file '${token_file}' \\
        ${api_key_arg} \\
        --output ${fhir_file.baseName}.upload.json
    """

    stub:
    "echo '{}' > ${fhir_file.baseName}.upload.json"
}

workflow UPLOAD_FHIR {
    take:
    validated_fhir_files

    main:
    if (!params.fhir_server_url) {
        log.info "UPLOAD_FHIR: no --fhir_server_url configured, skipping upload"
        results = Channel.empty()
    }
    else {
        UPLOAD_TO_FHIR(validated_fhir_files)
        results = UPLOAD_TO_FHIR.out.upload_result
    }

    emit:
    results = results
}
