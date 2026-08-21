#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process VALIDATE_FHIR {
    tag "${fhir_file.baseName}"
    publishDir "${params.results_dir}/fhir_validated", mode: 'copy'

    input:
    path(fhir_file)
    path(validator)

    output:
    path "${fhir_file.baseName}.validation.txt", emit: validation_report
    path "${fhir_file.baseName}.validation.summary.txt", emit: validation_summary
    path "${fhir_file}", emit: validated_fhir

    script:
    """
    java -jar ${validator} \\
        ${fhir_file} \\
        -version 4.0.1 \\
        -ig hl7.fhir.uv.genomics-reporting#current \\
        > ${fhir_file.baseName}.validation.txt 2>&1 || true

    sed 's/\\x1b\\[[0-9;]*m//g' ${fhir_file.baseName}.validation.txt \\
        > plain.txt

    n_error=\$(grep -cE '^[[:space:]]*Error' plain.txt || true)
    n_warn=\$(grep -cE '^[[:space:]]*Warning' plain.txt || true)

    {
        echo "FHIR validation summary: ${fhir_file.baseName}"
        echo "  profile   : hl7.fhir.uv.genomics-reporting#current (FHIR R4 4.0.1)"
        echo "  errors    : \${n_error}"
        echo "  warnings  : \${n_warn}"
        echo ""
        echo "Distinct error classes:"
        grep -E '^[[:space:]]*Error' plain.txt \\
            | sed -E 's/\\(line [0-9]+, col[0-9]+\\)//; s/Bundle\\.entry\\[[0-9]+\\]/Bundle.entry[N]/g; s/component\\[[0-9]+\\]/component[N]/g; s/codes = [^)]*/codes = .../' \\
            | sort | uniq -c | sort -rn || echo "  (none)"
    } > ${fhir_file.baseName}.validation.summary.txt

    echo "VALIDATE_FHIR ${fhir_file.baseName}: \${n_error} error(s), \${n_warn} warning(s)"
    """

    stub:
    """
    touch ${fhir_file.baseName}.validation.txt
    touch ${fhir_file.baseName}.validation.summary.txt
    """
}

workflow VALIDATE {
    take:
    fhir_json_files

    main:
    validator_path = file(params.fhir_validator)
    VALIDATE_FHIR(fhir_json_files.flatten(), validator_path)

    emit:
    validation_report  = VALIDATE_FHIR.out.validation_report
    validation_summary = VALIDATE_FHIR.out.validation_summary
    validated_fhir     = VALIDATE_FHIR.out.validated_fhir
}
