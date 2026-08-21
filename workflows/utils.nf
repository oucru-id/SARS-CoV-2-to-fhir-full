nextflow.enable.dsl = 2

process VERSIONS {
    publishDir "${params.results_dir}", mode: 'copy'

    output:
    path "software_versions.yml"

    script:
    """
    echo "pipeline:" > software_versions.yml
    echo "  name: sars_cov2_mutation_analysis" >> software_versions.yml
    echo "  version: ${params.version}" >> software_versions.yml
    echo "  nextflow: $nextflow.version" >> software_versions.yml

    echo "references:" >> software_versions.yml
    echo "  genome: \$(basename ${params.reference})" >> software_versions.yml
    echo "  contig: ${params.contig}" >> software_versions.yml
    echo "  gff: \$(basename ${params.gff})" >> software_versions.yml
    echo "  fhir_reference_accession: ${params.fhir_reference_accession}" >> software_versions.yml
    echo "  primer_scheme_default: ${params.primer_scheme}" >> software_versions.yml

    echo "processing_settings:" >> software_versions.yml
    echo "  variant_caller: ${params.variant_caller}" >> software_versions.yml
    echo "  consensus_caller: ${params.consensus_caller}" >> software_versions.yml
    echo "  min_coverage: ${params.min_coverage}" >> software_versions.yml
    echo "  min_allele_freq: ${params.min_allele_freq}" >> software_versions.yml
    echo "  min_base_quality: ${params.min_base_quality}" >> software_versions.yml
    echo "  host_removal: \$([ "${params.skip_host_removal}" = "true" ] && echo disabled || echo "hostile/${params.hostile_index}")" >> software_versions.yml
    echo "  markduplicates: \$([ "${params.skip_markduplicates}" = "true" ] && echo disabled || echo enabled)" >> software_versions.yml
    echo "  nanopore_min_q: ${params.nanopore_min_q}" >> software_versions.yml
    echo "  medaka_model_default: ${params.medaka_model}" >> software_versions.yml

    export BASE_DIR="${baseDir}"
    python3 $baseDir/scripts/get_versions.py >> software_versions.yml
    """

    stub:
    "touch software_versions.yml"
}
