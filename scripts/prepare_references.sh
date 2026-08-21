#!/usr/bin/env bash

set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${BASE_DIR}/data"
TOOLS_DIR="${BASE_DIR}/tools"

CONTIG="MN908947.3"

log() { printf '\n[prepare_references] %s\n' "$*"; }
have() { command -v "$1" >/dev/null 2>&1; }

mkdir -p "${DATA_DIR}" "${TOOLS_DIR}"

NFCORE_GENOME="https://github.com/nf-core/test-datasets/raw/viralrecon/genome/MN908947.3"

if [[ ! -f "${DATA_DIR}/${CONTIG}.fasta" ]]; then
    log "Downloading reference genome ${CONTIG}"
    curl -fsSL "${NFCORE_GENOME}/GCA_009858895.3_ASM985889v3_genomic.200409.fna.gz" \
        | gunzip -c > "${DATA_DIR}/${CONTIG}.fasta.tmp"
    awk -v c="${CONTIG}" 'NR==1 {print ">" c; next} {print}' \
        "${DATA_DIR}/${CONTIG}.fasta.tmp" > "${DATA_DIR}/${CONTIG}.fasta"
    rm -f "${DATA_DIR}/${CONTIG}.fasta.tmp"
else
    log "Reference genome already present, skipping"
fi

if [[ ! -f "${DATA_DIR}/${CONTIG}.gff" ]]; then
    log "Downloading reference GFF"
    curl -fsSL "${NFCORE_GENOME}/GCA_009858895.3_ASM985889v3_genomic.200409.gff.gz" \
        | gunzip -c > "${DATA_DIR}/${CONTIG}.gff.tmp"
    awk -v c="${CONTIG}" 'BEGIN{FS=OFS="\t"} /^#/ {print; next} {$1=c; print}' \
        "${DATA_DIR}/${CONTIG}.gff.tmp" > "${DATA_DIR}/${CONTIG}.gff"
    rm -f "${DATA_DIR}/${CONTIG}.gff.tmp"
else
    log "Reference GFF already present, skipping"
fi

log "Indexing reference"
samtools faidx "${DATA_DIR}/${CONTIG}.fasta"

GENOME_LEN=$(cut -f2 "${DATA_DIR}/${CONTIG}.fasta.fai")
if [[ "${GENOME_LEN}" != "29903" ]]; then
    echo "ERROR: expected a 29903 bp genome, got ${GENOME_LEN} bp" >&2
    exit 1
fi
log "Reference OK: ${CONTIG}, ${GENOME_LEN} bp"

if [[ ! -f "${DATA_DIR}/bowtie2_index/${CONTIG}.1.bt2" ]]; then
    log "Building bowtie2 index"
    mkdir -p "${DATA_DIR}/bowtie2_index"
    bowtie2-build --threads 4 "${DATA_DIR}/${CONTIG}.fasta" \
        "${DATA_DIR}/bowtie2_index/${CONTIG}" >/dev/null
else
    log "Bowtie2 index already present, skipping"
fi

ARTIC_NCOV="https://github.com/artic-network/artic-ncov2019/raw/master/primer_schemes/nCoV-2019"
VARSKIP="https://github.com/nebiolabs/VarSkip/raw/main"

declare -A SCHEME_URLS=(
    [artic_v3]="${ARTIC_NCOV}/V3/nCoV-2019.primer.bed"
    [artic_v4_1]="${ARTIC_NCOV}/V4.1/SARS-CoV-2.primer.bed"
    [artic_v5_3_2]="${ARTIC_NCOV}/V5.3.2/SARS-CoV-2.primer.bed"
    [neb_vsl1]="${VARSKIP}/neb_vsl1a.primer.bed"
    [neb_vss2]="${VARSKIP}/neb_vss2b.primer.bed"
)

mkdir -p "${DATA_DIR}/primer_schemes"

for scheme in "${!SCHEME_URLS[@]}"; do
    scheme_dir="${DATA_DIR}/primer_schemes/${scheme}"
    mkdir -p "${scheme_dir}"

    if [[ ! -s "${scheme_dir}/primer.bed" ]]; then
        log "Downloading primer scheme ${scheme}"
        curl -fsSL "${SCHEME_URLS[$scheme]}" > "${scheme_dir}/primer.bed.tmp"
        tr -d '\r' < "${scheme_dir}/primer.bed.tmp" \
        | awk -v c="${CONTIG}" 'BEGIN{FS=OFS="\t"}
            /^#/ || /^track/ || /^browser/ {next}
            NF < 4 {next}
            {
                strand = (NF >= 6 && ($6 == "+" || $6 == "-")) ? $6 : \
                         ($4 ~ /RIGHT/ ? "-" : "+")
                pool = (NF >= 5 && $5 != "") ? $5 : "1"
                print c, $2, $3, $4, pool, strand
            }' > "${scheme_dir}/primer.bed"
        rm -f "${scheme_dir}/primer.bed.tmp"
    fi

    if [[ ! -s "${scheme_dir}/primers.fasta" ]]; then
        log "Deriving primer FASTA for ${scheme}"
        bedtools getfasta -fi "${DATA_DIR}/${CONTIG}.fasta" \
            -bed "${scheme_dir}/primer.bed" -name -s \
            > "${scheme_dir}/primers.fasta"
    fi

    n_primers=$(grep -vc '^#' "${scheme_dir}/primer.bed" || true)
    log "  ${scheme}: ${n_primers} primers"
done

NEXTCLADE_TAG="${NEXTCLADE_DATASET_TAG:-}"
if [[ ! -f "${DATA_DIR}/nextclade_db/pathogen.json" ]]; then
    log "Fetching Nextclade SARS-CoV-2 dataset"
    mkdir -p "${DATA_DIR}/nextclade_db"
    if [[ -n "${NEXTCLADE_TAG}" ]]; then
        nextclade dataset get --name nextstrain/sars-cov-2 \
            --tag "${NEXTCLADE_TAG}" --output-dir "${DATA_DIR}/nextclade_db"
    else
        nextclade dataset get --name nextstrain/sars-cov-2 \
            --output-dir "${DATA_DIR}/nextclade_db"
    fi
else
    log "Nextclade dataset already present, skipping"
fi

SNPEFF_DATA="${DATA_DIR}/snpeff_db"
SNPEFF_CONFIG="${DATA_DIR}/snpEff.config"

if [[ ! -f "${SNPEFF_DATA}/${CONTIG}/snpEffectPredictor.bin" ]]; then
    log "Building snpEff database for ${CONTIG}"
    mkdir -p "${SNPEFF_DATA}/${CONTIG}"
    cp "${DATA_DIR}/${CONTIG}.gff"   "${SNPEFF_DATA}/${CONTIG}/genes.gff"
    cp "${DATA_DIR}/${CONTIG}.fasta" "${SNPEFF_DATA}/${CONTIG}/sequences.fa"

    cat > "${SNPEFF_CONFIG}" <<EOF
# Minimal snpEff config for the SARS-CoV-2 pipeline.
data.dir = ${SNPEFF_DATA}
${CONTIG}.genome : SARS-CoV-2 (Wuhan-Hu-1)
${CONTIG}.chromosomes : ${CONTIG}
EOF

    snpEff build -gff3 -v -noCheckCds -noCheckProtein \
        -config "${SNPEFF_CONFIG}" -dataDir "${SNPEFF_DATA}" "${CONTIG}"
else
    log "snpEff database already present, skipping"
fi

if [[ ! -s "${TOOLS_DIR}/fhir-validator.jar" ]]; then
    log "Downloading HL7 FHIR validator"
    curl -fsSL -o "${TOOLS_DIR}/fhir-validator.jar" \
        "https://github.com/hapifhir/org.hl7.fhir.core/releases/latest/download/validator_cli.jar"
else
    log "FHIR validator already present, skipping"
fi

log "Done. Reference artefacts are under ${DATA_DIR}"
