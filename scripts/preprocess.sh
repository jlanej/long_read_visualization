#!/usr/bin/env bash
# preprocess.sh - Long-read pre-processing pipeline
#
# Aligns haplotype assemblies to a reference genome and reads to assemblies,
# then builds coordinate-mapping indices for linked IGV.js visualization.
#
# Required tools: minimap2, samtools, bgzip, tabix, python3
set -euo pipefail

# ── Defaults ────────────────────────────────────────────────────────────────
THREADS=4
READ_TYPE="ont"
CRAM_REF=""
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="${SCRIPT_DIR}/../src"

# ── Usage ───────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Pre-process long-read data for multi-panel IGV.js visualization.

Required:
  -s, --sample-dir DIR    Sample directory containing CRAM and assembly FASTAs
  -r, --reference  FILE   Target reference genome FASTA
  -o, --output-dir DIR    Output directory

Optional:
  -t, --threads    INT    Number of threads [${THREADS}]
  --read-type      STR    Read technology: ont | hifi [${READ_TYPE}]
  --cram-ref       FILE   Reference FASTA used to encode the CRAM
                          (required when different from --reference)

Expected sample directory layout:
  SAMPLE/
    SAMPLE.*.cram          (long-read CRAM)
    SAMPLE.*.cram.crai
    SAMPLE_hap1*.fa.gz     (haplotype 1 assembly)
    SAMPLE_hap2*.fa.gz     (haplotype 2 assembly)
EOF
    exit 1
}

# ── Argument parsing ────────────────────────────────────────────────────────
[[ $# -eq 0 ]] && usage

SAMPLE_DIR=""
REFERENCE=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -s|--sample-dir) SAMPLE_DIR="$2"; shift 2 ;;
        -r|--reference)  REFERENCE="$2";  shift 2 ;;
        -o|--output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        -t|--threads)    THREADS="$2";    shift 2 ;;
        --read-type)     READ_TYPE="$2";  shift 2 ;;
        --cram-ref)      CRAM_REF="$2";   shift 2 ;;
        -h|--help)       usage ;;
        *)               echo "Unknown option: $1" >&2; usage ;;
    esac
done

[[ -z "${SAMPLE_DIR}" ]] && { echo "ERROR: --sample-dir is required" >&2; usage; }
[[ -z "${REFERENCE}" ]]  && { echo "ERROR: --reference is required"  >&2; usage; }
[[ -z "${OUTPUT_DIR}" ]] && { echo "ERROR: --output-dir is required" >&2; usage; }

# ── Validate read type ──────────────────────────────────────────────────────
case "${READ_TYPE}" in
    ont)  MM2_READ_PRESET="map-ont"  ;;
    hifi) MM2_READ_PRESET="map-hifi" ;;
    *)    echo "ERROR: --read-type must be 'ont' or 'hifi'" >&2; exit 1 ;;
esac

# ── Discover input files ────────────────────────────────────────────────────
SAMPLE_DIR="$(cd "${SAMPLE_DIR}" && pwd)"
SAMPLE_NAME="$(basename "${SAMPLE_DIR}")"

echo "=== Long-read pre-processing pipeline ==="
echo "Sample:     ${SAMPLE_NAME}"
echo "Sample dir: ${SAMPLE_DIR}"
echo "Reference:  ${REFERENCE}"
echo "Output:     ${OUTPUT_DIR}"
echo "Threads:    ${THREADS}"
echo "Read type:  ${READ_TYPE} (minimap2 -x ${MM2_READ_PRESET})"
echo ""

find_single_match() {
    local dir="$1" pattern="$2" label="$3"
    local matches
    matches=$(find "${dir}" -maxdepth 1 -name "${pattern}" | head -2)
    local count
    count=$(echo "${matches}" | grep -c . || true)
    if [[ ${count} -eq 0 ]]; then
        echo "ERROR: no ${label} found matching ${pattern} in ${dir}" >&2
        exit 1
    fi
    echo "${matches}" | head -1
}

CRAM=$(find_single_match "${SAMPLE_DIR}" "*.cram" "CRAM file")
HAP1=$(find_single_match "${SAMPLE_DIR}" "*hap1*.fa*" "hap1 assembly")
HAP2=$(find_single_match "${SAMPLE_DIR}" "*hap2*.fa*" "hap2 assembly")

echo "CRAM: ${CRAM}"
echo "Hap1: ${HAP1}"
echo "Hap2: ${HAP2}"
echo ""

# ── Prepare output directory ────────────────────────────────────────────────
mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR="$(cd "${OUTPUT_DIR}" && pwd)"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ── Helper: index FASTA if needed ───────────────────────────────────────────
ensure_faidx() {
    local fa="$1"
    if [[ ! -f "${fa}.fai" ]]; then
        log "Indexing ${fa} ..."
        samtools faidx "${fa}"
    fi
}

# ── Step 1: Index references and assemblies ─────────────────────────────────
log "Step 1: Ensuring FASTA indices exist"
ensure_faidx "${REFERENCE}"
ensure_faidx "${HAP1}"
ensure_faidx "${HAP2}"

# ── Step 2: Align assemblies to reference (BAM) ────────────────────────────
align_asm_to_ref() {
    local asm="$1" label="$2"
    local bam="${OUTPUT_DIR}/${SAMPLE_NAME}_${label}_to_ref.bam"
    if [[ -f "${bam}" ]]; then
        log "  ${bam} exists, skipping"
    else
        log "  Aligning ${label} to reference → ${bam}"
        minimap2 -a --eqx -x asm5 -t "${THREADS}" "${REFERENCE}" "${asm}" \
            | samtools sort -@ "${THREADS}" -o "${bam}"
        samtools index -@ "${THREADS}" "${bam}"
    fi
    echo "${bam}"
}

log "Step 2: Aligning assemblies to reference"
HAP1_BAM=$(align_asm_to_ref "${HAP1}" "hap1")
HAP2_BAM=$(align_asm_to_ref "${HAP2}" "hap2")

# ── Step 3: Produce PAF for coordinate mapping ─────────────────────────────
align_asm_paf() {
    local asm="$1" label="$2"
    local paf="${OUTPUT_DIR}/${SAMPLE_NAME}_${label}_to_ref.paf"
    if [[ -f "${paf}" ]]; then
        log "  ${paf} exists, skipping"
    else
        log "  Generating PAF for ${label} → ${paf}"
        minimap2 --eqx -c -x asm5 -t "${THREADS}" "${REFERENCE}" "${asm}" > "${paf}"
    fi
    echo "${paf}"
}

log "Step 3: Generating PAF alignments for coordinate mapping"
HAP1_PAF=$(align_asm_paf "${HAP1}" "hap1")
HAP2_PAF=$(align_asm_paf "${HAP2}" "hap2")

# ── Step 4: Build coordinate mapping indices ────────────────────────────────
build_mapping() {
    local paf="$1" label="$2"
    local prefix="${OUTPUT_DIR}/${SAMPLE_NAME}_${label}_to_ref"
    local json="${prefix}.mapping.json.gz"
    if [[ -f "${json}" ]]; then
        log "  ${json} exists, skipping"
    else
        log "  Building mapping index for ${label}"
        python3 "${SRC_DIR}/coordinate_mapper.py" build -p "${paf}" -o "${prefix}"
        # Create tabix-indexed BED
        local bed="${prefix}.mapping.bed"
        if [[ -f "${bed}" ]]; then
            bgzip -f "${bed}"
            tabix -p bed "${bed}.gz"
        fi
    fi
}

log "Step 4: Building coordinate mapping indices"
build_mapping "${HAP1_PAF}" "hap1"
build_mapping "${HAP2_PAF}" "hap2"

# ── Step 5: Extract reads from CRAM ────────────────────────────────────────
FASTQ="${OUTPUT_DIR}/${SAMPLE_NAME}_reads.fastq.gz"

if [[ -f "${FASTQ}" ]]; then
    log "Step 5: ${FASTQ} exists, skipping read extraction"
else
    log "Step 5: Extracting reads from CRAM"
    CRAM_REF_OPT=""
    if [[ -n "${CRAM_REF}" ]]; then
        CRAM_REF_OPT="--reference ${CRAM_REF}"
    fi
    # shellcheck disable=SC2086
    if command -v pigz &>/dev/null; then
        samtools fastq -@ "${THREADS}" ${CRAM_REF_OPT} "${CRAM}" \
            | pigz -p "${THREADS}" > "${FASTQ}"
    else
        samtools fastq -@ "${THREADS}" ${CRAM_REF_OPT} "${CRAM}" \
            | gzip > "${FASTQ}"
    fi
fi

# ── Step 6: Align reads to hap1 and hap2 ───────────────────────────────────
align_reads_to_asm() {
    local asm="$1" label="$2"
    local bam="${OUTPUT_DIR}/${SAMPLE_NAME}_reads_to_${label}.bam"
    if [[ -f "${bam}" ]]; then
        log "  ${bam} exists, skipping"
    else
        log "  Aligning reads to ${label} → ${bam}"
        minimap2 -a -x "${MM2_READ_PRESET}" -t "${THREADS}" "${asm}" "${FASTQ}" \
            | samtools sort -@ "${THREADS}" -o "${bam}"
        samtools index -@ "${THREADS}" "${bam}"
    fi
}

log "Step 6: Aligning reads to hap1 and hap2"
align_reads_to_asm "${HAP1}" "hap1"
align_reads_to_asm "${HAP2}" "hap2"

# ── Done ────────────────────────────────────────────────────────────────────
log "Pipeline complete. Outputs in ${OUTPUT_DIR}/"
echo ""
echo "Output files:"
ls -lh "${OUTPUT_DIR}/"
echo ""
echo "Query coordinate mappings with:"
echo "  python3 ${SRC_DIR}/coordinate_mapper.py query \\"
echo "    -i ${OUTPUT_DIR}/${SAMPLE_NAME}_hap1_to_ref.mapping.json.gz \\"
echo "    -r chr1:1000000-2000000"
