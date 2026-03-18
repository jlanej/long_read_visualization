#!/usr/bin/env bash
# generate_toy_dataset.sh - Generate a minimal toy dataset for integration testing.
#
# This script wraps src/generate_toy_dataset.py and is meant to be run AFTER
# a full preprocess.sh pipeline run.  It selects large deletions from the
# NA21110 SV callset bundled in resources/, uses the coordinate-mapping
# indices produced by the pipeline to find the corresponding assembly regions,
# then extracts small subsets of the reference, assemblies, and reads.
#
# Required tools: samtools, bgzip, python3
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="${SCRIPT_DIR}/../src"
RESOURCES_DIR="${SCRIPT_DIR}/../resources"

# ── Defaults ────────────────────────────────────────────────────────────────
NUM_VARIANTS=10
MIN_SIZE=1000
PADDING=50000

# ── Usage ───────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Generate a minimal toy dataset for integration testing by extracting regions
around large deletions from NA21110.

This script must be run AFTER a full preprocess.sh pipeline run because it
requires the coordinate-mapping JSON indices.

Required:
  --pipeline-output DIR  Output directory from a previous preprocess.sh run
  --hap1           FILE  Haplotype 1 assembly FASTA (.fa.gz)
  --hap2           FILE  Haplotype 2 assembly FASTA (.fa.gz)
  --reference      FILE  Reference genome FASTA
  --cram           FILE  Long-read CRAM file
  -o, --output-dir DIR   Output directory for the toy dataset

Optional:
  --vcf            FILE  SV callset VCF [bundled NA21110 VCF]
  --cram-ref       FILE  Reference FASTA used to encode the CRAM
  -n, --num-variants INT Number of variants to select [${NUM_VARIANTS}]
  --min-size       INT   Minimum deletion size in bp [${MIN_SIZE}]
  --padding        INT   Padding around each deletion [${PADDING}]
  --sample-name    STR   Sample name (auto-detected from CRAM filename)
  -h, --help             Show this help message
EOF
    exit 1
}

# ── Argument parsing ────────────────────────────────────────────────────────
[[ $# -eq 0 ]] && usage

PIPELINE_OUTPUT=""
HAP1=""
HAP2=""
CRAM=""
REFERENCE=""
VCF=""
CRAM_REF=""
OUTPUT_DIR=""
SAMPLE_NAME=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pipeline-output) PIPELINE_OUTPUT="$2"; shift 2 ;;
        --hap1)            HAP1="$2";            shift 2 ;;
        --hap2)            HAP2="$2";            shift 2 ;;
        --cram)            CRAM="$2";            shift 2 ;;
        --reference)       REFERENCE="$2";       shift 2 ;;
        --vcf)             VCF="$2";             shift 2 ;;
        --cram-ref)        CRAM_REF="$2";        shift 2 ;;
        -o|--output-dir)   OUTPUT_DIR="$2";      shift 2 ;;
        -n|--num-variants) NUM_VARIANTS="$2";    shift 2 ;;
        --min-size)        MIN_SIZE="$2";        shift 2 ;;
        --padding)         PADDING="$2";         shift 2 ;;
        --sample-name)     SAMPLE_NAME="$2";     shift 2 ;;
        -h|--help)         usage ;;
        *)                 echo "Unknown option: $1" >&2; usage ;;
    esac
done

# ── Validate required arguments ─────────────────────────────────────────────
[[ -z "${PIPELINE_OUTPUT}" ]] && { echo "ERROR: --pipeline-output is required" >&2; usage; }
[[ -z "${HAP1}" ]]            && { echo "ERROR: --hap1 is required" >&2; usage; }
[[ -z "${HAP2}" ]]            && { echo "ERROR: --hap2 is required" >&2; usage; }
[[ -z "${CRAM}" ]]            && { echo "ERROR: --cram is required" >&2; usage; }
[[ -z "${REFERENCE}" ]]       && { echo "ERROR: --reference is required" >&2; usage; }
[[ -z "${OUTPUT_DIR}" ]]      && { echo "ERROR: --output-dir is required" >&2; usage; }

# ── Derive sample name if not provided ──────────────────────────────────────
if [[ -z "${SAMPLE_NAME}" ]]; then
    SAMPLE_NAME="$(basename "${CRAM}" | cut -d. -f1)"
fi

# ── Locate mapping indices ──────────────────────────────────────────────────
HAP1_INDEX="${PIPELINE_OUTPUT}/${SAMPLE_NAME}_hap1_to_ref.mapping.json.gz"
HAP2_INDEX="${PIPELINE_OUTPUT}/${SAMPLE_NAME}_hap2_to_ref.mapping.json.gz"

[[ -f "${HAP1_INDEX}" ]] || {
    echo "ERROR: hap1 mapping index not found: ${HAP1_INDEX}" >&2
    echo "       Run preprocess.sh first, then re-run this script." >&2
    exit 1
}
[[ -f "${HAP2_INDEX}" ]] || {
    echo "ERROR: hap2 mapping index not found: ${HAP2_INDEX}" >&2
    echo "       Run preprocess.sh first, then re-run this script." >&2
    exit 1
}

# ── Default VCF ─────────────────────────────────────────────────────────────
if [[ -z "${VCF}" ]]; then
    VCF="${RESOURCES_DIR}/NA21110.shapeit5-phased-callset_final-vcf.phased.vcf.gz"
fi
[[ -f "${VCF}" ]] || { echo "ERROR: VCF not found: ${VCF}" >&2; exit 1; }

# ── Build Python command ────────────────────────────────────────────────────
CRAM_REF_OPT=""
if [[ -n "${CRAM_REF}" ]]; then
    CRAM_REF_OPT="--cram-ref ${CRAM_REF}"
fi

echo "=== Generating toy dataset ==="
echo "Sample:          ${SAMPLE_NAME}"
echo "Pipeline output: ${PIPELINE_OUTPUT}"
echo "VCF:             ${VCF}"
echo "Hap1 index:      ${HAP1_INDEX}"
echo "Hap2 index:      ${HAP2_INDEX}"
echo "Num variants:    ${NUM_VARIANTS}"
echo "Output:          ${OUTPUT_DIR}"
echo ""

# shellcheck disable=SC2086
python3 "${SRC_DIR}/generate_toy_dataset.py" \
    --vcf "${VCF}" \
    --hap1-index "${HAP1_INDEX}" \
    --hap2-index "${HAP2_INDEX}" \
    --hap1 "${HAP1}" \
    --hap2 "${HAP2}" \
    --reference "${REFERENCE}" \
    --cram "${CRAM}" \
    --output-dir "${OUTPUT_DIR}" \
    --num-variants "${NUM_VARIANTS}" \
    --min-size "${MIN_SIZE}" \
    --padding "${PADDING}" \
    ${CRAM_REF_OPT}

echo ""
echo "Toy dataset ready in ${OUTPUT_DIR}/"
echo ""
echo "Run the pipeline on the toy dataset:"
echo "  apptainer run \\"
echo "      --bind \"\${PWD}:/work\" \\"
echo "      docker://ghcr.io/jlanej/long_read_visualization:main \\"
echo "      --hap1       /work/${OUTPUT_DIR}/toy_hap1.fa.gz \\"
echo "      --hap2       /work/${OUTPUT_DIR}/toy_hap2.fa.gz \\"
echo "      --reference  /work/${OUTPUT_DIR}/toy_reference.fa.gz \\"
echo "      --cram       /work/${OUTPUT_DIR}/toy_reads.bam \\"
echo "      --output-dir /work/${OUTPUT_DIR}/pipeline_output \\"
echo "      --ont"
