#!/usr/bin/env bash
# preprocess.sh - Long-read pre-processing pipeline
#
# Aligns haplotype assemblies to a reference genome and reads to assemblies,
# then builds coordinate-mapping indices for linked IGV.js visualization.
#
# Required tools: minimap2, samtools, bgzip, tabix, python3, curl
set -euo pipefail

# ── Defaults ────────────────────────────────────────────────────────────────
THREADS=4
READ_TYPE="ont"
CRAM_REF=""
REMAP_CRAM_TO_REFERENCE=0
REF_CACHE_DIR="${PWD}/references"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="${SCRIPT_DIR}/../src"

# Known reference genomes (name → URL of bgzipped/gzipped FASTA)
declare -A GENOME_URLS
GENOME_URLS["hg38"]="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
GENOME_URLS["hg19"]="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
GENOME_URLS["chm13v2.0"]="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
GENOME_URLS["grch38"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"

# ── Usage ───────────────────────────────────────────────────────────────────
usage() {
    local known_genomes
    known_genomes="${!GENOME_URLS[*]}"
    cat <<EOF
Usage: $(basename "$0") [options]

Pre-process long-read data for multi-panel IGV.js visualization.

Required inputs:
  --hap1        FILE   Haplotype 1 assembly FASTA (.fa or .fa.gz)
  --hap2        FILE   Haplotype 2 assembly FASTA (.fa or .fa.gz)
  --cram        FILE   Long-read CRAM file
  -o, --output-dir DIR Output directory

Read type (controls minimap2 alignment preset):
  --ont                Oxford Nanopore reads  (minimap2 -x map-ont)  [default]
  --hifi               PacBio HiFi reads      (minimap2 -x map-hifi)

Reference (one of):
  -r, --reference FILE Target reference genome FASTA
  --genome        STR  Download a known reference if not already cached.
                       Supported: ${known_genomes}
  --ref-dir       DIR  Cache directory for downloaded references
                       [${REF_CACHE_DIR}]

Optional:
  -t, --threads   INT  Number of threads [${THREADS}]
  --cram-ref      FILE Reference FASTA used to encode the CRAM
                        (required when different from --reference)
  --remap-cram-to-reference
                       Remap input CRAM to --reference/--genome, then tag HP
                       in-place on the remapped CRAM (no extra *_reads.hp.cram)
  -h, --help           Show this help message
EOF
    exit 1
}

# ── Argument parsing ────────────────────────────────────────────────────────
[[ $# -eq 0 ]] && usage

HAP1=""
HAP2=""
CRAM=""
REFERENCE=""
GENOME=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --hap1)          HAP1="$2";          shift 2 ;;
        --hap2)          HAP2="$2";          shift 2 ;;
        --cram)          CRAM="$2";          shift 2 ;;
        -r|--reference)  REFERENCE="$2";     shift 2 ;;
        --genome)        GENOME="$2";        shift 2 ;;
        --ref-dir)       REF_CACHE_DIR="$2"; shift 2 ;;
        -o|--output-dir) OUTPUT_DIR="$2";    shift 2 ;;
        -t|--threads)    THREADS="$2";       shift 2 ;;
        --ont)           READ_TYPE="ont";    shift   ;;
        --hifi)          READ_TYPE="hifi";   shift   ;;
        --read-type)     READ_TYPE="$2";     shift 2 ;;
        --cram-ref)      CRAM_REF="$2";      shift 2 ;;
        --remap-cram-to-reference) REMAP_CRAM_TO_REFERENCE=1; shift ;;
        -h|--help)       usage ;;
        *)               echo "Unknown option: $1" >&2; usage ;;
    esac
done

[[ -z "${HAP1}" ]]       && { echo "ERROR: --hap1 is required" >&2; usage; }
[[ -z "${HAP2}" ]]       && { echo "ERROR: --hap2 is required" >&2; usage; }
[[ -z "${CRAM}" ]]       && { echo "ERROR: --cram is required" >&2; usage; }
[[ -z "${OUTPUT_DIR}" ]] && { echo "ERROR: --output-dir is required" >&2; usage; }

# Validate that input files exist
[[ -f "${HAP1}" ]] || { echo "ERROR: hap1 file not found: ${HAP1}" >&2; exit 1; }
[[ -f "${HAP2}" ]] || { echo "ERROR: hap2 file not found: ${HAP2}" >&2; exit 1; }
[[ -f "${CRAM}" ]] || { echo "ERROR: CRAM file not found: ${CRAM}" >&2; exit 1; }
[[ -n "${REFERENCE}" && ! -f "${REFERENCE}" ]] && \
    { echo "ERROR: reference file not found: ${REFERENCE}" >&2; exit 1; }

if [[ -z "${REFERENCE}" && -z "${GENOME}" ]]; then
    echo "ERROR: one of --reference or --genome is required" >&2
    usage
fi

# ── Validate read type ──────────────────────────────────────────────────────
case "${READ_TYPE}" in
    ont)  MM2_READ_PRESET="map-ont"  ;;
    hifi) MM2_READ_PRESET="map-hifi" ;;
    *)    echo "ERROR: --read-type must be 'ont' or 'hifi'" >&2; exit 1 ;;
esac

# ── Helpers ─────────────────────────────────────────────────────────────────
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }

# Print and execute a command, so the exact invocation is always visible.
run() {
    log "CMD: $*"
    "$@"
}

ensure_faidx() {
    local fa="$1"
    if [[ ! -f "${fa}.fai" ]]; then
        log "Indexing ${fa} ..."
        run samtools faidx "${fa}"
    fi
}

# ── Genome download ──────────────────────────────────────────────────────────
download_reference() {
    local genome="$1"
    local cache_dir="$2"

    if [[ -z "${GENOME_URLS[${genome}]+x}" ]]; then
        echo "ERROR: Unknown genome '${genome}'." >&2
        echo "       Supported genomes: ${!GENOME_URLS[*]}" >&2
        exit 1
    fi

    local url="${GENOME_URLS[${genome}]}"
    local filename
    filename="$(basename "${url}")"
    local ref_path="${cache_dir}/${filename}"

    mkdir -p "${cache_dir}"

    if [[ -f "${ref_path}" ]]; then
        log "Cached reference found: ${ref_path}"
    else
        log "Downloading ${genome} reference from ${url} ..."
        run curl -fL --progress-bar -o "${ref_path}" "${url}"
        log "Download complete: ${ref_path}"
    fi

    echo "${ref_path}"
}

# ── Resolve reference (download if needed) ──────────────────────────────────
if [[ -n "${GENOME}" && -z "${REFERENCE}" ]]; then
    REFERENCE=$(download_reference "${GENOME}" "${REF_CACHE_DIR}")
fi

# ── Resolve input paths to absolute paths ───────────────────────────────────
# (done after file-existence checks so errors are clear)
HAP1="$(cd "$(dirname "${HAP1}")" && pwd)/$(basename "${HAP1}")"
HAP2="$(cd "$(dirname "${HAP2}")" && pwd)/$(basename "${HAP2}")"
CRAM="$(cd "$(dirname "${CRAM}")" && pwd)/$(basename "${CRAM}")"
REFERENCE="$(cd "$(dirname "${REFERENCE}")" && pwd)/$(basename "${REFERENCE}")"

# Derive sample name from the CRAM filename (everything before the first dot).
# Example: NA21110.t2t.cram → NA21110 ; sample.hifi.cram → sample
SAMPLE_NAME="$(basename "${CRAM}" | cut -d. -f1)"

echo "=== Long-read pre-processing pipeline ==="
echo "Sample:     ${SAMPLE_NAME}"
echo "Hap1:       ${HAP1}"
echo "Hap2:       ${HAP2}"
echo "CRAM:       ${CRAM}"
echo "Reference:  ${REFERENCE}"
echo "Output:     ${OUTPUT_DIR}"
echo "Threads:    ${THREADS}"
echo "Read type:  ${READ_TYPE} (minimap2 -x ${MM2_READ_PRESET})"
echo ""

# ── Prepare output directory ────────────────────────────────────────────────
mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR="$(cd "${OUTPUT_DIR}" && pwd)"

# ── Step 1: Index references and assemblies ─────────────────────────────────
log "Step 1: Ensuring FASTA indices exist"
ensure_faidx "${REFERENCE}"
ensure_faidx "${HAP1}"
ensure_faidx "${HAP2}"

# ── Step 2: Align assemblies to reference (BAM) ────────────────────────────
align_asm_to_ref() {
    local asm="$1" label="$2"
    local bam="${OUTPUT_DIR}/${SAMPLE_NAME}_${label}_to_ref.bam"
    if [[ -f "${bam}" && -f "${bam}.bai" ]]; then
        log "  ${bam} + index exist, skipping"
    else
        log "  Aligning ${label} to reference → ${bam}"
        log "CMD: minimap2 -a --eqx -x asm5 -t ${THREADS} ${REFERENCE} ${asm} | samtools sort -@ ${THREADS} -o ${bam}"
        minimap2 -a --eqx -x asm5 -t "${THREADS}" "${REFERENCE}" "${asm}" \
            | samtools sort -@ "${THREADS}" -o "${bam}"
        run samtools index -@ "${THREADS}" "${bam}"
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
        local paf_tmp="${paf}.tmp"
        run minimap2 --eqx -c -x asm5 -t "${THREADS}" "${REFERENCE}" "${asm}" > "${paf_tmp}"
        mv "${paf_tmp}" "${paf}"
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
    local bed_gz="${prefix}.mapping.bed.gz"
    local bed_tbi="${prefix}.mapping.bed.gz.tbi"
    if [[ -f "${json}" && -f "${bed_gz}" && -f "${bed_tbi}" ]]; then
        log "  ${label} mapping index exists (json + bed.gz + tbi), skipping"
    else
        log "  Building mapping index for ${label}"
        run python3 "${SRC_DIR}/coordinate_mapper.py" build -p "${paf}" -o "${prefix}"
        # Create tabix-indexed BED
        local bed="${prefix}.mapping.bed"
        if [[ -f "${bed}" ]]; then
            run bgzip -f "${bed}"
            run tabix -p bed "${bed}.gz"
        fi
    fi
}

log "Step 4: Building coordinate mapping indices"
build_mapping "${HAP1_PAF}" "hap1"
build_mapping "${HAP2_PAF}" "hap2"

# ── Resolve CRAM reference for samtools ─────────────────────────────────────
# Prefer an explicit --cram-ref; fall back to the primary --reference / --genome
# so that samtools can always decode the CRAM file.
CRAM_REF_OPT=""
if [[ -n "${CRAM_REF}" ]]; then
    CRAM_REF_OPT="--reference ${CRAM_REF}"
elif [[ -n "${REFERENCE}" ]]; then
    CRAM_REF_OPT="--reference ${REFERENCE}"
fi

# ── Step 5–6: Align reads to hap1 and hap2 (streaming from CRAM) ──────────
#
# Instead of writing an intermediate FASTQ file to disk (which can exceed
# 100–200 GB for whole-genome long-read datasets), we stream reads directly
# from the CRAM into minimap2 using process substitution.  This trades a
# second CRAM decode for vastly reduced disk I/O and storage.
align_reads_to_asm() {
    local asm="$1" label="$2"
    local bam="${OUTPUT_DIR}/${SAMPLE_NAME}_reads_to_${label}.bam"
    if [[ -f "${bam}" && -f "${bam}.bai" ]]; then
        log "  ${bam} + index exist, skipping"
    else
        log "  Aligning reads to ${label} → ${bam}"
        # shellcheck disable=SC2086
        log "CMD: minimap2 -a -x ${MM2_READ_PRESET} -t ${THREADS} ${asm} <(samtools fastq -@ ${THREADS} ${CRAM_REF_OPT} ${CRAM}) | samtools sort -@ ${THREADS} -o ${bam}"
        # shellcheck disable=SC2086
        minimap2 -a -x "${MM2_READ_PRESET}" -t "${THREADS}" "${asm}" \
            <(samtools fastq -@ "${THREADS}" ${CRAM_REF_OPT} "${CRAM}") \
            | samtools sort -@ "${THREADS}" -o "${bam}"
        run samtools index -@ "${THREADS}" "${bam}"
    fi
}

log "Step 5: Aligning reads to hap1 and hap2"
align_reads_to_asm "${HAP1}" "hap1"
align_reads_to_asm "${HAP2}" "hap2"

# ── Step 5b: Assign haplotype tags (HP) via competitive alignment scoring ──
#
# 1. Compare the minimap2 Alignment Score (AS:i:) for every read between the
#    hap1 and hap2 BAMs → compute HP:i:1/2/0 assignments.
# 2. Tag the pipeline-generated reads_to_hap{1,2}.bam files in-place (used
#    by the haplotype-space assembly panels in IGV).
# 3a. (Default) Create a new HP-tagged CRAM (${SAMPLE_NAME}_reads.hp.cram)
#     derived from the original --cram input (which is never modified).
# 3b. (Optional --remap-cram-to-reference) Remap input CRAM to the selected
#     reference and tag HP in-place on the remapped CRAM.
#
# HP:i:1 = hap1 wins  |  HP:i:2 = hap2 wins  |  HP:i:0 = ambiguous
# IGV natively groups, sorts, and colours reads by the HP tag.
# Idempotent: re-running produces the same assignments and output files.
HAP1_BAM="${OUTPUT_DIR}/${SAMPLE_NAME}_reads_to_hap1.bam"
HAP2_BAM="${OUTPUT_DIR}/${SAMPLE_NAME}_reads_to_hap2.bam"
READS_HP_CRAM="${OUTPUT_DIR}/${SAMPLE_NAME}_reads.hp.cram"

# Choose the best available reference for CRAM encoding
HP_REF="${CRAM_REF:-${REFERENCE}}"

log "Step 5b: Assigning haplotype tags (HP) based on competitive alignment scores"
if [[ "${REMAP_CRAM_TO_REFERENCE}" -eq 1 ]]; then
    if [[ -f "${READS_HP_CRAM}" && -f "${READS_HP_CRAM}.crai" ]]; then
        log "  Remapped HP-tagged CRAM exists, skipping"
    else
        log "  Remapping input CRAM to selected reference: ${READS_HP_CRAM}"
        run samtools view -@ "${THREADS}" -T "${REFERENCE}" -C -o "${READS_HP_CRAM}" "${CRAM}"

        # shellcheck disable=SC2086
        log "CMD: python3 ${SRC_DIR}/assign_haplotypes.py --hap1-bam ${HAP1_BAM} --hap2-bam ${HAP2_BAM} --reads-in ${READS_HP_CRAM} --reads-in-place --reference ${REFERENCE}"
        run python3 "${SRC_DIR}/assign_haplotypes.py" \
            --hap1-bam      "${HAP1_BAM}" \
            --hap2-bam      "${HAP2_BAM}" \
            --reads-in      "${READS_HP_CRAM}" \
            --reads-in-place \
            --reference     "${REFERENCE}"
    fi
elif [[ -f "${READS_HP_CRAM}" && -f "${READS_HP_CRAM}.crai" ]]; then
    log "  HP-tagged CRAM exists, skipping"
else
    # shellcheck disable=SC2086
    log "CMD: python3 ${SRC_DIR}/assign_haplotypes.py --hap1-bam ${HAP1_BAM} --hap2-bam ${HAP2_BAM} --reads-in ${CRAM} --reads-out ${READS_HP_CRAM} --reference ${HP_REF}"
    # shellcheck disable=SC2086
    run python3 "${SRC_DIR}/assign_haplotypes.py" \
        --hap1-bam  "${HAP1_BAM}" \
        --hap2-bam  "${HAP2_BAM}" \
        --reads-in  "${CRAM}" \
        --reads-out "${READS_HP_CRAM}" \
        --reference "${HP_REF}"
fi

# ── Step 6: Align reference to assemblies (ref-on-asm cross-check BAMs) ────
#
# These BAMs are coordinate-sorted in *assembly* space so the reference
# sequence can be loaded as a second track in any assembly-space genome
# browser (e.g. IGV.js).  They serve as a reciprocal cross-check for the
# assembly-to-reference alignments produced in Step 2: every alignment
# block that appears in the hap→ref BAM should have a complementary block
# in the ref→hap BAM, making it straightforward to spot spurious mappings
# or missed alignments.
#
# Preset: asm5  (≥99% identity, same as Step 2 — query and target are
#               simply swapped so the output is sorted on assembly coords)
# Flag:   --eqx  extended CIGAR (=/X instead of M) for precise
#               match/mismatch visibility in the browser
align_ref_to_asm() {
    local asm="$1" label="$2"
    local bam="${OUTPUT_DIR}/${SAMPLE_NAME}_ref_to_${label}.bam"
    if [[ -f "${bam}" && -f "${bam}.bai" ]]; then
        log "  ${bam} + index exist, skipping"
    else
        log "  Aligning reference to ${label} → ${bam}"
        log "CMD: minimap2 -a --eqx -x asm5 -t ${THREADS} ${asm} ${REFERENCE} | samtools sort -@ ${THREADS} -o ${bam}"
        minimap2 -a --eqx -x asm5 -t "${THREADS}" "${asm}" "${REFERENCE}" \
            | samtools sort -@ "${THREADS}" -o "${bam}"
        run samtools index -@ "${THREADS}" "${bam}"
    fi
}

log "Step 6: Aligning reference to hap1 and hap2 (assembly-space cross-check BAMs)"
align_ref_to_asm "${HAP1}" "hap1"
align_ref_to_asm "${HAP2}" "hap2"

# ── Step 7: Align hap2 to hap1 and hap1 to hap2 (cross-haplotype BAMs) ────
#
# These BAMs show one haplotype aligned to the other in the target's
# coordinate space.  They let the user see directly how the two
# haplotypes relate without going through the reference.
#
# Preset: asm5  (≥99% identity — both assemblies come from the same
#               individual so most loci are nearly identical)
# Flag:   --eqx  extended CIGAR (=/X)
align_hap_to_hap() {
    local query_asm="$1" query_label="$2"
    local target_asm="$3" target_label="$4"
    local bam="${OUTPUT_DIR}/${SAMPLE_NAME}_${query_label}_to_${target_label}.bam"
    if [[ -f "${bam}" && -f "${bam}.bai" ]]; then
        log "  ${bam} + index exist, skipping"
    else
        log "  Aligning ${query_label} to ${target_label} → ${bam}"
        log "CMD: minimap2 -a --eqx -x asm5 -t ${THREADS} ${target_asm} ${query_asm} | samtools sort -@ ${THREADS} -o ${bam}"
        minimap2 -a --eqx -x asm5 -t "${THREADS}" "${target_asm}" "${query_asm}" \
            | samtools sort -@ "${THREADS}" -o "${bam}"
        run samtools index -@ "${THREADS}" "${bam}"
    fi
}

log "Step 7: Aligning hap2 to hap1 and hap1 to hap2 (cross-haplotype BAMs)"
align_hap_to_hap "${HAP2}" "hap2" "${HAP1}" "hap1"
align_hap_to_hap "${HAP1}" "hap1" "${HAP2}" "hap2"

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
