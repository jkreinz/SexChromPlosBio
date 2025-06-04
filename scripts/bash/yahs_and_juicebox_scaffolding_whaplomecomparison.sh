#!/bin/bash

# Hi-C Assembly Pipeline using YaHS and Juicer
# Author: Julia Kreiner
# Description: Pipeline for scaffolding genome assemblies using Hi-C data

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

# =============================================================================
# CONFIGURATION - Update these paths as needed
# =============================================================================

# Base directories
BASE_DIR="/ohta2/julia.kreiner/hifi_hic_malegenomes"
ASSEMBLY_DIR="${BASE_DIR}/Nune_assembly"
YAHS_DIR="${BASE_DIR}/yahs"
SOFTWARE_DIR="/ohta2/julia.kreiner/software"

# Sample prefix
PREFIX="Ama_N5"

# Tool paths
YAHS="${YAHS_DIR}/yahs"
JUICER_YAHS="${YAHS_DIR}/juicer"
JUICER_TOOLS="${SOFTWARE_DIR}/juicer_tools_1.22.01.jar"

# Input files
HAP1_FASTA="${ASSEMBLY_DIR}/Ama_N5_HiC_s40_hap1.fasta"
HAP2_FASTA="${ASSEMBLY_DIR}/Ama_N5_HiC_s40_hap2.fasta"
HAP1_BAM="${ASSEMBLY_DIR}/yahs/hap1/filt_bams/Ama_N5_HiC_DSN_filtsortrgdedup.bam"
HAP2_BAM="${ASSEMBLY_DIR}/yahs/hap2/filt_bams/Ama_N5_HiC_DSN_filtsortrgdedup.bam"

# YaHS output files
BIN_FILE="${ASSEMBLY_DIR}/yahs/hap2/yahs.out.bin"
SCAFFOLDS_AGP="${ASSEMBLY_DIR}/yahs/hap2/yahs.out_scaffolds_final.agp"
CONTIG_FAI="${BASE_DIR}/193_assembly/Ama_193_HiC_s40_hap2.fasta.fai"

# Final assembly files
FINAL_HAP1="${BASE_DIR}/final_assemblies/Atub_193_hap1.fasta"
FINAL_HAP2="${BASE_DIR}/final_assemblies/Atub_193_hap2.fasta"
ORIGINAL_HAP2_CONTIGS="${BASE_DIR}/hifiasm_run2_hets45/hifi_output/19_3_ATUB_s45.asm.hic.hap2.p_ctg.fasta"

# =============================================================================
# STEP 1: RUN YAHS SCAFFOLDING
# =============================================================================

echo "Step 1: Running YaHS scaffolding..."

# Scaffold haplotype 1
echo "  Scaffolding haplotype 1..."
"${YAHS}" "${HAP1_FASTA}" "${HAP1_BAM}"

# Scaffold haplotype 2
echo "  Scaffolding haplotype 2..."
"${YAHS}" "${HAP2_FASTA}" "${HAP2_BAM}"

# =============================================================================
# STEP 2: CONVERT TO JUICER FORMAT
# =============================================================================

echo "Step 2: Converting to Juicer format..."

# Convert HiC alignment to Juicer format
echo "  Converting alignments to Juicer format..."
("${JUICER_YAHS}" pre "${BIN_FILE}" "${SCAFFOLDS_AGP}" "${CONTIG_FAI}" | \
 sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | \
 awk 'NF' > alignments_sorted.txt.part) && \
(mv alignments_sorted.txt.part alignments_sorted.txt)

# Create scaffold size file
echo "  Creating scaffold size file..."
awk '{print $1 "\t" $2}' "${CONTIG_FAI}" > "${PREFIX}.scafsize"

# =============================================================================
# STEP 3: CREATE HI-C CONTACT MATRIX
# =============================================================================

echo "Step 3: Creating Hi-C contact matrix..."

# Create standard Hi-C contact matrix
echo "  Creating standard contact matrix..."
(java -jar -Xmx32G "${JUICER_TOOLS}" pre -j 10 \
 alignments_sorted.txt out.hic.part "${PREFIX}.scafsize") && \
(mv out.hic.part out.hic)

# Create Juicebox Assembly Tools (JBAT) compatible format
echo "  Creating JBAT-compatible format..."
"${JUICER_YAHS}" pre -a -o out_JBAT "${BIN_FILE}" "${SCAFFOLDS_AGP}" "${CONTIG_FAI}" > out_JBAT.log 2>&1

(java -jar -Xmx64G "${JUICER_TOOLS}" pre -j 10 \
 out_JBAT.txt out_JBAT.hic.part \
 <(cat out_JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}')) && \
(mv out_JBAT.hic.part out_JBAT.hic)

# =============================================================================
# STEP 4: POST-PROCESS MANUAL CURATION (if review.assembly exists)
# =============================================================================

echo "Step 4: Processing manual curation (if available)..."

if [[ -f "out_JBAT.review.assembly" ]]; then
    echo "  Found manual curation file, applying changes..."
    "${JUICER_YAHS}" post -o out_JBAT \
        out_JBAT.review.assembly \
        out_JBAT.liftover.agp \
        "${ORIGINAL_HAP2_CONTIGS}"
else
    echo "  No manual curation file found (out_JBAT.review.assembly)"
    echo "  Use Juicebox to manually curate, then rerun this section"
fi

# =============================================================================
# STEP 5: PREPARE ASSEMBLIES FOR COMPARISON
# =============================================================================

echo "Step 5: Preparing assemblies for comparison..."

# Sort and rename scaffolds by length for comparison
echo "  Sorting and renaming scaffolds..."
if [[ -f "${FINAL_HAP1}" ]]; then
    seqkit sort --by-length --reverse "${FINAL_HAP1}" | \
        seqkit replace --pattern '.+' --replacement 'Scaffold_{nr}' > \
        Atub_193_hap1_renamed.fasta
fi

if [[ -f "${FINAL_HAP2}" ]]; then
    seqkit sort --by-length --reverse "${FINAL_HAP2}" | \
        seqkit replace --pattern '.+' --replacement 'Scaffold_{nr}' > \
        Atub_193_hap2_renamed.fasta
fi

# =============================================================================
# STEP 6: SCAFFOLD ORIENTATION CORRECTION
# =============================================================================

echo "Step 6: Correcting scaffold orientations..."

# Create list of scaffolds to reverse (based on comparison analysis)
cat > scafstorev << 'EOF'
Scaffold_4
Scaffold_9
Scaffold_10
Scaffold_11
Scaffold_12
Scaffold_13
EOF

# Note: Based on D-GENIES comparison results:
# hap2 scaf7 = scaf4 (reverse)
# hap2 scaf4 = scaf5
# hap2 scaf10 = scaf7
# hap2 scaf13 = scaf9 (reverse)
# hap2 scaf11 = scaf10 (reverse)
# hap2 scaf5 = scaf11 (reverse)
# hap2 scaf9 = scaf12 (reverse)
# hap2 scaf12 = scaf13 (reverse)
# hap2 scaf16 = scaf14
# hap2 scaf14 = scaf16

# Extract scaffolds that need to be reversed, based on alignment of haplotypes to each other.
echo "  Extracting scaffolds to reverse..."
if [[ -f "Atub_193_hap2_renamed.fasta" ]]; then
    ./faSomeRecords Atub_193_hap2_renamed.fasta scafstorev scafstorev.fasta
    
    # Reverse complement the selected scaffolds
    echo "  Reverse complementing selected scaffolds..."
    revseq scafstorev.fasta > reversed.fasta
    
    # Extract scaffolds in correct orientation
    echo "  Extracting correctly oriented scaffolds..."
    ./faSomeRecords -exclude Atub_193_hap2_renamed.fasta scafstorev scafsinrightorder.fasta
    
    # Combine and sort final assembly
    echo "  Creating final ordered assembly..."
    cat reversed.fasta scafsinrightorder.fasta | \
        seqkit sort -n > Atub_193_hap2_ordered_v2.fasta
fi

# =============================================================================
# CLEANUP AND SUMMARY
# =============================================================================

echo "Pipeline completed successfully!"
echo ""
echo "Output files:"
echo "  - out.hic: Standard Hi-C contact matrix"
echo "  - out_JBAT.hic: Juicebox-compatible contact matrix"
echo "  - Atub_193_hap1_renamed.fasta: Renamed haplotype 1"
echo "  - Atub_193_hap2_ordered_v2.fasta: Final corrected haplotype 2"
echo ""
echo "Next steps:"
echo "  1. Use D-GENIES (https://dgenies.toulouse.inra.fr/) for assembly comparison"
echo "  2. Load out_JBAT.hic in Juicebox for manual curation if needed"
echo "  3. Validate assembly quality with additional QC tools"

# Clean up intermediate files (optional)
# rm -f alignments_sorted.txt.part out.hic.part out_JBAT.hic.part
# rm -f scafstorev.fasta reversed.fasta scafsinrightorder.fasta
