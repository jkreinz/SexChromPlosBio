#!/bin/bash

# Define input and output variables
input_vcf="$1"
output_dir="/ohta2/julia.kreiner/hifi_hic_malegenomes/remapping/commongarden_snps"
site_type="$2"

# Define directories
LOG_DIR="/ohta2/julia.kreiner/hifi_hic_malegenomes/remapping/commongarden_snps"

# Define threads and resources
threads=8
#mem_mb=16000 # 8 * 2000
#time="06:00:00"

# Create output directory if it doesn't exist
#mkdir -p "$output_dir"

# Define log file
log="$LOG_DIR/bcftools_split_variants_${site_type}.log"

# Define temporary directory
prefix=`echo $input_vcf | sed 's/_unfiltsnps.vcf.gz//g'`
mkdir -p $output_dir/tmp/$prefix
tmpdir="$output_dir/tmp/$prefix"

echo $prefix

if [ "$site_type" = "invariant" ]; then
    # Filter invariant sites
    bcftools view --threads "$threads" -O z --include 'N_ALT = 0' "$input_vcf" | \
    bcftools filter -i 'F_MISSING <= 0.25' -Oz >  "$output_dir/${prefix}_filt${site_type}.vcf.gz" && sleep 5 && tabix $output_dir/${prefix}_filt${site_type}.vcf.gz
elif [ "$site_type" = "snps" ]; then
    # Filter SNPs
   # (bcftools view --threads "$threads" -O v --types "$site_type" "$input_vcf" | \
   # vcfallelicprimitives --keep-info --keep-geno | \
   # bcftools view --threads "$threads" --types "$site_type" --min-alleles 2 --max-alleles 2 | \
   # bcftools sort -O z -T "$tmpdir" -o "$output_dir/${prefix}_unfilt${site_type}.vcf.gz" )
   (bcftools filter -i 'F_MISSING <= 0.25' $output_dir/${prefix}_unfilt${site_type}.vcf.gz | \
                bcftools filter -i 'QUAL >= 30' | \
                bcftools filter -i 'AB >= 0.25 & AB <= 0.75 | AB <= 0.01' | \
                bcftools filter -i 'SAF > 0 & SAR > 0' | \
                bcftools filter -i 'MQM >=30 & MQMR >= 30' | \
                bcftools filter -i '((PAIRED > 0.05) & (PAIREDR > 0.05) & (PAIREDR / PAIRED < 1.75 ) & (PAIREDR / PAIRED > 0.25)) | ((PAIRED < 0.05) & (PAIREDR < 0.05))' |
                bcftools filter -o "$output_dir/${prefix}_filtered${site_type}.vcf.gz" -O z -i '((AF > 0) & (AF < 1))' && sleep 5 && tabix "$output_dir/${prefix}_filtered${site_type}.vcf.gz")
fi
