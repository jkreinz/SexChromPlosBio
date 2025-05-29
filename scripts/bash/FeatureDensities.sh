#!/bin/bash

#Script for getting the densities of different types of repeats in 10kb windows across the genome
#################
#get densities
#################
reference_fai=/ohta2/julia.kreiner/hifi_hic_malegenomes/final_assemblies/Atub_193_hap2.fasta.fai
#I am going to use the split gff file from EDTA for this:
RM_out=/ohta2/julia.kreiner/hifi_hic_malegenomes/EDTA/193/hap2/Atub_193_hap2.fasta.mod.EDTA.intact.gff3
ref_anno=/ohta2/julia.kreiner/hifi_hic_malegenomes/final_assemblies/Atub_193_hap2.all.gff
pregenome="193_2"
bedopds=/ohta2/julia.kreiner/hifi_hic_malegenomes/densities/bin

#get start and end positions from fasta index
awk '{print $1 "\t" 0 "\t" $2}' $reference_fai > ${pregenome}_startend.fai.tmp
#chop scaffolds onto 10000 bp windows
bedops --chop 100000 ${pregenome}_startend.fai.tmp > ${pregenome}_chopped.txt.tmp
#sort the scaffold windows
bedtools sort -i ${pregenome}_chopped.txt.tmp > ${pregenome}_chopped_sorted.txt.tmp

#get genes
awk '($3 == "gene") {print $1 "\t" $4-1 "\t" $5}' $ref_anno | bedtools sort > ${pregenome}_genes.bed

#get TEs
awk '{print $1 "\t" $4-1 "\t" $5}' $RM_out | bedtools sort > ${pregenome}_TEs.bed

#find overlap between reference window file and data position file
bin/bedmap --delim '\t' --echo --count --bases-uniq-f ${pregenome}_chopped_sorted.txt.tmp ${pregenome}_genes.bed > ${pregenome}_genes_100kbdensity.bed &
bin/bedmap --delim '\t' --echo --count --bases-uniq-f ${pregenome}_chopped_sorted.txt.tmp ${pregenome}_TEs.bed > ${pregenome}_TEs_100kbdensity.bed &
