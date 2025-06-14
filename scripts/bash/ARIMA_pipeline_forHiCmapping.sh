#! /bin/bash

##############################################
# ARIMA GENOMICS MAPPING PIPELINE 02/08/2019 #
##############################################

#Below find the commands used to map HiC data.

#Replace the variables at the top with the correct paths for the locations of files/programs on your system.

#This bash script will map one paired end HiC dataset (read1 & read2 fastqs). Feel to modify and multiplex as you see fit to work with your volume of samples and system.

##########################################
# Commands #
##########################################

SRA='Ama_NAT_HiC'
LABEL='s45'
BWA='bwa'
SAMTOOLS='/usr/bin/samtools'
IN_DIR='/ohta2/julia.kreiner/hifi_hic_malegenomes/HiC/Ama_N5_HiC_DSN/'
REF='/ohta2/julia.kreiner/hifi_hic_malegenomes/hifiasm_run2_hets45/19_3_ATUB_s45.asm.hic.hap1.p_ctg.fasta'
FAIDX='$REF.fai'
PREFIX='19_3_ATUB_s45.asm.hic.hap1.p_ctg.fasta'
RAW_DIR='/ohta2/julia.kreiner/hifi_hic_malegenomes/hifiasm_run2_hets45/bams'
FILT_DIR='/ohta2/julia.kreiner/hifi_hic_malegenomes/hifiasm_run2_hets45/bams'
FILTER='/ohta2/julia.kreiner/hifi_hic_malegenomes/scripts/filter_five_end.pl'
COMBINER='/ohta2/julia.kreiner/hifi_hic_malegenomes/scripts/two_read_bam_combiner.pl'
STATS='/ohta2/julia.kreiner/hifi_hic_malegenomes/scripts/get_stats.pl'
PICARD='/ohta1/julia.kreiner/software/picard.jar'
TMP_DIR='/ohta2/julia.kreiner/hifi_hic_malegenomes/hifiasm_run2_hets45/tmp'
PAIR_DIR='/ohta2/julia.kreiner/hifi_hic_malegenomes/hifiasm_run2_hets45/bams'
REP_DIR='/ohta2/julia.kreiner/hifi_hic_malegenomes/hifiasm_run2_hets45/bams'
REP_LABEL=$LABEL\_rep1
MERGE_DIR='/ohta2/julia.kreiner/hifi_hic_malegenomes/hifiasm_run2_hets45/bams/final'
MAPQ_FILTER=10
CPU=12

echo "### Step 0: Check output directories exist & create them as needed"
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR
[ -d $REP_DIR ] || mkdir -p $REP_DIR
[ -d $MERGE_DIR ] || mkdir -p $MERGE_DIR

echo "### Step 0: Index reference" # Run only once! Skip this step if you have already generated BWA index files
#$BWA index -a bwtsw -p $PREFIX $REF

echo "### Step 1.A: FASTQ to BAM (1st)"
#$BWA mem -t $CPU $REF $IN_DIR/$SRA\_1.fq.gz | $SAMTOOLS view -@ $CPU -Sb - > $RAW_DIR/$SRA\_1.bam

echo "### Step 1.B: FASTQ to BAM (2nd)"
#$BWA mem -t $CPU $REF $IN_DIR/$SRA\_2.fq.gz | $SAMTOOLS view -@ $CPU -Sb - > $RAW_DIR/$SRA\_2.bam

echo "### Step 2.A: Filter 5' end (1st)"
$SAMTOOLS view -h $RAW_DIR/$SRA\_1.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/$SRA\_1.bam

echo "### Step 2.B: Filter 5' end (2nd)"
$SAMTOOLS view -h $RAW_DIR/$SRA\_2.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/$SRA\_2.bam

echo "### Step 3A: Pair reads & mapping quality filter"
perl $COMBINER ${SRA}_filt_1.bam ${SRA}_filt_2.bam $SAMTOOLS $MAPQ_FILTER | $SAMTOOLS view -bS -t $FAIDX - | $SAMTOOLS sort -@ $CPU -o ${SRA}_filtsort.bam -

echo "### Step 3.B: Add read group"
java -Xmx4G -Djava.io.tmpdir=temp/ -jar $PICARD AddOrReplaceReadGroups INPUT=${SRA}_filtsort.bam OUTPUT=${SRA}_filtsortrg.bam ID=$SRA LB=$SRA SM=$LABEL PL=ILLUMINA PU=none

###############################################################################################################################################################
###                                           How to Accommodate Technical Replicates                                                                       ###
### This pipeline is currently built for processing a single sample with one read1 and read2 fastq file.                                                    ###
### Technical replicates (eg. one library split across multiple lanes) should be merged before running the MarkDuplicates command.                          ###
### If this step is run, the names and locations of input files to subsequent steps will need to be modified in order for subsequent steps to run correctly.###
### The code below is an example of how to merge technical replicates.                                                                                      ###
###############################################################################################################################################################
#       REP_NUM=X #number of the technical replicate set e.g. 1
#       REP_LABEL=$LABEL\_rep$REP_NUM
#       INPUTS_TECH_REPS=('bash' 'array' 'of' 'bams' 'from' 'replicates') #BAM files you want combined as technical replicates
#   example bash array - INPUTS_TECH_REPS=('INPUT=A.L1.bam' 'INPUT=A.L2.bam' 'INPUT=A.L3.bam')
#       java -Xmx8G -Djava.io.tmpdir=temp/ -jar $PICARD MergeSamFiles $INPUTS_TECH_REPS OUTPUT=$TMP_DIR/$REP_LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT

echo "### Step 4: Mark duplicates"
java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $PICARD MarkDuplicates INPUT=${SRA}_filtsortrg.bam OUTPUT=${SRA}_filtsortrgdedup.bam METRICS_FILE=$REP_DIR/metrics.$REP_LABEL.txt TMP_DIR=$TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

$SAMTOOLS index ${SRA}_filtsortrgdedup.bam

perl $STATS ${SRA}_filtsortrgdedup.bam > ${SRA}_filtsortrgdedup.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"

#########################################################################################################################################
###                                       How to Accommodate Biological Replicates                                                    ###
### This pipeline is currently built for processing a single sample with one read1 and read2 fastq file.                              ###
### Biological replicates (eg. multiple libraries made from the same sample) should be merged before proceeding with subsequent steps.###
### The code below is an example of how to merge biological replicates.                                                               ###
#########################################################################################################################################
#
#       INPUTS_BIOLOGICAL_REPS=('bash' 'array' 'of' 'bams' 'from' 'replicates') #BAM files you want combined as biological replicates
#   example bash array - INPUTS_BIOLOGICAL_REPS=('INPUT=A_rep1.bam' 'INPUT=A_rep2.bam' 'INPUT=A_rep3.bam')
#
#       java -Xmx8G -Djava.io.tmpdir=temp/ -jar $PICARD MergeSamFiles $INPUTS_BIOLOGICAL_REPS OUTPUT=$MERGE_DIR/$LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT
#
#       $SAMTOOLS index $MERGE_DIR/$LABEL.bam

# perl $STATS $MERGE_DIR/$LABEL.bam > $MERGE_DIR/$LABEL.bam.stats

# echo "Finished Mapping Pipeline through merging Biological Replicates"
