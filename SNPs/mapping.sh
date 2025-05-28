# 1. Map with bwa
ref="/ohta2/julia.kreiner/hifi_hic_malegenomes/final_assemblies/Atub_193_hap2.fasta"
refout="193_2"
threads="4"
picard=/ohta/tyler.kent/Software/picard-2.8.1/picard.jar
prefix=$1
indir=/ohta2/julia.kreiner/hifi_hic_malegenomes/remapping/trimmedreads


bwa mem -t $threads $ref ${indir}/${prefix}_R1_trimmed.fastq.gz ${indir}/${prefix}_R2_trimmed.fastq.gz -R $(echo "@RG\tID:${prefix}\tCN:NOVOGENE\tPL:ILLUMINA\tPM:NOVASEQ.S4\tSM:${prefix}") | samtools view -hb -o ${prefix}_${refout}.bam -  2> mapping.log

# 2. Sort with samtools
samtools view -@ $threads -hF 4 ${prefix}_${refout}.bam | samtools sort -T /ohta2/julia.kreiner/hifi_hic_malegenomes/remapping/output/tmp -@ $threads - > ${prefix}_${refout}_sorted.bam

# 3. Mark Duplicates

java -Xmx15G -Djava.io.tmpdir=~/tmp -jar $picard MarkDuplicates I=${prefix}_${refout}_sorted.bam O=${prefix}_${refout}.dd.bam M=${prefix}_${refout}.dedup.log AS=true

#4. index bam
samtools index -@ $threads ${prefix}_${refout}.dd.bam
