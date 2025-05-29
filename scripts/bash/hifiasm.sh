prefix="Nune"
pathtohic="/ohta2/julia.kreiner/hifi_hic_malegenomes/HiC/Ama_N5_HiC_DSN/"
prefix="Ama_N5_HiC_DSN"
hifipath="/ohta2/julia.kreiner/hifi_hic_malegenomes/HiFi_Data/Nune_rundata"
hifi_pre="Nune5_all"

#run hifiasm with both hic data and pacbio data
hifiasm -o ${prefix}.asm -t32 --h1 ${pathtohic}/${prefix}_1.fq.gz --h2 ${pathtohic}/${prefix}_2.fq.gz ${hifipath}/${hifi_pre}.fastq.gz &

#convert each resulting haplome from gfa to fasta format
awk '/^S/{print ">"$2"\n"$3}' Ama_N5_HiC_DSN_s40.asm.hic.hap1.p_ctg.gfa | fold > Ama_N5_HiC_s40_hap1.fasta &

awk '/^S/{print ">"$2"\n"$3}' Ama_N5_HiC_DSN_s40.asm.hic.hap2.p_ctg.gfa | fold > Ama_N5_HiC_s40_hap2.fasta

