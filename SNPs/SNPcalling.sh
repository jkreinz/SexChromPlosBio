#!/bin/bash
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=60
#SBATCH --time=24:00:00
#SBATCH --job-name freebayes_parallel
cd $SLURM_SUBMIT_DIR

#turn off implicit threading
export OMP_NUM_THREADS=1

#module load nixpkgs/16.09
#module load gcc/7.3.0

module load CCEnv
module load nixpkgs/16.09
module load parallel/20160722
module load StdEnv/2020
module load freebayes/1.3.6

#run on multiple nodes
HOSTS=$(scontrol show hostnames $SLURM_NODELIST | tr '\n' ,)

hap=1
ref=193
refdir=/scratch/w/wrighste/kreinerj/2023/refs
bamdir=/scratch/w/wrighste/kreinerj/2023/expanded_snpcalling/bams
vcfdir=/scratch/w/wrighste/kreinerj/2023/expanded_snpcalling

#mkdir -p ${vcfdir}/${ref}_${hap}

#/home/w/wrighste/kreinerj/freebayes/scripts/fasta_generate_regions.py ${refdir}/Atub_${ref}_hap${hap}.fasta 100000 > ${refdir}/${ref}_${hap}_100kbregions.txt

#python /scratch/w/wrighste/kreinerj/2023/expanded_snpcalling/checkifcomplete.py  > checked_lenient
#bash parallel_bash.sh output.txt > checkedall
#grep ": Complete" checkedall | sed 's|/scratch/w/wrighste/kreinerj/2023/expanded_snpcalling/193_1/||g' | awk -F":" '{print $1":"$2}' | sed 's/193_1_//g' | sed 's/.vcf//g' > done1931
#grep -v -Ff done2 ${refdir}/${ref}_${hap}_100kbregions.txt > redo2


# iterate over regions using gnu parallel to dispatch jobs
cat ${vcfdir}/redo2 | parallel --env OMP_NUM_THREADS,PATH,LD_LIBRARY_PATH --joblog slurm-$SLURM_JOBID.log -k -j $SLURM_NTASKS_PER_NODE -S $HOSTS --wd $PWD "freebayes -f ${refdir}/Atub_${ref}_hap${hap}.fasta --use-best-n-alleles 2 --report-monomorphic ${bamdir}/*_${ref}_${hap}.dd.bam > ${vcfdir}/${ref}_${hap}/${ref}_${hap}_{}.vcf" --region {}
