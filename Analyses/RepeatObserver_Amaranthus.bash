############################################
# Run the Amaranthus genome
# April 2024
# Cassandra Elphinstone
###########################################

cd /home/celphin/scratch/repeats/auto_script/

# download the Setup_Run_Repeats.sh from github
wget https://raw.githubusercontent.com/celphin/RepeatOBserverV1/main/Setup_Run_Repeats.sh
chmod +x Setup_Run_Repeats.sh
dos2unix Setup_Run_Repeats.sh

# install the R package
# already done here
install.packages("devtools")
library(devtools)
install_github("celphin/RepeatOBserverV1") #to install the package
# Select 1:All to install all the required packages
library(RepeatOBserverV1) # to load the package

#------------------------------
# copy over the Amaranthus haplotypes
scp *.fasta celphin@beluga.computecanada.ca:/home/celphin/scratch/repeats/auto_script/

#---------------------------------
# run
cat << EOF > Auto_Atub193_hap1.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Atub193 -f Atub_193_hap1.fasta -h H1 -c 20 -m 191000M -g FALSE

EOF

sbatch Auto_Atub193_hap1.sh

#---------------------------------
# run
cat << EOF > Auto_Atub193_hap2.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Atub193 -f Atub_193_hap2.fasta -h H2 -c 20 -m 191000M -g FALSE

EOF

sbatch Auto_Atub193_hap2.sh

#---------------------------------
# run
cat << EOF > Auto_AtubNat_hap1.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i AtubNat -f Atub_Nat_hap1.fasta -h H1 -c 20 -m 191000M -g FALSE

EOF

sbatch Auto_AtubNat_hap1.sh

#---------------------------------
# run
cat << EOF > Auto_AtubNat_hap2.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i AtubNat -f Atub_Nat_hap2.fasta -h H2 -c 20 -m 191000M -g FALSE

EOF

sbatch Auto_AtubNat_hap2.sh

#---------------------------------
# run
cat << EOF > Auto_AtubNune5_hap1.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i AtubNune5 -f Atub_Nune5_hap1.fasta -h H1 -c 20 -m 191000M -g FALSE

EOF

sbatch Auto_AtubNune5_hap1.sh

#---------------------------------
# run
cat << EOF > Auto_AtubNune5_hap2.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i AtubNune5 -f Atub_Nune5_hap2.fasta -h H2 -c 20 -m 191000M -g FALSE

EOF

sbatch Auto_AtubNune5_hap2.sh


#------------------------------
sq
          # JOBID     USER      ACCOUNT           NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON)
       # 46658246  celphin def-rieseber Auto_AtubNune5  PD   23:00:00     1   20        N/A 191000M  (Priority)
       # 46658247  celphin def-rieseber Auto_AtubNune5  PD   23:00:00     1   20        N/A 191000M  (Priority)
       # 46658248  celphin def-rieseber Auto_AtubNat_h  PD   23:00:00     1   20        N/A 191000M  (Priority)
       # 46658249  celphin def-rieseber Auto_AtubNat_h  PD   23:00:00     1   20        N/A 191000M  (Priority)
       # 46658250  celphin def-rieseber Auto_Atub193_h  PD   23:00:00     1   20        N/A 191000M  (Priority)
       # 46658251  celphin def-rieseber Auto_Atub193_h  PD   23:00:00     1   20        N/A 191000M  (Priority)

#------------------------------
# copy over 

# to cedar shared
scp -r ./Atub*_H*-AT/Summary_output celphin@cedar.computecanada.ca:/home/celphin/projects/def-rieseber/Dryas_shared_data/Amaranthus

# to local
scp celphin@cedar.computecanada.ca:/home/celphin/projects/def-rieseber/Dryas_shared_data/Amaranthus/*.png .
scp -r celphin@beluga.computecanada.ca:/home/celphin/scratch/repeats/auto_script/output_chromosomes/Atub*_H*-AT/Summary_output/spectra/spectra_parts_35-2000 .

##########################################
# copy over outgroup genomes from google drive link
cd /home/celphin/scratch/repeats/Amaranthus

# rename
mv Amacr_genome.fasta Amacr.fasta
mv Amaranthus_hypochondriacus.faa Amahypochondriacus.fasta 
mv Amaranthus_palmeri.faa Amapalmeri.fasta 
mv A.tricolor.red.chromosomes.fa Amatricolor.fasta

#------------------------------------------
# run
cat << EOF > Auto_Amacr.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=11:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=128000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Amacr -f Amacr.fasta -h H0 -c 15 -m 128000M -g FALSE

EOF

sbatch Auto_Amacr.sh

#----------------------------------------
cat << EOF > Auto_Amahypochondriacus.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=11:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=128000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Amahypochondriacus -f Amahypochondriacus.fasta -h H0 -c 15 -m 128000M -g FALSE

EOF

sbatch Auto_Amahypochondriacus.sh

#-----------------------------------------------------
cat << EOF > Auto_Amatricolor.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=11:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=128000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Amatricolor -f Amatricolor.fasta -h H0 -c 15 -m 128000M -g FALSE

EOF

sbatch Auto_Amatricolor.sh

#------------------------------------------
cat << EOF > Auto_Amapalmeri.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=11:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=128000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Amapalmeri -f Amapalmeri.fasta -h H0 -c 15 -m 128000M -g FALSE

EOF

sbatch Auto_Amapalmeri.sh


###############################
# copy over summary files
cp -v -u -r /home/celphin/scratch/repeats/Amaranthus/output_chromosomes/*/Summary_output /home/celphin/scratch/repeats/Amaranthus/Summary


###############################
# search for centromeric repeats
cd /home/celphin/scratch/repeats/auto_script/EMBOSS_repeats

# To run this script EMBOSS6.6.0 is required. 
module load nixpkgs/16.09  intel/2018.3  emboss/6.6.0

# download the script
wget https://raw.githubusercontent.com/celphin/RepeatOBserverV1/main/repeat_seq_finder.sh
chmod +x repeat_seq_finder.sh
dos2unix repeat_seq_finder.sh

# uncomment seqkit load
nano repeat_seq_finder.sh

# direct the script to your chromosome files: "/home/celphin/scratch/repeats/auto_script/input_chromosomes/NerLuet_H0-AT/chromosome_files"
# Choose the Chromosome fasta file that you want to look at the repeat in: "NerLuet_H0-AT_Chr1part01.fasta"
# Define the bp range: 13200000 13500000
# Define the repeat length you are looking for: 4 (will work for repeats up to 600bp long)
# Define the number of bp per line in your fasta file: 60 (default for RepeatOBserverV1)

./repeat_seq_finder.sh "/home/celphin/scratch/repeats/auto_script/input_chromosomes/Atub193_H2-AT/chromosome_files" \
"Atub193_H2-AT_Chr8part01.fasta" 28000000 29000000 600 60 50

#Atub193_H2-AT_Chr8part01.fasta_596_914279-915429
CACGCTTAACTGCGGAGTTCTGATGGGATCCGTTGCATTAGTGCTGGTATGATCGCACCTATTCATTCACTGTAACAATTTTTTTATATGCGCAATGCCGCCCATTAAAAACTCCAAGCATTCCAAAGCTCGTAAATTCGCAACCGAGACCCCGATCTCGAGCCTCCTTCTTCCCTAATTCTGTTTTTTCCCCGTCCCGGTATTGTCCGAATCACTAAAGAGTCCGCCGGCATTAAAAGCCCGGTGAACTCGCAACAAAAAAGCGACAAAATATTTAAGAAATCCGCGCAACACTTGGAAATCAAGAGGCCATCCATGCCAGGCACGCTTCCGCGGGGAGTTCCGACGAGGTCCGGTGCAAATTCTGCATTCATGAACGCACCGATGCACGCTTCTTTCGAACAATTCGTTCGTATGCGCATCGACCCACCTCAAGGGATCCTAATTTTCGAGGGCTCGTAAATTCGAGGTCCGGACCCGGTTGCGAGTTATAATTTTAAAATTAAATTATTTCCAACTAAAATATTATTTTTGTATTCACTAACATGTCCTTTACTCGAAAATAGCGCAACGAATTCGAAACAAAAAAAAAATTA
AAAAACGGGAAAAAGGGGTGCAACACGAGGACTTCCCAGGAGGTCACCCATCCTAGTACTACTCTCGCCCTGTAACAATTTTTTTATATGCGCAATGCCGCCCATTAAAAACTCCAAGCATTCCAAAGCTCGTAAATTCGCAACCGAGACCCCGATCTCGAGCCTCCTTCTTCCCTAATTCTGTTTTTTCCCCGTCCCGGTATTGTCCGAATCACTAAAGAGTCCGCCGGCATTAAAAGCCCGGTGAACTCGCAACAAAAAAGCGACAAAATATTTAAGAAATCCGCGCAACACTTGGAAATCAAGAGGCCATCCATGCCAGGCACGCTTCCGCGGGGAGTTCCGACGCGGTCCGGTGCAAATTCTGCATTCATGAACGCACCGATGCACGCTTCTTTCGAACAATTCGTTCGTATGCGCATCGACCCACCTCAAGGGATCCTAATTTTCGAGGGCTCGTAAATTCGAGGTCGGAACCCGGTTGCGAGTTATAATTTTAAAATTAATTTATTTCCAACTAAAATATTATTTTTGTATTCACTAACATGCCCTTTACTCGAAAATAGCGCAACGAATTCGAAACAAAAAAAAATTAA
AAAACGGG

#Atub193_H2-AT_Chr8part01.fasta_595_914805-915481
CTAACATGTCCTTTACTCGAAAATAGCGCAACGAATTCGAAACAAAAAAAAAATTAAAAAACGGGAAAAAGGGGTGCAACACGAGGACTTCCCAGGAGGTCACCCATCCTAGTACTACTCTCGCCCTGTAACAATTTTTTTATATGCGCAATGCCGCCCATTAAAAACTCCAAGCATTCCAAAGCTCGTAAATTCGCAACCGAGACCCCGATCTCGAGCCTCCTTCTTCCCTAATTCTGTTTTTTCCCCGTCCCGGTATTGTCCGAATCACTAAAGAGTCCGCCGGCATTAAAAGCCCGGTGAACTCGCAACAAAAAAGCGACAAAATATTTAAGAAATCCGCGCAACACTTGGAAATCAAGAGGCCATCCATGCCAGGCACGCTTCCGCGGGGAGTTCCGACGCGGTCCGGTGCAAATTCTGCATTCATGAACGCACCGATGCACGCTTCTTTCGAACAATTCGTTCGTATGCGCATCGACCCACCTCAAGGGATCCTAATTTTCGAGGGCTCGTAAATTCGAGGTCGGAACCCGGTTGCGAGTTATAATTTTAAAATTAATTTATTTCCAACTAAAATATTATTTTTGTATTC
ACTAACATGCCCTTTACTCGAAAATAGCGCAACGAATTCGAAACAAAAAAAAATTAAAAAACGGGAAAAAGGGGTGCAACACGAGGACTTCCCAGGAGGTCACCCATCCTAGTACTACTCTCGCC


######################################
