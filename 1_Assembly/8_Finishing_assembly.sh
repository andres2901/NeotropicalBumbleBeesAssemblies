#!/bin/bash

# ###### Parameters area for resource request to SLURM ############################
#
#SBATCH --job-name=Finishing	#Job Name.
#SBATCH -p medium			#Partition/queue in which to run the job.
#SBATCH -N 1				#Number of nodes required.
#SBATCH -n 1				#Parallel task.
#SBATCH --mem=128G	#Memory per CPU
#SBATCH --time=168:00:00			#Maximum run time	
#SBATCH -o finishing.o%j			#Output file name
#
########################################################################################

fastq_path=$1 #path to the fastq file
assembly_path=$2 $path to assembly file
output_name=$2 #name for first output

# ################## Module Loading Area ############################################

module load anaconda
source activate ntlink

########################################################################################

#scaffolding and gap filling

ntLink scaffold gap_fill target=$assembly_path reads=$fastq_path k=40 w=500 t=5

mv ${assembly_path}.k40.w500.z1000.ntLink.scaffold.fa ntlink_round.fa

# racon polishing

conda deactivate
source activate minimap2

minimap2 -x map-ont --secondary=no -t 8 ntlink_round.fa $fastq_path > minimap.paf

conda deactivate
source activate racon

racon -q 20 -t 8 $fastq_path minimap.paf ntlink_round.fa > ntlink_racon.fa

#medaka
conda deactivate
source activate medaka

medaka_consensus -m r1041_e82_400bps_sup_v5.0.0 -i $fastq_path -d ntlink_racon.fa -t 8 -o medaka_Final_Assembly

rm minimap.paf
########################################################################################
