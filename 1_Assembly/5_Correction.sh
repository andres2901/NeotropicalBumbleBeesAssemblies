#!/bin/bash

# ###### Parameters area for resource request to SLURM ############################
#
#SBATCH --job-name=Correction	#Job Name.
#SBATCH -p medium			#Partition/queue in which to run the job.
#SBATCH -N 1				#Number of nodes required.
#SBATCH -n 1				#Parallel task.
#SBATCH --mem=128G	#Memory per CPU
#SBATCH --time=168:00:00			#Maximum run time	
#SBATCH -o correction.o%j			#Output file name
#
########################################################################################

fastq_path=$1 #path to the fastq file
assembly_path=$2 #path to the assembly file
dir_name=$3 #name for medaka directory
draft_name=$4 #name for the output

# ################## Module Loading Area ############################################

module load anaconda
source activate medaka

# medaka correction

medaka_consensus -m r1041_e82_400bps_sup_v5.0.0 -i $fastq_path -d $assembly_path -t 8 -o $dir_name

#filter by length of contigs

conda deactivate
source activate seqkit

seqkit seq -m 50000 ${dir_name}/consensus.fasta > $drat_name

#verify completeness of the first draft

conda deactivate
source activate miniBUSCO

compleasm run -t 8 -a $drat_name -l hymenoptera -L /hpcfs/home/ciencias_biologicas/af.lizcano/mb_downloads/ -o compleasm_first_draft

#verify quality of first draft

conda deactivate
module load quast

quast -o quast_first_draft $drat_name
########################################################################################
