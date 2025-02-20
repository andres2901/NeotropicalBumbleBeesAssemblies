#!/bin/bash

# ###### Parameters area for resource request to SLURM ############################
#
#SBATCH --job-name=Assembly	#Job Name.
#SBATCH -p medium			#Partition/queue in which to run the job.
#SBATCH -N 1				#Number of nodes required.
#SBATCH -n 1				#Parallel task.
#SBATCH --mem=128G	#Memory per CPU
#SBATCH --time=168:00:00			#Maximum run time	
#SBATCH -o assembly.o%j			#Output file name
#
########################################################################################

file_path=$1 #path to the fastq file
dir_name=$2 #name for output directory

# ################## Module Loading Area ############################################

module load anaconda
source activate Flye

########################################################################################

flye --nano-raw $file_path --out-dir $dir_name --threads 32 --iterations 2 --genome-size 250m --no-alt-contigs

########################################################################################
