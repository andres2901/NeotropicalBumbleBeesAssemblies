#!/bin/bash

# ###### Parameters area for resource request to SLURM ############################
#
#SBATCH --job-name=InitialQuality	#Job Name.
#SBATCH -p medium			#Partition/queue in which to run the job.
#SBATCH -N 1				#Number of nodes required.
#SBATCH -n 1				#Parallel task.
#SBATCH --mem=32G	#Memory per CPU
#SBATCH --time=2:00:00			#Maximum run time	
#SBATCH -o initial_quality.o%j			#Output file name
#
########################################################################################

file_path=$1 #path to the fastq file

# ################## Module Loading Area ############################################

module load anaconda
source activate nano

########################################################################################

NanoPlot -t 8 --tsv_stats --fastq $file_path

########################################################################################
