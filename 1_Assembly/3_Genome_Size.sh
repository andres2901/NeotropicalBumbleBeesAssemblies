#!/bin/bash

# ###### Parameters area for resource request to SLURM ############################
#
#SBATCH --job-name=InitialQuality	#Job Name.
#SBATCH -p short			#Partition/queue in which to run the job.
#SBATCH -N 1				#Number of nodes required.
#SBATCH -n 1				#Parallel task.
#SBATCH --mem=32G	#Memory per CPU
#SBATCH --time=2:00:00			#Maximum run time	
#SBATCH -o genomeSize.o%j			#Output file name
#
########################################################################################

file_path=$1 #path to the fastq file

# ################## Module Loading Area ############################################

module load anaconda
source activate kat

########################################################################################

kat hist -o k21 -t 8 -m 21 $file_path
kat hist -o k27 -t 8 -m 27 $file_path
kat hist -o k31 -t 8 -m 31 $file_path
