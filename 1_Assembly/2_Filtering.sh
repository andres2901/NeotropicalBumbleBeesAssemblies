#!/bin/bash

# ###### Parameters area for resource request to SLURM ############################
#
#SBATCH --job-name=Filtering	#Job Name.
#SBATCH -p medium			#Partition/queue in which to run the job.
#SBATCH -N 1				#Number of nodes required.
#SBATCH -n 1				#Parallel task.
#SBATCH --mem=256	#Memory per CPU
#SBATCH --time=72:00:00			#Maximum run time	
#SBATCH -o filtering.o%j			#Output file name
#
########################################################################################

file_path=$1 #path to the fastq file
output_name=$2 #name for first output

# ################## Module Loading Area ############################################

module load anaconda
source activate porechop

########################################################################################

porechop --discard_middle -t 10 -i $file_path -o ${output_name}.fastq.gz

source deactivate
conda activate seqkit

seqkit rmdup ${output_name}.fastq.gz | gzip -c > reads_rmdup.fastq.gz #remove duplicates if exist
seqkit seq -m 200 -Q 20 reads_rmdup.fastq.gz | gzip -c > Dataset1.fastq.gz # generate first dataset
seqkit seq -m 10000 Dataset1.fastq.gz | gzip -c > Dataset2.fastq.gz # generate second dataset
########################################################################################
