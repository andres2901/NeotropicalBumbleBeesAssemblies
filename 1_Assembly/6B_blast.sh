#!/bin/bash

# ###### Parameters area for resource request to SLURM ############################
#
#SBATCH --job-name=Blast	#Job Name.
#SBATCH -p medium			#Partition/queue in which to run the job.
#SBATCH -N 1				#Number of nodes required.
#SBATCH -n 1				#Parallel task.
#SBATCH --mem=256G	#Memory per CPU
#SBATCH --time=168:00:00			#Maximum run time	
#SBATCH -o blast.o%j			#Output file name
#
########################################################################################

file_path=$1 #path to the assembly file

# ################## Module Loading Area ############################################

module load anaconda
source activate blast

########################################################################################

blastn -db /hpcfs/home/ciencias_biologicas/af.lizcano/nt/core_nt -query $file_path -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -num_threads 16 -out blast.out

########################################################################################
