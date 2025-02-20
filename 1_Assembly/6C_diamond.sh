#!/bin/bash

# ###### Parameters area for resource request to SLURM ############################
#
#SBATCH --job-name=Diamond	#Job Name.
#SBATCH -p medium			#Partition/queue in which to run the job.
#SBATCH -N 1				#Number of nodes required.
#SBATCH -n 1				#Parallel task.
#SBATCH --mem=128G	#Memory per CPU
#SBATCH --time=168:00:00			#Maximum run time	
#SBATCH -o diamond.o%j			#Output file name
#
########################################################################################

file_path=$1 #path to the assembly file

# ################## Module Loading Area ############################################

module load anaconda
source activate diamond

########################################################################################

diamond blastx --query $file_path --db /hpcfs/home/ciencias_biologicas/af.lizcano/uniprot/reference_proteomes.dmnd --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --sensitive --max-target-seqs 1 --evalue 1e-25 --threads 16 > diamond.out
########################################################################################
