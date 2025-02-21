#!/bin/bash

# ###### Parameters area for resource request to SLURM ############################
#
#SBATCH --job-name=Coverage	#Job Name.
#SBATCH -p medium			#Partition/queue in which to run the job.
#SBATCH -N 1				#Number of nodes required.
#SBATCH -n 1				#Parallel task.
#SBATCH --mem=128G	#Memory per CPU
#SBATCH --time=168:00:00			#Maximum run time	
#SBATCH -o coverage.o%j			#Output file name
#
########################################################################################

fastq_path=$1 #path to the fastq file
assembly_path=$2 #path to the assembly file

# ################## Module Loading Area ############################################

module load anaconda
source activate minimap2

########################################################################################

minimap2 -x map-ont -a --sam-hit-only --secondary=no -t 8 $assembly_path $fastq_path > minimap_coverage.sam

conda deactivate
source activate samtools

samtools sort -@16 -O BAM -o coverage_sorted.bam minimap_coverage.sam
rm minimap_coverage*
########################################################################################
