#!/bin/bash

# ###### Parameters area for resource request to SLURM ############################
#
#SBATCH --job-name=NonCoding	#Job Name.
#SBATCH -p medium			#Partition/queue in which to run the job.
#SBATCH -N 1				#Number of nodes required.
#SBATCH -n 1				#Parallel task.
#SBATCH --mem=64G	#Memory per CPU
#SBATCH --time=72:00:00			#Maximum run time	
#SBATCH -o nonCoding.o%j			#Output file name
#
########################################################################################

file_path=$1 #path to the assembly file
output_name=$2 #name for output

# ################## Module Loading Area ############################################

module load anaconda
source activate trnascan

########################################################################################

tRNAscan-SE -o ${output_name}_tRNAscan.tsv -f ${output_name}_tRNAscanStructure -m ${output_name}_tRNAscanStatistics.tsv -j ${output_name}_tRNAscan.gff3 ${file_path}.fa

## run infernal

conda deactivate
source activate infernal

cmscan --rfam --cut_ga --nohmmonly --tblout ${output_name}_infernal.tblout --fmt 2 --clanin /hpcfs/home/ciencias_biologicas/af.lizcano/database/Rfam/Rfam.clanin /hpcfs/home/ciencias_biologicas/af.lizcano/database/Rfam/Rfam.cm ${file_path}.fa > ${output_name}_infernal.cmscan
########################################################################################
