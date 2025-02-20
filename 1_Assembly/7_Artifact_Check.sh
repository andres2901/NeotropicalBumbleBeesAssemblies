#!/bin/bash
#!/bin/bash

# ###### Parameters area for resource request to SLURM ############################
#
#SBATCH --job-name=ArtifactCheck	#Job Name.
#SBATCH -p medium			#Partition/queue in which to run the job.
#SBATCH -N 1				#Number of nodes required.
#SBATCH -n 1				#Parallel task.
#SBATCH --mem=128G	#Memory per CPU
#SBATCH --time=12:00:00			#Maximum run time	
#SBATCH -o artifactCheck.o%j			#Output file name
#
########################################################################################

file_path=$1 #path to the assembly file

# ################## Module Loading Area ############################################

module load anaconda
source activate miniBUSCO

########################################################################################

compleasm run -t 8 -a $file_path -l hymenoptera -L /hpcfs/home/ciencias_biologicas/af.lizcano/mb_downloads/ -o artifact_check
########################################################################################
