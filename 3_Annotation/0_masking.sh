#!/bin/bash

# ###### Parameters area for resource request to SLURM ############################
#
#SBATCH --job-name=Masking	#Job Name.
#SBATCH -p long			#Partition/queue in which to run the job.
#SBATCH -N 1				#Number of nodes required.
#SBATCH -n 1				#Parallel task.
#SBATCH --mem=256G	#Memory per CPU
#SBATCH --time=360:00:00			#Maximum run time	
#SBATCH -o masking.o%j			#Output file name
#
########################################################################################

file_path=$1 #path to the assembly file
dir_name=$2 #name for the TransposonUltimate directory
output_prefix=$3 #prefix for the masked output assembly

# ################## Module Loading Area ############################################

module load anaconda
source activate transposon_annotation_tools_env

#Create project

reasonaTE -mode createProject -projectFolder workspace -projectName $dir_name -inputFasta $file_path

#run different tools

reasonaTE -mode annotate -projectFolder workspace -projectName $dir_name -tool helitronScanner
reasonaTE -mode annotate -projectFolder workspace -projectName $dir_name -tool ltrHarvest
reasonaTE -mode annotate -projectFolder workspace -projectName $dir_name -tool mitefind
reasonaTE -mode annotate -projectFolder workspace -projectName $dir_name -tool mitetracker
reasonaTE -mode annotate -projectFolder workspace -projectName $dir_name -tool sinefind
reasonaTE -mode annotate -projectFolder workspace -projectName $dir_name -tool sinescan xxxxx -k 8
reasonaTE -mode annotate -projectFolder workspace -projectName $dir_name -tool tirvish
reasonaTE -mode annotate -projectFolder workspace -projectName $dir_name -tool transposonPSI
reasonaTE -mode annotate -projectFolder workspace -projectName $dir_name -tool NCBICDD1000

#run repeatMasker

conda deactivate
source activate repeatMasker

/hpcfs/home/ciencias_biologicas/af.lizcano/RepeatMasker/RepeatMasker -species Hymenoptera -e abblast -dir workspace/${dir_name}/repMasker/ -pa 4 workspace/${dir_name}/sequence.fasta

#run repeatModeler

conda deactivate
source activate repeatModeler

/hpcfs/home/ciencias_biologicas/af.lizcano/RepeatModeler-2.0.5/BUildDatabase -name sequence_index workspace/${dir_name}/sequence.fasta
mv sequence_index* workspace/${dir_name}/repeatmodel/

/hpcfs/home/ciencias_biologicas/af.lizcano/RepeatModeler-2.0.5/RepeatModeler -database workspace/${dir_name}/repeatmodel/sequence_index -threads 8

#run Must

conda deactivate
source activate must

perl /hpcfs/home/ciencias_biologicas/af.lizcano/MUST.r2-4-002.Release/MUST_Pipe.pl workspace/${dir_name}/sequence.fasta workspace/${dir_name}/must/result.txt workspace/${dir_name}/must/temp/ 8 8 50 2 30 100 600 50 0.80

#parse tools results

conda deactivate
source activate transposon_annotation_tools_env

reasonaTE -mode parseAnnotations -projectFolder workspace -projectName ${dir_name}

#generate final results

conda deactivate
source activate transposon_annotation_reasonaTE

reasonaTE -mode pipeline -projectFolder workspace -projectName ${dir_name}
reasonaTE -mode statistics -projectFolder workspace -projectName ${dir_name}

#mask assembly

conda deactivate
source activate bedtools

bedtools maskfasta -soft -fi $file_path -bed workspace/${dir_name}/finalResults/FinalAnnotations_TransposonMask.gff3 -fo $output_prefix
########################################################################################
