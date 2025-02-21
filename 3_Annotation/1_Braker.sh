#!/bin/bash

# ###### Parameters area for resource request to SLURM ############################
#
#SBATCH --job-name=Braker	#Job Name.
#SBATCH -p medium			#Partition/queue in which to run the job.
#SBATCH -N 1				#Number of nodes required.
#SBATCH -n 1				#Parallel task.
#SBATCH --mem=128G	#Memory per CPU
#SBATCH --time=168:00:00			#Maximum run time	
#SBATCH -o braker.o%j			#Output file name
#
########################################################################################

assembly_path=$1 #path to the assembly file
output_name=$2 #name for output directory
rnaseq_dir=$3 #path of the directory of rnaseq data
species_name=$4 #name for the species in Braker and prefix for the IDs in the final gff3 file

# ################## Module Loading Area ############################################

module load anaconda
source activate braker3

########################################################################################

export GENEMARK_PATH=/hpcfs/home/ciencias_biologicas/af.lizcano/GeneMark-ETP/bin/
export PROTHINT_PATH=/hpcfs/home/ciencias_biologicas/af.lizcano/GeneMark-ETP/bin/gmes/ProtHint/bin/

braker.pl --gff3 --threads 8 --workingdir=$output_name --genome=$assembly_path --rnaseq_sets_ids=ERR11837462,SRR12527964,SRR28005379,SRR6148376,SRR6148372,SRR6148369,SRR6148366 --rnaseq_sets_dir=$rnaseq_dir --prot_seq=/hpcfs/home/ciencias_biologicas/af.lizcano/Gencore/Proyecto_Bombus/hymenoptera_proteins.fa --species=$species_name

conda deactivate
source activate agat

agat_sp_filter_incomplete_gene_coding_models.pl --gff ${output_name}/braker.gff3 --fasta $assembly_path -o temp_file1.gff #remove IDs with incomplete gene models
agat_sp_filter_by_ORF_size.pl --gff temp_file1.gff -o temp_file2.gff #removes genes with ORF length <100 aa
agat_sp_manage_IDs.pl --gff temp_file2.gff --prefix ${species_name}_ -o ${output_name}.gff #change the IDs of the final gff file

rm temp_file*
########################################################################################
