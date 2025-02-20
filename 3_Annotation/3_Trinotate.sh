#!/bin/bash

# ###### Parameters area for resource request to SLURM ############################
#
#SBATCH --job-name=Trinotate	#Job Name.
#SBATCH -p medium			#Partition/queue in which to run the job.
#SBATCH -N 1				#Number of nodes required.
#SBATCH -n 1				#Parallel task.
#SBATCH --mem=128G	#Memory per CPU
#SBATCH --time=124:00:00			#Maximum run time	
#SBATCH -o trinotate.o%j			#Output file name
#
########################################################################################

assembly_path=$1 #path to the assembly file
output_name=$2 #name for output

# ################## Module Loading Area ############################################

module load anaconda
source activate trinotate

########################################################################################

blastx -query ${SPECIE}.transcripts.cdna.fa -db /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate2/uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > ${SPECIE}_cdna_uniprot.blastx 2> ${SPECIE}_cdna_uniprot.blastx.err 

blastp -query ${SPECIE}.proteins.fa -db /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate2/uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > ${SPECIE}_proteins_uniprot.blastp 2> ${SPECIE}_proteins_uniprot.blastp.err 

hmmscan --cpu 8 --noali --domtblout ${SPECIE}_pfam_HMMer.out /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate2/Pfam-A.hmm ${SPECIE}.proteins.fa > ${SPECIE}_pfam_HMMer.log 

/hpcfs/home/ciencias_biologicas/af.lizcano/tmhmm-2.0c/bin/tmhmm --short < ${SPECIE}.proteins.fa > ${SPECIE}_tmhmm.out 2> ${SPECIE}_tmhmm.log 

conda deactivate
source activate signalp

signalp6 --fastafile ${SPECIE}.proteins.fa --output_dir ${SPECIE}_signalp --organism eukarya

########################################################################################
