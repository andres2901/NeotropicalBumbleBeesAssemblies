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

PREFIX=$1 #prefix of files for annotation
output_name=$2 #name for output

# ################## Module Loading Area ############################################

module load anaconda
source activate trinotate

########################################################################################

blastx -query ${PREFIX}.transcripts.cdna.fa -db /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate2/uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > ${PREFIX}_cdna_uniprot.blastx 2> ${PREFIX}_cdna_uniprot.blastx.err 

blastp -query ${PREFIX}.proteins.fa -db /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate2/uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > ${PREFIX}_proteins_uniprot.blastp 2> ${PREFIX}_proteins_uniprot.blastp.err 

hmmscan --cpu 8 --noali --domtblout ${PREFIX}_pfam_HMMer.out /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate2/Pfam-A.hmm ${PREFIX}.proteins.fa > ${PREFIX}_pfam_HMMer.log 

/hpcfs/home/ciencias_biologicas/af.lizcano/tmhmm-2.0c/bin/tmhmm --short < ${PREFIX}.proteins.fa > ${PREFIX}_tmhmm.out 2> ${PREFIX}_tmhmm.log 

conda deactivate
source activate signalp

signalp6 --fastafile ${PREFIX}.proteins.fa --output_dir ${PREFIX}_signalp --organism eukarya

/hpcfs/home/ciencias_biologicas/af.lizcano/Trinotate-Trinotate-v4.0.2/Trinotate /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate2/TrinotateBoilerplate.sqlite init --gene_trans_map ${PREFIX}_geneTranscriptMap.txt --transcript_fasta ${PREFIX}_cdna.fa --transdecoder_pep ${PREFIX}_ProteinsTrinotate.fa
/hpcfs/home/ciencias_biologicas/af.lizcano/Trinotate-Trinotate-v4.0.2/Trinotate /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate/Trinotate.sqlite LOAD_swissprot_blastp ${PREFIX}_proteins_uniprot.blastp
/hpcfs/home/ciencias_biologicas/af.lizcano/Trinotate-Trinotate-v4.0.2/Trinotate /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate/Trinotate.sqlite LOAD_swissprot_blastx ${PREFIX}_cdna_uniprot.blastx
/hpcfs/home/ciencias_biologicas/af.lizcano/Trinotate-Trinotate-v4.0.2/Trinotate /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate/Trinotate.sqlite LOAD_pfam ${PREFIX}_pfam_HMMer.out
/hpcfs/home/ciencias_biologicas/af.lizcano/Trinotate-Trinotate-v4.0.2/Trinotate /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate/Trinotate.sqlite LOAD_tmhmm ${PREFIX}_tmhmm.out
/hpcfs/home/ciencias_biologicas/af.lizcano/Trinotate-Trinotate-v4.0.2/Trinotate /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate/Trinotate.sqlite report > ${PREFIX}_trinotateReport.txt
/hpcfs/home/ciencias_biologicas/af.lizcano/Trinotate-Trinotate-v4.0.2/util/report_summary/trinotate_report_summary.pl ${PREFIX}_trinotateReport.txt ${PREFIX}_trinotateReport

########################################################################################
