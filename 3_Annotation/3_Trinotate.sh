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

assembly_file=$1 #path of assembly file
gff_file=$2 #path of gff file
PREFIX=$3 #prefix of files for annotation

# ################## Module Loading Area ############################################

module load anaconda
source activate trinotate

########################################################################################

#generate cdna and protein file

/hpcfs/home/ciencias_biologicas/af.lizcano/Trinotate-Trinotate-v4.0.2/util/Trinotate_GTF_or_GFF3_annot_prep.pl --annot $gff_file --genome_fa $assembly_file --out_prefix ${PREFIX}

#run necessary tools for annotation

blastx -query ${PREFIX}.transcripts.cdna.fa -db /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate2/uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > ${PREFIX}_cdna_uniprot.blastx 2> ${PREFIX}_cdna_uniprot.blastx.err 

blastp -query ${PREFIX}.proteins.fa -db /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate2/uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > ${PREFIX}_proteins_uniprot.blastp 2> ${PREFIX}_proteins_uniprot.blastp.err 

hmmscan --cpu 8 --noali --domtblout ${PREFIX}_pfam_HMMer.out /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate2/Pfam-A.hmm ${PREFIX}.proteins.fa > ${PREFIX}_pfam_HMMer.log 

/hpcfs/home/ciencias_biologicas/af.lizcano/tmhmm-2.0c/bin/tmhmm --short < ${PREFIX}.proteins.fa > ${PREFIX}_tmhmm.out 2> ${PREFIX}_tmhmm.log 

conda deactivate
source activate signalp

signalp6 --fastafile ${PREFIX}.proteins.fa --output_dir ${PREFIX}_signalp --organism eukarya

conda deactivate
source activate trinotate

#run eggnoggmapper

Trinotate --db /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate2/TrinotateBoilerplate.sqlite --init --gene_trans_map ${PREFIX}_gene-to-trans-map --transcript_fasta ${PREFIX}.transcripts.cdna.fa --transdecoder_pep ${PREFIX}.proteins.fa

#load results

Trinotate --db /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate2/TrinotateBoilerplate.sqlite --CPU 8 --transcript_fasta ${PREFIX}.transcripts.cdna.fa --transdecoder_pep ${PREFIX}.proteins.fa --run "EggnogMapper"
Trinotate --db /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate2/TrinotateBoilerplate.sqlite --LOAD_swissprot_blastp ${PREFIX}_proteins_uniprot.blastp
Trinotate --db /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate2/TrinotateBoilerplate.sqlite --LOAD_pfam ${PREFIX}_pfam_HMMer.out
Trinotate --db /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate2/TrinotateBoilerplate.sqlite --LOAD_signalp ${PREFIX}_signalp/output.gff3
Trinotate --db /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate2/TrinotateBoilerplate.sqlite --LOAD_EggnogMapper <file>
Trinotate --db /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate2/TrinotateBoilerplate.sqlite --LOAD_tmhmmv2 ${PREFIX}_tmhmm.out
Trinotate --db /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate2/TrinotateBoilerplate.sqlite --LOAD_swissprot_blastx ${PREFIX}_cdna_uniprot.blastx

#generate reports

Trinotate --db /hpcfs/home/ciencias_biologicas/af.lizcano/database/Trinotate2/TrinotateBoilerplate.sqlite --report > ${PREFIX}_trinotateReport.txt
/hpcfs/home/ciencias_biologicas/af.lizcano/Trinotate-Trinotate-v4.0.2/util/report_summary/trinotate_report_summary.pl ${PREFIX}_trinotateReport.txt ${PREFIX}_trinotateReport
########################################################################################
