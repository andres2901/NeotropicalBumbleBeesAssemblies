# ###### Parameters area for resource request to SLURM ############################
#
#SBATCH --job-name=Tree	#Job Name.
#SBATCH -p medium			#Partition/queue in which to run the job.
#SBATCH -N 1				#Number of nodes required.
#SBATCH -n 1				#Parallel task.
#SBATCH --mem=256G	#Memory per CPU
#SBATCH --time=72:00:00			#Maximum run time	
#SBATCH -o tree.o%j			#Output file name
#
########################################################################################

genes_path=$1 #path to the directory where multifasta of genes are stored

# ################## Module Loading Area ############################################

module load anaconda
source activate MAFFT

########################################################################################

#Alignment
ls ${genes_path}/genes* | awk -F '.' '{print $1}' | while read line #the multifasta files should start whit "genes" and follow by the name of the gene
do
mafft --thread -1 --genafpair --maxiterate 1000 ${line}.fa > aligned_${line}.fa
done

#Trimmed alignment

conda deactivate

ls aligned* | while read line
do
/hpcfs/home/ciencias_biologicas/af.lizcano/Gblocks_0.91b/Gblocks $line -t=d -b5=h -e=.fa
done

rm *.htm
mkdir Alignments
mv *.fa.fa Alignments

#Tree
source activate Iqtree
iqtree2 -p Alignments/ --prefix concat -B 1000 -merit AICc

mkdir Trees
mv concat* Trees/
