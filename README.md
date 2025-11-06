# Neotropical Bumble bees assembly

Bioinformatics analysis pipeline till annotation for:

>Lizcano-Salas et al. (2025). First genome assemblies of Neotropical *Thoracobombus* Bumble Bees *Bombus pauloensis* and *Bombus pullatus*. bioRxiv. https://doi.org/10.1101/2025.10.13.682240

The workflow was written for and run on Universidad de los Andes HPC which uses SLURM worload manager. The bash scripts should be modified based on the specified system were want to be run and the name and path of each file.

## 0 Basecalling

Initially, Basecalling was performed on the GridION machine using Dorado v0.7.3, followed by SAMtools conversion. The commands used were:

```
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.0.0 <pod5_dir> > <ouput_name>.bam #basecalling
samtools fastq <input_name>.bam > <output_name>.fastq #convert format
```

## 1 Assembly

<div align="center">
	<img src="https://github.com/andres2901/NeotropicalBumbleBeesAssemblies/blob/main/0_images/assembly.svg" alt="Assembly workflow" width="500" height="700">
</div>

Using the basecalled files generated for each species, we execute the scripts within the `1_assembly` directory in the following order:

* `1_Initial_Quality.sh`: Determines the quality of the raw reads.
* `2_Filtering.sh`: Removes adapters and barcodes, filters reads based on quality and length, and generates the two datasets described in the paper.
* `3_Genome_Size.sh`: Estimates genome size based on read k-mer analysis.
* `4_Assembly.sh`: Assembles the genome.
* `5_Correction.sh`: Performs polishing and verifies the statistics of the initial draft assembly.
* `6A_Coverage.sh`: Determines coverage and depth.
* `6B_Blast.sh`: Performs BLAST searches on the contigs.
* `6C_diamond.sh`: Performs DIAMOND searches on the contigs.
* `7_Artifact_Check.sh`: Runs Compleasm for manual artifact checking.
* `8_Finishing_Assembly.sh`: Performs gap filling and subsequent polishing of the refined contigs.

Scripts beginning with '6' were used for BlobToolKit contamination checks. The BlobToolKit dataset was manually generated using the following commands:

```
blobtools create --fasta <assembly>.fasta --taxdump /hpcfs/home/ciencias_biologicas/af.lizcano/taxdump <dataset_directory_name> #create dataset
blobtools add --hits <blast_ouput> --hits <diamon_output> --taxrule bestsumorder --taxdump ~/blobtoolkit/taxdump <dataset_directory_name> #add blast and diamond hits
blobtools add --cov <coverage>.bam <dataset_directory_name> #add coverage data
```

Following the BlobToolKit check, a second manual contamination screening was performed using the Foreign Contamination Screen (FCS, https://github.com/ncbi/fcs) tool suite, specifically the FCS-GX function on the Galaxy platform (https://usegalaxy.org/).

Artifact checking was conducted manually, utilizing the results of the 7_Artifact_Check.sh script and dot plot analysis with Gepard (https://github.com/univieCUBE/gepard).

## 2 Tree reconstruction

<div align="center">
	<img src="https://github.com/andres2901/NeotropicalBumbleBeesAssemblies/blob/main/0_images/tree_reconstruction.svg" alt="Tree reconstruction workflow" width="500" height="700">
</div>

Following the final assembly from the previous section, a phylogenetic analysis was conducted to verify the origin of the assemblies. Genes were initially retrieved manually using the following commands:

```
makeblastdb -in <assembly_file> -dbtype nucl -out <database_name> #create database for blast
blastn -db <database_name> -query <genes_file> -max_target_seqs 1 -out <ouput_name> #run blast
samtools faidx <assembly> <region_of_interest> > <output_name> #extract region of interest, this is performed individually per region
revseq <input_path> <ouput_path> #create the reverse complement when necessary
```

Following the gene extraction, multifasta files were manually created, including sequences from other *Thoracobombus* species and an outgroup. The `Tree_reconstruction.sh` script, located in the `2_Tree_reconstruction` directory, was then executed using these multifasta files.

## 3 Annotation

<div align="center">
	<img src="https://github.com/andres2901/NeotropicalBumbleBeesAssemblies/blob/main/0_images/annotation.svg" alt="Annotation workflow" width="500" height="700">
</div>

Using the assemblies generated in the [first section](#1-Assembly), genome annotation was performed by executing the scripts within the `3_Annotation` directory in the following order:

* 0_masking.sh: Identifies transposable elements and masks the assembly.
* 1_Braker.sh: Predicts genes.
* 2_Trinotate.sh: Annotates the function of predicted genes.
* 3_non-coding_annotation.sh: Annotates non-coding RNAs.

Transposable element identification using LTR_Pred for `TransposonUltimate`pipeline was performed manually, prior to running the `0_masking.sh` script. The command used was:

```
source activate udocker_prod #activate conda environment with LTR_Pred
udocker run ltrpred #run the udocker container
mkdir ltrpred_data

scp <source username@IP>:<path_to_assembly> ltrpred_data/sequence.fa

R
LTRpred:: LTRpred(genome.file = "ltrpred_data/sequence.fa", cores=16) #run ltrpred
quit()

scp -r sequence_ltrpred <destination username@IP>:<path_to_directory>
```

Some results, with contig names updated for NCBI submission, are stored in the `3_Annotation` folder.
