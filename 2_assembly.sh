#!/bin/bash

help_message () {
	echo ""
	echo "Usage: sbatch 2_assembly.sh <analysis_folder> <assembly_name> <forward_files> <reverse_files> "
	echo "Options:"
	echo ""
	echo "	<analysis_folder>      folder of study"
	echo "	<assembly_name>         name to give to assembly"
    echo " <forward_files>  path of forward files to assemble (list of several files separated by comma)"
	echo "	<reverse_files>		reverse files to assemble (list of several files separated by comma)"
	echo "";}


module load quast/
module load bowtie2/
module load megahit/
module load prodigal/
module load blast/

#arguments
nom_analyse=$1
meta_id=$2
reads_1=$3
reads_2=$4

assembly_dir=$nom_analyse/2_assembly
geneContig_dir=$nom_analyse/4_genes_contigs
metaquast_dir=$assembly_dir/metaquast
METAWRAP=/shared/software/miniconda/pkgs/metawrap-1.2-0/bin/metawrap-scripts

###############


echo "Running assembly"
megahit -1 ${reads_1} -2 ${reads_2} -o $assembly_dir/${meta_id} -t $SLURM_CPUS_PER_TASK

#The contigs are sorted by length and their naming is changed to resemble that of SPAdes (including the contig ID, length, and coverage)
echo "Changing names of megahit to better format"
${METAWRAP}/fix_megahit_contig_naming.py $assembly_dir/${meta_id}/final.contigs.fa 0 > $assembly_dir/${meta_id}/final_renamed.contigs.fa

echo "Assessing the quality of assembly with QUAST"
mkdir -p $metaquast_dir/${meta_id}
metaquast.py --max-ref-number 0 --min-contig 0 --memory-efficient -o $metaquast_dir/${meta_id} --labels ${meta_id} $assembly_dir/${meta_id}/final_renamed.contigs.fa

echo "Indexing contigs file for mapping"
bowtie2-build $assembly_dir/${meta_id}/final_renamed.contigs.fa $assembly_dir/assembly_index/${meta_id}_index

echo "Doing protein-coding gene prediction with (meta)prodigal FOR PROKARYOTES"
prodigal -i $assembly_dir/${meta_id}/final_renamed.contigs.fa -o $geneContig_dir/prodigal/genes_${meta_id}.gff -f gff -p meta -a $geneContig_dir/prodigal/proteins_${meta_id}.faa -d $geneContig_dir/prodigal/nucleotides_${meta_id}.fa

echo "Generating a blast DB with proteins in contigs"
makeblastdb -in $geneContig_dir/prodigal/proteins_${meta_id}.faa -parse_seqids -dbtype prot -out $geneContig_dir/blastDB/blastdb${meta_id}

#LANCER MULTIQC APRÈS AVOIR OBTENU TOUS LES ASSEMBLAGES