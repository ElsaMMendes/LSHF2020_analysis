#!/bin/bash

module load bowtie2/
module load samtools/

nom_analyse_reads=$1
nom_analyse_index=$2
meta_id=$3
folder_results_path=$4

cleaning_dir=$nom_analyse_reads/1_cleaning
replicate_file=$nom_analyse_reads/replicates_file.txt
index_path_prefix=$nom_analyse_index/2_assembly/assembly_index/${meta_id}_index
geneContig_dir=$nom_analyse_index/4_genes_contigs
salmon_path=$geneContig_dir/salmon_quant

echo "aligns reads to the contigs for each sample/library"
for rep in $(cat $replicate_file); do
  bowtie2 -1 $cleaning_dir/${meta_id}-${rep}/${meta_id}-${rep}_P1.fastq.gz -2 $cleaning_dir/${meta_id}-${rep}/${meta_id}-${rep}_P2.fastq.gz -x $index_path_prefix  | samtools view -Sb | samtools sort > $folder_results_path/${meta_id}-${rep}_mapped.bam
  samtools index $folder_results_path/${meta_id}-${rep}_mapped.bam
done

echo "contig coverage for each site  (replicates combined)"
samtools coverage $folder_results_path/${meta_id}*.bam -o $folder_results_path/samtools/${meta_id}_contigs_coverage.txt

module load salmon/
#INDEX
salmon index -t $geneContig_dir/prodigal/nucleotides_${meta_id}.fa -i $salmon_path/${meta_id}/${meta_id}_index -k 31

for rep in $(cat $replicate_file); do
  salmon quant -i $salmon_path/${meta_id}/${meta_id}_index -l IU -1 $cleaning_dir/trimmomatic_QC/${meta_id}-${rep}/${meta_id}-${rep}_P1.fastq.gz -2 $cleaning_dir/trimmomatic_QC/${meta_id}-${rep}/${meta_id}-${rep}_P2.fastq.gz --validateMappings -o $salmon_path/${meta_id}/${meta_id}/${meta_id}-${rep}_quantified
done


