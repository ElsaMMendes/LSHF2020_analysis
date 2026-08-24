#!/bin/bash

nom_analyse=$1
contigs_fasta=$2
database_folder=$3
taxonomy_folder=$4
predicted_proteins_fasta=$5
alignment_file=$6
site=$7

geneContig_dir=$nom_analyse/4_genes_contigs

module load cat/

cat ${alignment_file} | cut -f 1-12 > $geneContig_dir/taxonomy_contigs/${site}_cut_diamond.m8
tail -n +2 $geneContig_dir/taxonomy_contigs/${site}_cut_diamond.m8 > $geneContig_dir/taxonomy_contigs/${site}_cut2_diamond.m8

CAT contigs -c ${contigs_fasta} -d ${database_folder} -t ${taxonomy_folder} -p ${predicted_proteins_fasta} --out_prefix $geneContig_dir/taxonomy_contigs/${site}_CAT -a $geneContig_dir/taxonomy_contigs/${site}_cut2_diamond.m8

CAT add_names -i $geneContig_dir/taxonomy_contigs/${site}_CAT.contig2classification.txt -o $geneContig_dir/taxonomy_contigs/${site}_CAT.contig2classification.official_names.txt -t ${taxonomy_folder}  --exclude_scores --only_official

CAT summarise -c ${contigs_fasta} -i $geneContig_dir/taxonomy_contigs/${site}_CAT.contig2classification.official_names.txt -o $geneContig_dir/taxonomy_contigs/${site}_CAT.summary.txt

rm $geneContig_dir/taxonomy_contigs/${site}_cut2_diamond.m8
rm $geneContig_dir/taxonomy_contigs/${site}_cut_diamond.m8
