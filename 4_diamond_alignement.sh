#!/bin/bash

nom_analyse=$1
protein_file=$2
site=$3

geneContig_dir=$nom_analyse/4_genes_contigs

module load diamond/

fmt="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle"
echo "$fmt" > $geneContig_dir/protein_align_nr_diamond/head_${site}.m8

diamond blastp -p 16 -d /shared/bank/nr/current/diamond/nr.dmnd -q ${protein_file} -o $geneContig_dir/protein_align_nr_diamond/${site}_aln_diamond.nohead.m8 -f 6 $fmt

cat $geneContig_dir/protein_align_nr_diamond/head_${site}.m8 $geneContig_dir/protein_align_nr_diamond/${site}_aln_diamond.nohead.m8 > $geneContig_dir/protein_align_nr_diamond/${site}_aln_diamond.m8

rm $geneContig_dir/protein_align_nr_diamond/${site}_aln_diamond.nohead.m8
rm $geneContig_dir/protein_align_nr_diamond/head_${site}.m8

