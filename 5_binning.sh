#!/bin/bash

module load metawrap/
module load metabat2/
module load perl/
module load concoct/

nom_analyse=$1
meta_id=$2

assembly_dir=$nom_analyse/2_assembly
mapping_dir=$nom_analyse/3_mapping
binning_dir=$nom_analyse/5_MAG_recovery/INITIAL_BINNING

########################################################################################################
########################   BINNING WITH 3 TOOLS ########################
########################################################################################################

echo "Binning with concoct"
mkdir -p $binning_dir/concoct_bins/${meta_id}/
cut_up_fasta.py $assembly_dir/${meta_id}/final_renamed.contigs.fa -c 10000 -o 0 --merge_last -b $binning_dir/concoct_bins/${meta_id}/contigs_10K.bed > $binning_dir/concoct_bins/${meta_id}/contigs_10K.fa
concoct_coverage_table.py $binning_dir/concoct_bins/${meta_id}/contigs_10K.bed $mapping_dir/${meta_id}*.bam > $binning_dir/concoct_bins/${meta_id}/coverage_table.tsv

concoct --composition_file $binning_dir/concoct_bins/${meta_id}/contigs_10K.fa --coverage_file $binning_dir/concoct_bins/${meta_id}/coverage_table.tsv --basename $binning_dir/concoct_bins/${meta_id}/bins --threads $SLURM_CPUS_PER_TASK --length_threshold 1000
merge_cutup_clustering.py $binning_dir/concoct_bins/${meta_id}/bins_clustering_gt1000.csv > $binning_dir/concoct_bins/${meta_id}/clustering_merge.csv
mkdir -p $binning_dir/concoct_bins/${meta_id}/bins_concoct
extract_fasta_bins.py $assembly_dir/${meta_id}/final_renamed.contigs.fa $binning_dir/concoct_bins/${meta_id}/clustering_merge.csv --output_path $binning_dir/concoct_bins/${meta_id}/bins_concoct

cd $binning_dir/concoct_bins/${meta_id}/bins_concoct/
for filename in *.fa ; do mv "${filename}" "concoct.${filename}" ; done

echo "Performing binning with metabat2"
mkdir -p $binning_dir/metabat2_bins/${meta_id}
jgi_summarize_bam_contig_depths --outputDepth $binning_dir/metabat2_bins/${meta_id}_depth.txt $mapping_dir/${meta_id}*.bam
metabat2 -i $assembly_dir/${meta_id}/final_renamed.contigs.fa -a $binning_dir/metabat2_bins/${meta_id}_depth.txt -o $binning_dir/metabat2_bins/${meta_id}/metabat2 -m 1500 --numThreads $SLURM_CPUS_PER_TASK

echo "Binning of contigs with maxbin2"
mkdir -p $binning_dir/maxbin2_bins/${meta_id}
cut -f1,3- $binning_dir/metabat2_bins/${meta_id}_depth.txt > $binning_dir/maxbin2_bins/${meta_id}_depth.txt
run_MaxBin.pl -contig $assembly_dir/${meta_id}/final_renamed.contigs.fa -out $binning_dir/maxbin2_bins/${meta_id}/maxbin2 -abund $binning_dir/maxbin2_bins/${meta_id}_depth.txt -min_contig_length 1000 -thread $SLURM_CPUS_PER_TASK
mkdir $binning_dir/maxbin2_bins/${meta_id}/bins_maxbin2
cd $binning_dir/maxbin2_bins/${meta_id}/
for file in maxbin2*.fasta ; do mv "${file}" "bins_maxbin2/${file%.fasta}.fa" ; done
