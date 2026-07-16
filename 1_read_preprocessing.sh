#!/bin/bash

help_message () {
	echo ""
	echo "Usage: sbatch 1_read_preprocessing.sh [options] -o output_dir -1 raw_dir -a adapter.fa "
	echo "Options:"
	echo ""
	echo "	-o STR          output directory (analysis folder)"
	echo "	-1 STR          raw data directory"
	echo "	-a STR		adapter file (and path)"
	echo "";}


########################################################################################################
###############    LOADING IN THE PARAMETERS  AND MAKING SURE EVERYTHING IS SET UP     #################
########################################################################################################

#output=$1

# default params
nom_analyse="false"; seq_data="false"; adapter="false"

# loop through input params
while true; do
	case "$1" in
		-o) nom_analyse=$2; shift 2;;
		-1) seq_data=$2; shift 2;;
		-a) adapter=$2; shift 2;;
		-h | --help) help_message; exit 1; shift 1;;
		--) help_message; exit 1; shift; break ;;
		*) break;;
	esac
done

# check if all parameters are entered, OU logique
if [ "$nom_analyse" = "false" ] || [ "$seq_data" = "false" ] || [ "$adapter" = "false" ]; then
	help_message; exit 1
fi

########################################################################################################
################                    READ CLEANING AND REPORTS                   ########################
########################################################################################################

module load fastqc/
module load trimmomatic/
module load multiqc/

###########
mkdir -p $nom_analyse # S'IL N'EXISTE PAS SINON PAS CREER
cleaning_dir=$nom_analyse/1_cleaning
multiqc_folder=$cleaning_dir/multiqc
mkdir -p $cleaning_dir $cleaning_dir/trimmomatic_QC #si n'existe pas

ls $seq_data/*_R[1].fastq.gz|sed -e 's#.*/##' -e 's#_R[1].fastq.gz##' >> $nom_analyse/samples_file.txt
samples_file=$nom_analyse/samples_file.txt

ls $seq_data/*.gz |sed -e 's#.*/##' | cut -d '-' -f1 | uniq >> $nom_analyse/sites_file.txt
if grep -q "fastq.gz" $nom_analyse/sites_file.txt; then
    rm sites_file.txt
else
    site=$nom_analyse/sites_file.txt
fi

ls $seq_data/*_R[1].fastq.gz|sed -e 's#.*/*-##' -e 's#_R[1].fastq.gz##' | sort | uniq >> $nom_analyse/replicates_file.txt
if grep -q "/" $nom_analyse/replicates_file.txt; then
    rm $nom_analyse/replicates_file.txt
else
    replicate=$nom_analyse/replicates_file.txt
fi
#for F in analyse_test/7_metabat_bin/bins_dir_fasta/*.fa ; do mv "${F}" "${F%.fa}.fasta" ; done

if [ ! -s $samples_file ]; then error "$samples_file file does not exist. Exiting..."; fi

#################

echo "DOING QUALITY CONTROL BEFORE CLEANING"
for sample in $(cat $samples_file); do
  mkdir -p $cleaning_dir/trimmomatic_QC/$sample $cleaning_dir/trimmomatic_QC/$sample/QC-before
  fastqc "$seq_data"/"$sample"_R1.fastq.gz "$seq_data"/"$sample"_R2.fastq.gz  -o $cleaning_dir/trimmomatic_QC/$sample/QC-before
done
wait
multiqc $cleaning_dir/trimmomatic_QC/*/QC-before -o $multiqc_folder/QC_before


echo "DOING THE CLEANING"
for sample in $(cat $samples_file); do
  trimmomatic PE -threads $SLURM_CPUS_PER_TASK -phred33 "$seq_data"/"$sample"_R1.fastq.gz "$seq_data"/"$sample"_R2.fastq.gz \
   "$cleaning_dir"/trimmomatic_QC/"$sample"/"$sample"_P1.fastq.gz "$cleaning_dir"/trimmomatic_QC/"$sample"/"$sample"_U1.fastq.gz \
   "$cleaning_dir"/trimmomatic_QC/"$sample"/"$sample"_P2.fastq.gz "$cleaning_dir"/trimmomatic_QC/"$sample"/"$sample"_U2.fastq.gz \
   ILLUMINACLIP:"$adapter":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:20
done
wait

echo "DOING QC AFTER CLEANING"
for sample in $(cat $samples_file); do
  mkdir -p $cleaning_dir/trimmomatic_QC/$sample/QC-after
  fastqc "$cleaning_dir"/"$sample"/"$sample"_P1.fastq.gz "$cleaning_dir"/"$sample"/"$sample"_P2.fastq.gz  -o "$cleaning_dir"/trimmomatic_QC/"$sample"/QC-after
done
wait
multiqc $cleaning_dir/trimmomatic_QC/*/QC-after -o $multiqc_folder/QC_after

#################

echo "READ QC PIPELINE COMPLETE!!!"