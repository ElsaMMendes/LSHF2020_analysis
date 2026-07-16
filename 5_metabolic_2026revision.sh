#!/bin/bash

#SBATCH -A feob_iron_oxidizing_bacteria_
#SBATCH --job-name=metabolic
#SBATCH --mem 220GB
#SBATCH --cpus-per-task=50
#SBATCH --partition long
#SBATCH --output=/shared/projects/feob_iron_oxidizing_bacteria_/jobLog_metabolic-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elsa.mendes@mio.osupytheas.fr


path_genomes=$1
#[path_to_folder_with_genome_files_of_a_previous_run] 
output_name=$2
#[output_directory_of_a_previous_run]
txt_location_pe=$3

#export PERL5LIB=$HOME/perl5/lib/perl5:$HOME/perl5/lib/perl5/x86_64-linux-thread-multi:$PERL5LIB
#perl -MStatistics::Descriptive -e 'print "Perl module OK\n"'
#conda create -n metabolic -c bioconda perl perl-statistics-descriptive hmmer diamond prodigal
#conda activate metabolic

#export GTDBTK_DATA_PATH=/shared/bank/gtdbtk/release226/

#module load hmmer/
#module load prodigal/
#module load sambamba/
#module load bamtools/
#module load coverm/
#module load r/3.6.3
#module load diamond/
#module load samtools/
#module load bowtie2/
#module load gtdbtk/
#module load perl/

eval "$(conda shell.bash hook)"
conda activate /shared/projects/feob_iron_oxidizing_bacteria_/scripts/METABOLIC_v4.0_v2
#first run:
#perl ../METABOLIC/METABOLIC-C.pl -t 50 -m-cutoff 0.75 -in-gn $path_genomes -kofam-db full -o $output_name -r $txt_location_pe 
#redoing cause depth files not generated
#perl ../METABOLIC/METABOLIC-C.2nd_run.pl -in-gn $path_genomes -r $txt_location_pe  -o $output_name -2nd-run true -2nd-run-suffix 2nd_run_test
#obtain files for class
perl ../METABOLIC/METABOLIC-C.2nd_run.pl -in-gn $path_genomes -r $txt_location_pe -o $output_name -2nd-run true -2nd-run-suffix MW_score_tax_class -depth-file All_gene_collections_mapped.depth.2nd_run_test.txt -tax class
#ontain files for family
perl ../METABOLIC/METABOLIC-C.2nd_run.pl -in-gn $path_genomes -r $txt_location_pe -o $output_name -2nd-run true -2nd-run-suffix MW_score_tax_family -depth-file All_gene_collections_mapped.depth.2nd_run_test.txt -tax family
#-tax default: "phylum"; other options: "class", "order", "family", "genus", "species", and "bin" (MAG itself)).

conda deactivate