#!/bin/bash

module load metawrap/
module load perl/

nom_analyse=$1
meta_id=$2
#parameters for dereplication:
comp=$3
contam=$4
refinement_dir=$nom_analyse/5_MAG_treatment/BIN_REFINEMENT

mapping_dir=$nom_analyse/3_mapping
binning_dir=$nom_analyse/5_MAG_recovery/INITIAL_BINNING
refinement_dir=$nom_analyse/5_MAG_recovery/BIN_REFINEMENT
mag_treatment_dir=$nom_analyse/5_MAG_treatment

echo "Refine bin sets and access quality of bins with checkm"
metawrap bin_refinement -o $refinement_dir/metawrap_${meta_id} -t $SLURM_CPUS_PER_TASK -A $binning_dir/metabat2_bins/${meta_id}/ -B $binning_dir/maxbin2_bins/${meta_id}/bins_maxbin2/ -C $binning_dir/concoct_bins/${meta_id}/bins_concoct/ -c 50 -x 10

module load coverm/
echo "calculing coverage of genomes"
coverm genome --bam-files $mapping_dir/${meta_id}*.bam --genome-fasta-directory $refinement_dir/metawrap_${meta_id}/metawrap_50_10_bins/ -x fa -o $mag_treatment_dir/coverM/${meta_id}_abundance -t 8  --min-covered-fraction 0

####################### REGROUPER les bins de tous les sites ensembles pour par exemple ensuite faire une derep ou autre
mkdir -p $refinement_dir/ALL_50_10_bins/
cp $refinement_dir/metawrap_${meta_id}/metawrap_50_10_bins/bin.* $refinement_dir/ALL_50_10_bins/
cd $refinement_dir/ALL_50_10_bins/
for filename in bin.*.fa ; do mv "${filename}" "${meta_id}_${filename}" ; done

################## dereplication

module load drep/
mkdir -p $refinement_dir/drep/drep_${meta_id}

#intra dreplication; fait aussi pour 95% identitité ( -sa 0.95 ) 
echo "genome,completeness,contamination" > $refinement_dir/drep/drep_${meta_id}/bins_stats_${meta_id}.csv
cut -f1,2,3 $refinement_dir/${meta_id}/metawrap_50_10_bins.stats | awk -F "\t" '{print $1".fa,"$2","$3}' >> $refinement_dir/drep/drep_${meta_id}/bins_stats_${meta_id}.csv
sed -i '2d' $refinement_dir/drep/drep_${meta_id}/bins_stats_${meta_id}.csv
dRep dereplicate $refinement_dir/drep/drep_${meta_id} -p $SLURM_CPUS_PER_TASK --S_algorithm fastANI -g $refinement_dir/${meta_id}/metawrap_50_10_bins/* --genomeInfo $refinement_dir/drep/drep_${meta_id}/bins_stats_${meta_id}.csv --completeness $comp --contamination $contam -pa 0.9 -sa 0.99

#for inter dereplication 
cp $refinement_dir/drep/drep_${meta_id}/bins_stats_${meta_id}.csv $refinement_dir/drep/drep_${meta_id}/bins_stats_${meta_id}_BIS.csv
sed -i '1d' $refinement_dir/drep/drep_${meta_id}/bins_stats_${meta_id}_BIS.csv
sed -i "s/bin/${meta_id}_bin/g" $refinement_dir/drep/drep_${meta_id}/bins_stats_${meta_id}_BIS.csv

#TO RUN ALONE AT THE END ; run inter duplication :
#echo "genome,completeness,contamination" > $refinement_dir/drep/drep_ALL/bins_stats_ALL.csv
#for site in $(cat $sites_file); do cat $refinement_dir/drep/drep_${site}/bins_stats_${site}_BIS.csv >> $refinement_dir/drep/drep_ALL/bins_stats_ALL.csv ; done
#dRep dereplicate $refinement_dir/drep/drep_ALL -p $SLURM_CPUS_PER_TASK --S_algorithm fastANI -g $refinement_dir/ALL_metawrap_50_10_bins/* --genomeInfo $refinement_dir/drep/drep_ALL/bins_stats_ALL.csv --completeness $comp --contamination $contam -pa 0.9 -sa 0.99


echo "***************QUALITY OF MAGS***************"
module load checkm2/
checkm2 predict --input $refinement_dir/metawrap_${meta_id}/metawrap_50_10_bins --output-directory ${refinement_dir}/checkM2_metawrap_${meta_id} -t 16 --extension fa --database_path /shared/bank/checkM2/version_3/CheckM2_database/uniref100.KO.1.dmnd 
#avec checkm: checkm lineage_wf -x fa META/ANALYSIS_MG/5_MAG_recovery/BIN_REFINEMENT/${site}_bins META/ANALYSIS_MG/5_MAG_recovery/BIN_REFINEMENT/checkM_ALL/checkM_${site} -t 12 --reduced_tree
#--general par défaut (--specific pour maximiser HQ bactériens mais pas trop si archées présentes)

awk -F'\t' '
NR==1 {print; next}
{
  if ($2 >= 90 && $3 <= 5) class="HQ";
  else if ($2 >= 50 && $3 <= 10) class="MQ";
  else if ($2 < 50 && $3 <= 10) class="LQ";
  else class="REJECT";
  print $0 "\t"class
}' ${refinement_dir}/checkM2_metawrap_${meta_id}/quality_report.tsv > ${refinement_dir}/checkM2_metawrap_${meta_id}/quality_report_classified.tsv

mkdir -p $refinement_dir/metawrap_${meta_id}/metawrap_50_10_bins_filtered

awk -F'\t' 'NR>1 && !($NF=="REJECT" || $NF=="LQ") {print $1}' ${refinement_dir}/checkM2_metawrap_${meta_id}/quality_report_classified.tsv > ${refinement_dir}/checkM2_metawrap_${meta_id}/mags_to_keep.txt

while read mag; do
    cp $refinement_dir/metawrap_${meta_id}/metawrap_50_10_bins/${mag}.fa \
       $refinement_dir/metawrap_${meta_id}/metawrap_50_10_bins_filtered/
done < ${refinement_dir}/checkM2_metawrap_${meta_id}/mags_to_keep.txt
