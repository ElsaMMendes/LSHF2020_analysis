#!/bin/bash

module load gtdbtk/
export GTDBTK_DATA_PATH=/shared/bank/gtdbtk/release207_v2/
module load prokka/
METAWRAP=/shared/software/miniconda/pkgs/metawrap-1.2-0/bin/metawrap-scripts

nom_analyse=$1
meta_id=$2

mag_treatment_dir=$nom_analyse/5_MAG_treatment
refinement_dir=$nom_analyse/5_MAG_recovery/BIN_REFINEMENT

echo "***************Taxonomically affiliates the bins with gtdb-tk***************"
module load gtdbtk/
export GTDBTK_DATA_PATH=/shared/bank/gtdbtk/release226/
gtdbtk classify_wf --genome_dir $refinement_dir/metawrap_${meta_id}/metawrap_50_10_bins_filteredBIS -x fa --out_dir $mag_treatment_dir/gtdbtk_v226/${meta_id} --pplacer_cpus $SLURM_CPUS_PER_TASK --cpus $SLURM_CPUS_PER_TASK


echo "annotate bins with prokka, important to obtain rRNAs & tRNAs - otherwise preger"
#for archaea
for bin_name in $(tail -n +2 $mag_treatment_dir/gtdbtk/${meta_id}/classify/gtdbtk.ar53.summary.tsv | cut -f 1 ); do
  ${METAWRAP}/shorten_contig_names.py $refinement_dir/metawrap_${meta_id}/metawrap_50_10_bins/${bin_name}.fa > $mag_treatment_dir/prokka/${meta_id}/tmp_${bin_name}.fa
  prokka --quiet --kingdom Archaea --gcode 11 --cpus $SLURM_CPUS_PER_TASK --outdir $mag_treatment_dir/prokka/${meta_id}/$bin_name --prefix $bin_name $mag_treatment_dir/prokka/${meta_id}/tmp_${bin_name}.fa
  rm $mag_treatment_dir/prokka/${meta_id}/tmp_${bin_name}.fa
done

#for bacteria
for bin_name in $(tail -n +2  $mag_treatment_dir/gtdbtk/${meta_id}/classify/gtdbtk.bac120.summary.tsv | cut -f 1 ); do
  ${METAWRAP}/shorten_contig_names.py $refinement_dir/metawrap_${meta_id}/metawrap_50_10_bins/${bin_name}.fa > $mag_treatment_dir/prokka/${meta_id}/tmp_${bin_name}.fa
  prokka --quiet --gcode 11 --cpus $SLURM_CPUS_PER_TASK --outdir $mag_treatment_dir/prokka/${meta_id}/$bin_name --prefix $bin_name $mag_treatment_dir/prokka/${meta_id}/tmp_${bin_name}.fa
  rm $mag_treatment_dir/prokka/${meta_id}/tmp_${bin_name}.fa
done
