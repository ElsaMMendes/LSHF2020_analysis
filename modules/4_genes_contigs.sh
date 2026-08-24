#!/bin/bash

module load eggnog-mapper/
module load barrnap/

nom_analyse=$1
meta_id=$2

assembly_dir=$nom_analyse/2_assembly
geneContig_dir=$nom_analyse/4_genes_contigs
protein_file=${geneContig_dir}/prodigal/proteins_${meta_id}.faa

echo "functional annotatation of proteins"
emapper.py -i ${protein_file} -o $geneContig_dir/eggNOG/${meta_id}_eggnog --output_dir eggnog -m diamond --data_dir /shared/bank/emapperdb/5.0.2/ --cpu 12 --target_orthologs one2one --dbmem --dmnd_ignore_warnings --evalue 0.001 --pident 40 --query_cover 20 --subject_cover 20 --itype proteins --pfam_realign none


module load kofamscan/
#cpu 40 for kofamscan
#to be done 1 time !!!!!!
#wget https://www.genome.jp/ftp/db/kofam/ko_list.gz 
#wget -P /shared/ifbstor1/projects/mesopelagic_carbon_pump/DB/ https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
#tar -zxvf /shared/ifbstor1/projects/mesopelagic_carbon_pump/DB/profiles.tar.gz
#gunzip ko_list.gz
mkdir -p $geneContig_dir/kofamscan
exec_annotation -o $geneContig_dir/kofamscan/${meta_id}_result_kofamscan.txt ${protein_file} --ko-list=/shared/ifbstor1/projects/mesopelagic_carbon_pump/DB/ko_list --profile=/shared/ifbstor1/projects/mesopelagic_carbon_pump/DB/profiles/ --cpu 40 --tmp-dir ${protein_file}_tmp

#comptage du nombre de gènes uniques (k###_) présents dans le fichier (qu’ils apparaissent en colonne 1 ou en colonne 2 )
echo "Le nombre de gènes uniques annotés par kofamscan:"
awk '
{
    # cas 1 : première colonne = "*", le gène est dans $2
    if ($1 == "*" && $2 ~ /^k[0-9]+_/) {
        print $2
    }
    # cas 2 : première colonne = un gène (k141…)
    else if ($1 ~ /^k[0-9]+_/) {
        print $1
    }
}' $geneContig_dir/kofamscan/${meta_id}_result_kofamscan.txt | sort -u | wc -l


grep "*" $geneContig_dir/kofamscan/${meta_id}_result_kofamscan.txt > $geneContig_dir/kofamscan/${meta_id}_bestHit_kofamscan.txt
#attention changer chemin input et outpul file en début de script
./filter_kofam_hits.py ${meta_id}
grep "*" $geneContig_dir/kofamscan/${meta_id}_filtered_kofamscan.txt > $geneContig_dir/kofamscan/${meta_id}_filteredBestHit_kofamscan.txt


echo "Searching for 16S/23S/5S"
mkdir -p $geneContig_dir/other_features/ribosomal
barrnap --outseq $geneContig_dir/other_features/ribosomal/${meta_id}.fasta $assembly_dir/${meta_id}/final_renamed.contigs.fa
#then separates 16S and from others to create ${meta_id}_16S.fasta and ${meta_id}others.fasta

cp $geneContig_dir/prodigal/proteins_${meta_id}.faa $assembly_dir/${meta_id}/

eval "$(conda shell.bash hook)"
conda activate fegenie-1.2
FeGenie.py -bin_dir $assembly_dir/${meta_id}/ -bin_ext faa -out $geneContig_dir/fegenie/${meta_id} -t $SLURM_CPUS_PER_TASK --orfs --meta
rm -f $assembly_dir/${meta_id}/proteins_${meta_id}.faa
conda deactivate
