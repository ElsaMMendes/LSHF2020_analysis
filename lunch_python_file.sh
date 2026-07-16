#!/bin/bash

#SBATCH -A feob_iron_oxidizing_bacteria_
#SBATCH --job-name=paths
#SBATCH --partition fast
#SBATCH --mem 20GB
#SBATCH --cpus-per-task=10
#SBATCH --output=/shared/projects/feob_iron_oxidizing_bacteria_//META/ANALYSIS_MG/4_genes_contigs/jobLog_pathway_article-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elsa.mendes@mio.osupytheas.fr

for file in metabolism_kegg_path_lists/*; do
    filename=$(basename "$file")
    ./keggPathway_abundanceV2-pathway.py CAP20-A CAP20 "$filename"
    ./keggPathway_abundanceV2-pathway.py CAP20-B CAP20 "$filename"
    ./keggPathway_abundanceV2-pathway.py CAP20-C CAP20 "$filename"

    ./keggPathway_abundanceV2-pathway.py NTE20-A NTE20 "$filename"
    ./keggPathway_abundanceV2-pathway.py NTE20-B NTE20 "$filename"
    ./keggPathway_abundanceV2-pathway.py NTE20-C NTE20 "$filename"

    ./keggPathway_abundanceV2-pathway.py LL20-A LL20 "$filename"
    ./keggPathway_abundanceV2-pathway.py LL20-B LL20 "$filename"
    ./keggPathway_abundanceV2-pathway.py LL20-C LL20 "$filename"

    ./keggPathway_abundanceV2-pathway.py Y320-A Y320 "$filename"
    ./keggPathway_abundanceV2-pathway.py Y320-B Y320 "$filename"
    ./keggPathway_abundanceV2-pathway.py Y320-C Y320 "$filename"
done


samples="CAP20-A CAP20-B CAP20-C NTE20-A NTE20-B NTE20-C LL20-A LL20-B LL20-C Y320-A Y320-B Y320-C"

# Header
#echo -ne "Category"
#for sample in $samples; do
    #echo -ne "\t$sample"
#done
#echo "" > matrix.tsv

# Pour chaque catégorie
#for file in metabolism_kegg_path_lists/*.txt; do
    
    #filename=$(basename "$file")
    #category=${filename%.txt}
    
    #echo -ne "$category" >> matrix.tsv
    
    #for sample in $samples; do
        
        #infile="aggregated_TPM_${sample}_${filename}"
        
        #sum=$(awk -F'\t' '
            #NR>1 {s+=$2}
            #END{print s+0}
        #' "$infile")
        
        #echo -ne "\t$sum" >> matrix.tsv
    #done
    
    #echo "" >> matrix.tsvl

#done




echo -ne "KEGG_Pathway" > matrix_energy_metabolism.tsv
for sample in $samples; do
    echo -ne "\t$sample" >> matrix_energy_metabolism.tsv
done
echo "" >> matrix_energy_metabolism.tsv

# Récupérer la liste des 9 KEGG pathways (à partir d’un des fichiers)
first_file=$(ls *energy_metabolism.txt | head -n1)

awk -F'\t' 'NR>1 {print $1}' "$first_file" | while read pathway; do
    
    echo -ne "$pathway" >> matrix_energy_metabolism.tsv
    
    for sample in $samples; do
        
        infile="aggregated_TPM_${sample}_energy_metabolism.txt"
        
        # Extraire la valeur TPM correspondant au pathway
        value=$(awk -F'\t' -v p="$pathway" '
            $1==p {print $2}
        ' "$infile")
        
        # Si pathway absent → mettre 0
        value=${value:-0}
        
        echo -ne "\t$value" >> matrix_energy_metabolism.tsv
    done
    
    echo "" >> matrix_energy_metabolism.tsv

done
