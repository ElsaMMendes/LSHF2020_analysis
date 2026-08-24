#!/usr/bin/env python3
import pandas as pd

#attention!!! changer la fin du script si on veut étudier par keggpathway (ko) et non par catégorie métabolique
import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Traitement de données pour un échantillon donné.")
    parser.add_argument("sample_name", help="Nom de l’échantillon (ex: CAP20-A)")
    parser.add_argument("site_name", help="ex: CAP20")
    parser.add_argument("pathfile_name", help="ex: amino_acid_metabolism.txt")
    return parser.parse_args()

args = parse_args()
sample = args.sample_name
site = args.site_name
pathfile = args.pathfile_name

# Charger les fichiers
ego_nog_file = f"../eggNOG/{site}_eggnog.emapper.annotations"
quant_file = f"../salmon_quant/{site}/{sample}_quantified/quant.sf"
carbo_file = f"metabolism_kegg_path_lists/{pathfile}"

# Lire les KO liés au métabolisme des glucides
with open(carbo_file, "r") as f:
    carbohydrate_kos = set(line.strip() for line in f)

# Lire test_eggNOG.txt en ignorant les 4 premières lignes
ego_nog = pd.read_csv(ego_nog_file, sep="\t", skiprows=4)

# Remplacer #query par query
ego_nog.columns = [col.replace("#query", "query") for col in ego_nog.columns]

# Garder uniquement les colonnes d'intérêt (query, KEGG_Pathway)
ego_nog = ego_nog[['query', 'KEGG_Pathway']].dropna()

# Séparer les KEGG_Pathway en plusieurs lignes
ego_nog = ego_nog.assign(KEGG_Pathway=ego_nog['KEGG_Pathway'].str.split(",")).explode('KEGG_Pathway')

# Filtrer les KO présents dans carbohydrate_metabolism.txt
ego_nog = ego_nog[ego_nog['KEGG_Pathway'].isin(carbohydrate_kos)]

# Charger quant.txt et sélectionner les colonnes nécessaires
quant = pd.read_csv(quant_file, sep="\t")
quant = quant[['Name', 'TPM']]

# Associer les TPM aux gènes trouvés
ego_nog = ego_nog.merge(quant, left_on='query', right_on='Name', how='left').drop(columns=['Name'])

# Éliminer les doublons : un gène ne peut apparaître qu'une seule fois
#ego_nog = ego_nog.drop_duplicates(subset=['query'])

#if fractionnement: 
#ego_nog['n_pathways'] = ego_nog.groupby('query')['KEGG_Pathway'].transform('count')
#ego_nog['TPM_fraction'] = ego_nog['TPM'] / ego_nog['n_pathways']
#result = ego_nog.groupby('KEGG_Pathway')['TPM_fraction'].sum() => a la place de la ligne suivante

# Agréger les TPM par KEGG_Pathway
result = ego_nog.groupby('KEGG_Pathway')['TPM'].sum().reset_index()

# Sauvegarder le résultat
result.to_csv(f"aggregated_TPM_{sample}_{pathfile}", sep="\t", index=False)

print("Traitement terminé. Résultats enregistrés")
