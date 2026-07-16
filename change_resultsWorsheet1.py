import pandas as pd

# Charger le fichier (adapte le séparateur si besoin)
df = pd.read_csv("/shared/projects/feob_iron_oxidizing_bacteria_/META/ANALYSIS_MG/5_MAG_treatment/RESULTS_METABOLIC-Y320/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet1.tsv", sep="\t")

# 1. Garder uniquement les colonnes "bin.*Hits"
bin_hits_cols = [col for col in df.columns if col.startswith("bin") and col.endswith("Hits")]

# Colonnes fixes à garder (tout sauf celles à supprimer)
cols_to_drop = ["Reaction", "Substrate", "Product", "Hmm detecting threshold"]
base_cols = [col for col in df.columns if col not in cols_to_drop and not col.startswith("bin")]

# Nouveau dataframe avec colonnes voulues
df = df[base_cols + bin_hits_cols]

# 2. Remplacer "None" par "0" dans les colonnes Hits
df[bin_hits_cols] = df[bin_hits_cols].apply(lambda col: col.str.replace("None", "0", regex=False))
#si df[bin_hits_cols] = df[bin_hits_cols].replace("None", "0") => None;None.none ne sera pas remplacer car considéré comme chaine entière 
# 3. Remplacer les espaces par des points dans les headers uniquement
df.columns = df.columns.str.replace(" ", ".", regex=False)

# Sauvegarde
df.to_csv("/shared/projects/feob_iron_oxidizing_bacteria_/META/ANALYSIS_MG/5_MAG_treatment/nb_genomes_per_pathway/matrices/Y320_ArticleGenesPathway_metabolic.txt", sep="\t", index=False)
