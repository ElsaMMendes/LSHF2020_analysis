import os
import pandas as pd
import random
import re


# === 1 ligne par HMM ===
#pas de doublons dans la colonne des profils HMM (sinon prob recherche redondance)
def normalize_matrix(input_file, matrix_file):
    dfInit = pd.read_csv(input_file, sep="\t", dtype=str)
    #garder dbtype str?
    
    # Identifier les colonnes MAG (.Hits)
    hits_cols = [col for col in dfInit.columns if col.endswith(".Hits")]
    
    expanded_rows = []

    # Parcours de chaque ligne avec index d'origine
    for idx, row in dfInit.iterrows():
        #for _, row in df.iterrows():, ignores l’index et tu n’utilises que row
        hmm_files = [h.strip() for h in str(row["Hmm.file"]).split(",")]
        #['K00823.hmm', 'K07250.hmm', 'K13524.hmm', 'K14268.hmm', 'K03918.hmm']
        n_hmms = len(hmm_files)

        # Vérifier la cohérence des longueurs dans toutes les colonnes *.Hits
        for col in hits_cols:
            col_hits = str(row[col]).split(";")
            #['0', '0', '0', '0', 'k141_429300_length_7756_cov_7.0000_1']
            if len(col_hits) != n_hmms:
                raise ValueError(f"Incohérence dans {col} à la ligne {idx} : {col_hits} vs HMMs {hmm_files}")

        # Génération de n lignes pour chaque HMM
        for i in range(n_hmms):
            new_row = row.copy()
            #on recopie le tableau et change/ajoute certaines
            new_row["Hmm.file"] = hmm_files[i]
            new_row["original_index"] = idx  # Ajout de l'index d'origine
            for col in hits_cols:
                #new_row[col] = row[col].split(";")[i].strip()
                new_row[col] = str(row[col]).split(";")[i].strip()
            expanded_rows.append(new_row)

    # Création de la DataFrame étendue
    df_expanded = pd.DataFrame(expanded_rows)
    df_expanded.to_csv(matrix_file, sep="\t", index=False)
    print(f"Fichier normalisé créé")
    return df_expanded


# -------------------------------
# PARSEUR DE FICHIER HMM
# -------------------------------
#quand virgule addition car veut dire que le mag a plusieurs copies/versions d'un gène
#exemple dictionnaire sortie: {"k141_123": 611.1, "k141_456": 607.1, ...}
def parse_hmm_scores(hmm_dir, hmm_base_name):
    path = os.path.join(hmm_dir, f"{hmm_base_name}.total.hmmsearch_result.txt")
    scores = {}
    if not os.path.exists(path):
        return scores
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            #continue veut dire on ignore, pas pris en compte (.strip() enlève les espaces inutiles)
            parts = line.strip().split()
            if len(parts) >= 6:
                gene_id = parts[0]
                try:
                    score = float(parts[5])
                    scores[gene_id] = score
                except ValueError:
                    continue
    return scores

# -------------------------------
# SUPPRESSION DES DOUBLONS PAR MAG
#gène présent plus d’une fois:garder uniquement celui avec le meilleur bitscore (ou un au hasard si égalité). Les autres sont remplacés par "0"
# -------------------------------
def resolve_duplicates(df, hits_columns, hmm_dir):
    for col in hits_columns:
        gene_to_rows = {}
        #{'k141_3370381_length_54765_cov_16.4124_15': [6], 'k141_295052_length_26580_cov_18.0000_22': [7]}
        
        for idx, val in df[col].items():
            #ici idx est bien à quelle ligne et pas prise en compte du "orginal_index"
            genes = [g.strip() for g in str(val).split(",") if g.startswith("k141")]
            for gene in genes:
                gene_to_rows.setdefault(gene, []).append(idx)
        #en gros cette boucle parcours la col, si 0 :[] s'affiche et si un gène est trouvé lors du parcours :il est ajouté avec sa position dans le dico
                
        for gene, indices in gene_to_rows.items():
            #si un gène présent qu'un fois sans doublons on passe outre
            if len(indices) <= 1:
                continue
            bitscores = []
            #[] pour chaque doublon: [(50.3, 235), (85.4, 236)]
            for idx in indices:
                hmm_base = df.at[idx, "Hmm.file"]
                #hmm_file.replace(".hmm", "")
                hmm_scores = parse_hmm_scores(hmm_dir, hmm_base)
                score = hmm_scores.get(gene, 0)
                #0 est la valeur a recup? faut faire .get car on récupère touts les gènes dans le hmm recherché
                bitscores.append((score, idx))
            max_score = max(score for score, _ in bitscores)
            #Le _ est une convention pour dire « je n’ai pas besoin de cette valeur » !!!!!!!!!!!!!!!!!!!!!!!!
            top_indices = [i for score, i in bitscores if score == max_score]
            keep_idx = random.choice(top_indices) if len(top_indices) > 1 else top_indices[0]
            #print(max_score,top_indices)= 1494.4 [133]
            for idx in indices:
                #print(idx, keep_idx) = 196 133 = 0; si à la ligne 196 ça dit de garder 133 alors mettre 0 à 196. si 133 133 c'est qu'on peut keep gene
                if idx != keep_idx:
                    val = str(df.at[idx, col])
                    genes = val.split(",")
                    cleaned = ["0" if g == gene else g for g in genes]
                    #['0', 'k141_7640239_length_32867_cov_18.3140_25']
                    df.at[idx, col] = ",".join(cleaned)
                    #0,k141_7640239_length_32867_cov_18.3140_25
    return df


def load_abundance_dicts(abundance_files):
    abundance_dicts = {}
    for rep, file in abundance_files.items():
        df = pd.read_csv(file, sep="\t")
        abundance_dicts[rep] = dict(zip(df["Name"], df["TPM"]))
    return abundance_dicts


def compute_sums(df, replicates):
    for rep in replicates:
        cols = [c for c in df.columns if c.endswith(f"_TPM_{rep}")]

        if len(cols) == 0:
            print(f"[WARNING] aucune colonne TPM trouvée pour {rep}")
            df[f"sum_TPM_{rep}"] = 0
        else:
            print(f"{rep}: {len(cols)} colonnes trouvées")
            df[f"sum_TPM_{rep}"] = df[cols].sum(axis=1)

    return df

def compute_abundance_list(gene_list, abundance_dict):
    if pd.isna(gene_list) or gene_list == "0":
        return 0
    return sum(
        abundance_dict.get(g.strip(), 0)
        for g in gene_list.split(",")
        if g.strip() != "0" and g.startswith("k141")
    )


def compute_all_abundances(df, hits_columns, abundance_dicts):
    all_new_cols = []

    for rep, abundance_dict in abundance_dicts.items():

        data = {
            f"{col}_TPM_{rep}": df[col].apply(
                lambda x: compute_abundance_list(x, abundance_dict)
            )
            for col in hits_columns
        }

        tmp_df = pd.DataFrame(data, index=df.index)
        all_new_cols.append(tmp_df)

    df = pd.concat([df] + all_new_cols, axis=1)

    return df

def drop_columns(df, cols):
    # supprimer colonnes .Hits
    df = df.drop(columns=cols)
    return df

def reaggregate(df, replicates):
    # colonnes d'abondance à sommer
    abundance_cols = [f"sum_TPM_{rep}" for rep in replicates]

    # colonnes TPM détaillées
    #tpm_cols = [c for c in df.columns if "_TPM_" in c]

    # autres colonnes (metadata)
    metadata_cols = [
        c for c in df.columns
        if c not in abundance_cols + ["original_index"]
        #if c not in abundance_cols + tpm_cols + ["original_index"]
    ]

    agg_dict = {}

    def merge_hmm(values):
        all_vals = []
        for v in values:
            all_vals.extend([x.strip() for x in str(v).split(",") if x != "0"])
        return ",".join(sorted(set(all_vals)))

    for col in metadata_cols:
        if col == "Hmm.file":
            agg_dict[col] = merge_hmm
        else:
            agg_dict[col] = "first"

    # sommer TPM
    #for col in tpm_cols:
        #agg_dict[col] = "sum"

    # sommer abondances finales
    for col in abundance_cols:
        agg_dict[col] = "sum"

    grouped = df.groupby("original_index", dropna=False).agg(agg_dict).reset_index()

    return grouped



# -------------------------------
# PIPELINE PRINCIPAL
# -------------------------------

def main():
    meta_id = "Y320"
    input_file=f"/shared/projects/feob_iron_oxidizing_bacteria_/META/ANALYSIS_MG/5_MAG_treatment/nb_genomes_per_pathway/matrices/{meta_id}_ArticleGenesPathway_metabolic.txt"
    matrix_file= f"/shared/projects/feob_iron_oxidizing_bacteria_/META/ANALYSIS_MG/5_MAG_treatment/nb_genomes_per_pathway/matrices/{meta_id}_ArticleGenesPathway_metabolic_normalized.txt"
    hmm_dir = f"/shared/projects/feob_iron_oxidizing_bacteria_/META/ANALYSIS_MG/5_MAG_treatment/RESULTS_METABOLIC-{meta_id}/intermediate_files/Hmmsearch_Outputs/"

    output_file_bin = f"/shared/projects/feob_iron_oxidizing_bacteria_/META/ANALYSIS_MG/5_MAG_treatment/nb_genomes_per_pathway/script_outputs/abund_{meta_id}_metabolic_detailed_bin_final.csv"
    output_file= f"/shared/projects/feob_iron_oxidizing_bacteria_/META/ANALYSIS_MG/5_MAG_treatment/nb_genomes_per_pathway/script_outputs/abund_{meta_id}_metabolic_expanded_final.csv"
    output_file2= f"/shared/projects/feob_iron_oxidizing_bacteria_/META/ANALYSIS_MG/5_MAG_treatment/nb_genomes_per_pathway/script_outputs/abund_{meta_id}_metabolic_reaggregated_final.csv"

    abundance_files = {
    "A": f"/shared/projects/feob_iron_oxidizing_bacteria_/META/ANALYSIS_MG/4_genes_contigs/salmon_quant/{meta_id}/{meta_id}-A_quantified/quant.sf",
    "B": f"/shared/projects/feob_iron_oxidizing_bacteria_/META/ANALYSIS_MG/4_genes_contigs/salmon_quant/{meta_id}/{meta_id}-B_quantified/quant.sf",
    "C": f"/shared/projects/feob_iron_oxidizing_bacteria_/META/ANALYSIS_MG/4_genes_contigs/salmon_quant/{meta_id}/{meta_id}-C_quantified/quant.sf"
}
    random.seed(42)

    # 1. créer une ligne par hmm pour homogénéisé l'analyse
    df_expanded = normalize_matrix(input_file, matrix_file)
    # garder index original pour reaggregation
    df = df_expanded.reset_index(drop=True)
    hits_columns = [col for col in df.columns if col.endswith(".Hits")]
    # 2. Supprimer redondances
    df = resolve_duplicates(df, hits_columns, hmm_dir)

    # charger abundances
    abundance_dicts = load_abundance_dicts(abundance_files)

    # calcul TPM par replicate
    df = compute_all_abundances(df, hits_columns, abundance_dicts)
    print(" Remplacement par abondance TPM OK")
    #print([c for c in df.columns if "_TPM_" in c][:10])
    #print(df.columns)
    #tab details per bin
    df_bin_clean = drop_columns(df.copy(), hits_columns)
    df_bin_clean.to_csv(output_file_bin, sep="\t", index=False)
    print("Fichier with details per bin generated")

    # calcul sommes A/B/C
    df = compute_sums(df, abundance_dicts.keys())
    #df, hits_columns,
    print("Somme par replicats calculated OK")

    # OUTPUT EXPANDED without itit columns bin.*.Hits and detail of TMP per bin
    bin_columns = [col for col in df.columns if col.startswith("bin.")]
    df_expanded = drop_columns(df.copy(),bin_columns)
    df_expanded.to_csv(output_file, sep="\t", index=False)
    print("Fichier résultats expanded disponible")

    # OUTPUT REAGGREGATED
    df_reaggregated = reaggregate(df_expanded, abundance_dicts.keys())
    df_reaggregated.to_csv(output_file2, sep="\t", index=False)
    print("Fichier résultats reaggregated disponible")

if __name__ == "__main__":
    main()

