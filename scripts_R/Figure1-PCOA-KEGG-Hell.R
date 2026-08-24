#PCOA Energy-KEGG

# Installation et chargement des packages
install.packages("ggrepel")
install.packages("ggforce")
install.packages("pheatmap")
install.packages("devtools")
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
if (!requireNamespace("edgeR", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("edgeR")
}
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

# Charger les bibliothèques nécessaires
library(readxl)
library(vegan)
library(ape)
library(ggplot2)
library(ggforce)
library(pheatmap)
library(dplyr)
library(edgeR)
library(openxlsx)
library(pairwiseAdonis)
library(cowplot)
library(patchwork)


# 1. Lire le fichier Excel
donnees <- read_excel("PCOA-Energy-KEGG-Hell.xlsx")

# 2. Extraire les noms des MAG (première colonne) et des échantillons (colonnes suivantes)
noms_mag <- donnees[[1]]  # Noms des MAG (première colonne)
noms_echantillons <- colnames(donnees)[-1]  # Noms des échantillons (colonnes 2 à n)

# 3. Créer la matrice d'abondances (lignes = MAG, colonnes = échantillons)
matrice_abondances <- as.matrix(donnees[, -1])  # Toutes les colonnes sauf la première
rownames(matrice_abondances) <- noms_mag  # Noms des MAG en lignes
colnames(matrice_abondances) <- noms_echantillons  # Noms des échantillons en colonnes

# Vérifier la structure
print(dim(matrice_abondances))
print(colnames(matrice_abondances))

# 4. Calculer la distance de Bray-Curtis (transposer pour que les échantillons soient en lignes)
distance <- vegan::vegdist(t(matrice_abondances), method = "bray")

# 5. Réaliser la PCoA
pcoa <- ape::pcoa(distance)

# 6. Extraire les coordonnées et assigner les noms des échantillons
coords <- pcoa$vectors
#rownames(coords) <- noms_echantillons  # Assigner les noms des échantillons aux coordonnées
rownames(coords) <- colnames(matrice_abondances)  # Assigner les noms des échantillons aux coordonnées



# Créer un data frame pour ggplot2
df_pcoa <- data.frame(
  PC1 = coords[, 1],
  PC2 = coords[, 2],
  Echantillon = rownames(coords)
)

# Extraire le nom du site (ex: "CAP" de "CAP-Rep1")
df_pcoa$Site <- sapply(strsplit(df_pcoa$Echantillon, "-"), function(x) x[1])

# Extraire le nom du site (ex: "CAP-Rep" de "CAP-Rep1")
#df_pcoa$Site <- gsub("-Rep[0-9]", "", df_pcoa$Site)

# Vérifier les noms des sites
print(unique(df_pcoa$Site))
print(df_pcoa$Echantillon)

# Extraire les valeurs propres
eigenvalues <- pcoa$values$Eigenvalues

# Afficher les valeurs propres
print(eigenvalues)

# Calculer la variance totale
total_variance <- sum(eigenvalues)

# Calculer les pourcentages de variance pour chaque axe
percentage_variance <- (eigenvalues / total_variance) * 100

# Afficher les pourcentages de variance pour les deux premiers axes
percentage_PC1 <- round(percentage_variance[1], 1)
percentage_PC2 <- round(percentage_variance[2], 1)

# Afficher les résultats
cat("PC1 explique", percentage_PC1, "% de la variance totale.\n")
cat("PC2 explique", percentage_PC2, "% de la variance totale.\n")

# Définir les symboles pour chaque site
symbols <- c("CAP" = 15, "NTE" = 16, "LL" = 17, "Y3" = 18)  # 15 = carré, 16 = cercle, 17 = triangle, 18 = losange

# Créer le graphique
p <- ggplot(df_pcoa, aes(x = PC1, y = PC2, color = Site, shape = Site)) +
  geom_point(color = "black",size = 2) +
  geom_text_repel(aes(label = Echantillon), size =2, box.padding = 0.5, color = "black") +  # Affiche les noms des échantillons
  ggforce::geom_mark_ellipse(aes(group = Site), level = 0.95, color = "black", linewidth = 0.2) +  # Ellipses de confiance en noir
  scale_shape_manual(values = symbols) +
  scale_color_grey(start = 0.2, end = 0.8) +  # Nuances de gris pour les couleurs
  labs(
    #title = "Analyse en Coordonnées Principales (PCoA)",
    x = paste0("PC1 (", percentage_PC1, "%)"),
    y = paste0("PC2 (", percentage_PC2, "%)"),
    shape = "Site",
    #color = "Site"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size =6),   # Taille des valeurs sur l'axe x
    axis.text.y = element_text(size = 6),   # Taille des valeurs sur l'axe y
    axis.title.x = element_text(size = 8), # Taille de la légende de l'axe x
    axis.title.y = element_text(size = 8)  # Taille de la légende de l'axe y
  )

# Ajouter le "B" en dehors du cadre
p_with_letterB<- ggdraw(p) +
  draw_label("B", x = 0.05, y = 0.97, size = 12, fontface = "bold")

# Afficher le graphique
print(p_with_letterB)

# Sauvegarder le graphique
ggsave("pcoa_ENERGY_KEGG-Hell.pdf", p_with_letterB, width = 8, height = 6, dpi = 300)

# Afficher la structure de l'objet pcoa
str(pcoa)



# Ouvrir une connexion vers un fichier texte
file_path <- "resultats_PCoA-Energy-KEGG-Hell.txt"
sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
cat("RESULTATS PcOA ENERGY-KEGG-Hell ")
cat("\n\n% Variance totale sur les 2 axes PC1 et PC2\n\n")
cat("PC1 explique", percentage_PC1, "% de la variance totale.\n")
cat("PC2 explique", percentage_PC2, "% de la variance totale.\n")
# Fermer la connexion pour revenir à la console
sink()  # Arrête la redirection

##Test statistique sur la PCOA

#Anosim

# Créer un vecteur de groupes (ex: "CAP20", "NTE20", etc.)
# Exemple : Si vos échantillons sont nommés "CAP20-A", "CAP20-B", etc.
groupes <- gsub("-.*", "", rownames(coords))  # Extrait "CAP20" de "CAP20-A"

# Exécuter l'ANOSIM
anosim_result <- anosim(distance, groupes, permutations = 999)
print(anosim_result)

#PERMANOVA

# Créer un data frame avec les groupes
groupes_df <- data.frame(Site = groupes)

# Exécuter la PERMANOVA
adonis_result <- adonis2(distance ~ Site, data = groupes_df, permutations = 999)
print(adonis_result)

# Ouvrir une connexion vers un fichier texte
file_path <- "resultats_PCoA-Energy-KEGG-Hell.txt"
sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
cat("\n\nTest satistique sur la PCoA\n\n ")
print(anosim_result)
print(adonis_result)
# Fermer la connexion pour revenir à la console
sink()  # Arrête la redirection


##Analyse de corrélation entre les abondances des MAG et les coordonnées de la PCoA.

# Extraire les coordonnées de la PCoA (déjà calculées)
coords <- pcoa$vectors
rownames(coords) <- colnames(matrice_abondances)  # Noms des échantillons

# Transposer la matrice d'abondances pour avoir les MAG en colonnes
#matrice_transposee <- t(matrice_abondances)


# Calculer les corrélations de Spearman entre les abondances des MAG et les axes de la PCoA
cor_pcoa_mag <- cor(t(matrice_abondances), coords, method = "spearman")

# Afficher les corrélations pour PC1
print(cor_pcoa_mag[, 1])

# Afficher les corrélations pour PC2
print(cor_pcoa_mag[, 2])

# Extraire les MAG les plus corrélés avec PC1 et PC2, avec leurs valeurs
top_mag_pc1_with_cor <- head(sort(abs(cor_pcoa_mag[, 1]), decreasing = TRUE), 5)
top_mag_pc2_with_cor <- head(sort(abs(cor_pcoa_mag[, 2]), decreasing = TRUE), 5)

# Afficher les résultats
print(top_mag_pc1_with_cor)
print(top_mag_pc2_with_cor)

# Ouvrir une connexion vers un fichier texte
file_path <- "resultats_PCoA-Energy-KEGG-Hell.txt"
sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
cat("\n\nTest de correlation Abondances/axes PcOA\n\n ")
cat("Axe PC1\n ")
print(top_mag_pc1_with_cor)
cat("Axe PC2\n ")
print(top_mag_pc2_with_cor)
# Fermer la connexion pour revenir à la console
sink()  # Arrête la redirection


### 5 métabolismes les plus abondants par site
# Extraire les noms des sites à partir des noms de colonnes
sites <- sub("-Rep[0-9]+$", "", colnames(matrice_abondances))

# Transformer la matrice en format long
df_long <- as.data.frame(matrice_abondances) %>%
  rownames_to_column(var = "Metabolisme") %>%  # Les noms de lignes deviennent la colonne MAG
  pivot_longer(
    cols = -Metabolisme,
    names_to = "Site_Replicat",
    values_to = "Abondance"
  ) %>%
  separate(Site_Replicat, into = c("Site", "Replicat"), sep = "-Rep", remove = FALSE)

# Vérifier que df_long est correctement créé
head(df_long)

# Calculer l'abondance moyenne par métabolisme et par site
abondance_moyenne <- df_long %>%
  group_by(Metabolisme, Site) %>%
  summarise(Abondance_Moyenne = mean(Abondance), .groups = "drop")

# Identifier les 5 métabolismes les plus abondants pour chaque site
top5_metabolismes_par_site <- abondance_moyenne %>%
  group_by(Site) %>%
  top_n(5, Abondance_Moyenne) %>%
  arrange(Site, desc(Abondance_Moyenne))

# Afficher les résultats
print(top5_metabolismes_par_site)

# Ouvrir une connexion vers un fichier texte
file_path <- "resultats_PCoA-Energy-KEGG-Hell.txt"
sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
cat("\n\n5 KEGG les plus abondants par site\n\n ")
print(top5_metabolismes_par_site)
# Fermer la connexion pour revenir à la console
sink()  # Arrête la redirection


# Fonction pour tester si les 5 métabolismes les plus abondants sont significativement plus abondants que les autres
test_significativite_metabolismes <- function(df_long, site) {
  # Filtrer les données pour le site spécifique
  df_site <- df_long %>% filter(Site == site)
  
  # Calculer l'abondance moyenne par métabolisme
  abondance_moyenne <- df_site %>%
    group_by(Metabolisme) %>%
    summarise(Abondance_Moyenne = mean(Abondance), .groups = "drop")
  
  # Identifier les 5 métabolismes les plus abondants
  top5 <- abondance_moyenne %>%
    top_n(5, Abondance_Moyenne) %>%
    pull(Metabolisme)
  
  # Créer une colonne pour indiquer si le métabolisme est dans le top 5
  df_site <- df_site %>%
    mutate(Top5 = ifelse(Metabolisme %in% top5, "Top5", "Autres"))
  
  # Test de normalité (Shapiro-Wilk)
  normalite <- df_site %>%
    group_by(Top5) %>%
    shapiro_test(Abondance)
  
  # Test d'homogénéité des variances (Levene)
  homogenite <- levene_test(df_site, Abondance ~ Top5)
  
  # Choisir le test approprié en fonction de la normalité et de l'homogénéité des variances
  if (all(normalite$p > 0.05) & homogenite$p > 0.05) {
    # Si les données sont normales et les variances homogènes, utiliser ANOVA + Tukey
    test_result <- df_site %>%
      anova_test(Abondance ~ Top5) %>%
      get_anova_table()
    
    posthoc_result <- df_site %>%
      tukey_hsd(Abondance ~ Top5)
    
    cat("\nTest ANOVA pour le site", site, ":\n")
    print(test_result)
    
    cat("\nTest post-hoc de Tukey pour le site", site, ":\n")
    print(posthoc_result)
    
    # Ouvrir une connexion vers un fichier texte
    file_path <- "resultats_PCoA-Energy-KEGG-Hell.txt"
    sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
    cat("\n\nTest de significativité abandnces ANOVA+Tukey\n\n ")
    cat("\nTest ANOVA pour le site", site, ":\n")
    print(test_result)
    cat("\nTest post-hoc de Tukey pour le site", site, ":\n")
    print(posthoc_result)
    # Fermer la connexion pour revenir à la console
    sink()  # Arrête la redirection
    
  } else {
    # Si les données ne sont pas normales ou les variances ne sont pas homogènes, utiliser Kruskal-Wallis + Dunn
    test_result <- df_site %>%
      kruskal_test(Abondance ~ Top5)
    
    posthoc_result <- df_site %>%
      dunn_test(Abondance ~ Top5, p.adjust.method = "bonferroni")
    
    cat("\nTest de Kruskal-Wallis pour le site", site, ":\n")
    print(test_result)
    
    cat("\nTest post-hoc de Dunn pour le site", site, ":\n")
    print(posthoc_result)
    
    # Ouvrir une connexion vers un fichier texte
    file_path <- "resultats_PCoA-Energy-KEGG-Hell.txt"
    sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
    cat("\n\nTest de significativité abandnces Kruskal Wallis+Dunn\n\n ")
    cat("\nTest de Kruskal-Wallis pour le site", site, ":\n")
    print(test_result)
    cat("\nTest post-hoc de Dunn pour le site", site, ":\n")
    print(posthoc_result)
    # Fermer la connexion pour revenir à la console
    sink()  # Arrête la redirection
    
  }
}

# Appliquer le test pour chaque site
sites_uniques <- unique(df_long$Site)
for (site in sites_uniques) {
  test_significativite_metabolismes(df_long, site)
}

#visualisation boxplot

# Fonction pour retirer les outliers
remove_outliers <- function(data, variable) {
  Q1 <- quantile(data[[variable]], 0.25, na.rm = TRUE)
  Q3 <- quantile(data[[variable]], 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  data %>%
    filter(.data[[variable]] >= lower_bound, .data[[variable]] <= upper_bound)
}

# Préparer les données pour chaque site
sites <- unique(df_long$Site)
plots <- list()

for (site in sites) {
  df_site <- df_long %>%
    filter(Site == site)
  
  abondance_moyenne_site <- df_site %>%
    group_by(Metabolisme) %>%
    summarise(Abondance_Moyenne = mean(Abondance), .groups = "drop")
  
  top5_site <- abondance_moyenne_site %>%
    top_n(5, Abondance_Moyenne) %>%
    pull(Metabolisme)
  
  # Ajouter la colonne Top5
  df_site <- df_site %>%
    mutate(Top5 = ifelse(Metabolisme %in% top5_site, "Top5", "Autres"))
  
  # Retirer les outliers
  df_site_no_outliers <- remove_outliers(df_site, "Abondance")
  
  # Calculer la valeur p pour le test de Kruskal-Wallis
  p_value <- df_site_no_outliers %>%
    kruskal_test(Abondance ~ Top5) %>%
    .$p
  
  p_text <- ifelse(p_value < 0.0001, "p < 0.0001", paste0("p = ", round(p_value, 4)))
  
  p <- ggplot(df_site_no_outliers, aes(x = Top5, y = Abondance, fill = Top5)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_grey(start = 0.7, end = 0.3) +
    labs(title = paste("Site", site),
         x = "Group",
         y = "Abundance") +
    annotate("text", x = 1.5, y = Inf, label = p_text, vjust = 2, hjust = 0.5) +
    theme_minimal() +
    theme(legend.position = "none",  # Retire la légende supérieure
          plot.title = element_text(hjust = 0.5))
  
  plots[[site]] <- p
}

# Combiner les graphiques en une seule feuille
combined_plot <- ggarrange(plotlist = plots, ncol = 2, nrow = 2, common.legend = FALSE)

# Afficher le graphique combiné
print(combined_plot)

# Sauvegarder le graphique
ggsave(filename = "abondance_top5_ENERGY_KEGG-Hell.pdf", plot = combined_plot, device = "pdf", width = 8, height = 6, dpi = 300)


### Analyse EdgR


# -------- Entrées --------
# matrice_abondances : matrice des comptages (lignes = features, colonnes = échantillons)
# groupes : vecteur facteur ou caractère des groupes des échantillons

# Exemple d'import :
# matrice_abondances <- read.csv("data_counts.csv", row.names=1)
# groupes <- factor(c("CAP", "CAP", "NTE", "NTE", "LL", "LL", "Y3", "Y3"))  

# 2. Fonction analyse edgeR avec correction Bonferroni

analyze_edgeR_pair_bonferroni <- function(counts_matrix, groupes_vector, groupe1, groupe2) {
  
  cat("Analyse edgeR Bonferroni pour :", groupe1, "vs", groupe2, "\n")
  
  samples_keep <- which(groupes_vector %in% c(groupe1, groupe2))
  counts_sub <- counts_matrix[, samples_keep]
  group_sub <- factor(groupes_vector[samples_keep], levels = c(groupe1, groupe2))
  
  dge <- DGEList(counts = counts_sub, group = group_sub)
  
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  
  if(nrow(dge) == 0) {
    warning("Aucun gène/pathway après filtrage pour cette comparaison.")
    return(NULL)
  }
  
  dge <- calcNormFactors(dge)
  design <- model.matrix(~ group_sub)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef=2)
  res <- topTags(qlf, n=Inf)$table
  
  res$Bonferroni <- p.adjust(res$PValue, method = "bonferroni")
  
  n_sig_bonf <- sum(res$Bonferroni < 0.05, na.rm = TRUE)
  cat("Nombre d'espèces/pathways significatifs (Bonferroni < 0.05) :", n_sig_bonf, "\n\n")
  
  # Ouvrir une connexion vers un fichier texte
  file_path <- "resultats_PCoA-Energy-KEGG-Hell.txt"
  sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
  cat("\n\nAbondances différentielle KEGG EdgeR\n\n ")
  cat("Nombre d'espèces/pathways significatifs (Bonferroni < 0.05) :", n_sig_bonf, "\n\n")
  # Fermer la connexion pour revenir à la console
  sink()  # Arrête la redirection
  
  return(res)
}

# 3. Analyse toutes paires de groupes

group_levels <- unique(groupes)
pairs <- combn(group_levels, 2, simplify = FALSE)
res_list_bonf <- list()

for(pair in pairs) {
  res_pair <- analyze_edgeR_pair_bonferroni(matrice_abondances, groupes, pair[1], pair[2])
  pair_name <- paste(pair, collapse = "_vs_")
  res_list_bonf[[pair_name]] <- res_pair
}

# 4. Affichage console

for(name in names(res_list_bonf)) {
  cat("\n=== Résultats pour :", name, "===\n")
  res_df <- res_list_bonf[[name]]
  if(is.null(res_df)) {
    cat("Aucun résultat disponible.\n")
  } else {
    print(head(res_df, n=20))
  }
}

# 5. Export dans un fichier Excel, une feuille par comparaison

wb <- createWorkbook()

for(name in names(res_list_bonf)) {
  res_df <- res_list_bonf[[name]]
  if(!is.null(res_df)) {
    # Remettre les noms de ligne dans une colonne
    res_df_out <- cbind(Feature = rownames(res_df), res_df)
    addWorksheet(wb, sheetName = substr(name, 1, 31))  
    # Limite Excel: nom de feuille max 31 chars, on tronque si besoin
    writeData(wb, sheet = substr(name, 1, 31), x = res_df_out)
  }
}

# Enregistrer le fichier Excel
output_filename <- "edgeR_Bonferroni_results_ENERGY-KEGG-Hell.xlsx"
saveWorkbook(wb, file = output_filename, overwrite = TRUE)

cat("\nTous les résultats ont été exportés dans :", output_filename, "\n")


### Tests statistiques sur PCOA et groupes

#Vérifier les groupes
table(groupes)

#distance en matrice
distance_matrix <- as.matrix(distance)  # ton objet distance mis en matrice

# Relancer la pairwise PERMANOVA
pairwise_results_Bonferroni <- pairwise.adonis(distance_matrix, groupes, p.adjust.m = "bonferroni", perm = 999)

cat("\navec Bonferroni:\n")
print(pairwise_results_Bonferroni)

# Ouvrir une connexion vers un fichier texte
file_path <- "resultats_PCoA-Energy-KEGG-Hell.txt"
sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
cat("\n\nTest de significativité différence groupes PcOA\n\n ")
cat("avec Bonferroni:\n")
print(pairwise_results_Bonferroni)
sink()  # Arrête la redirection

#autre ajustement FDR
pairwise_results_FDR <- pairwise.adonis(distance_matrix, groupes, p.adjust.m = "fdr", perm = 999)

cat("\navec FDR:\n")
print(pairwise_results_FDR)

# Ouvrir une connexion vers un fichier texte
file_path <- "resultats_PCoA-Energy-KEGG-Hell.txt"
sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
cat("\n\nTest de significativité différence groupes PcOA\n\n ")
cat("\navec FDR:\n")
print(pairwise_results_FDR)
# Fermer la connexion pour revenir à la console
sink()  # Arrête la redirection


#Groupe CAP et NTE car Proche sur la PcOA
groupes_fusion <- as.character(groupes)
groupes_fusion[groupes_fusion %in% c("NTE", "CAP")] <- "CAP_NTE"
groupes_fusion <- factor(groupes_fusion)

# Vérifie les groupes fusionnés
table(groupes_fusion)

# Relancer la pairwise PERMANOVA avec ce nouveau facteur
pairwise_results_fusion_Bonferroni <- pairwise.adonis(distance_matrix, groupes_fusion, p.adjust.m = "bonferroni", perm = 999)


cat("Fusion CAP-NTE avec bonferroni:\n")
print(pairwise_results_fusion_Bonferroni)

# Ouvrir une connexion vers un fichier texte
file_path <- "resultats_PCoA-Energy-KEGG-Hell.txt"
sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
cat("\n\nTest de significativité différence groupes PcOA\n\n ")
cat("Fusion CAP-NTE avec bonferroni:\n")
print(pairwise_results_fusion_Bonferroni)

# Fermer la connexion pour revenir à la console
sink()  # Arrête la redirection


#autre ajustement FDR
pairwise_results_fusion_FDR <- pairwise.adonis(distance_matrix, groupes_fusion, p.adjust.m = "fdr", perm = 999)

cat("Fusion CAP-NTE avec FDR:\n")
print(pairwise_results_fusion_FDR)

# Ouvrir une connexion vers un fichier texte
file_path <- "resultats_PCoA-Energy-KEGG-Hell.txt"
sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
cat("\n\nTest de significativité différence groupes PcOA\n\n ")
cat("Fusion CAP-NTE avec FDR:\n")
print(pairwise_results_fusion_FDR)

# Fermer la connexion pour revenir à la console
sink()  # Arrête la redirection

## visaulisation de plus d'axes de la PcoA

# 1. Calcul PCoA avec 4 axes
pcoa_result <- cmdscale(distance_matrix, k = 4, eig = TRUE)

# 2. Calcul % variance expliquée pour chaque axe
eig_vals <- pcoa_result$eig
variance_explained <- eig_vals / sum(eig_vals) * 100

# 3. Prépare dataframe pour ggplot
df_pcoa <- as.data.frame(pcoa_result$points)
colnames(df_pcoa) <- paste0("PC", 1:ncol(df_pcoa))
df_pcoa$group <- groupes_fusion  # Remplace par ton facteur de groupes

# 4. Fonction pour tracer deux axes donnés
plot_pcoa_axes <- function(df, x_axis, y_axis, variance, groups){
  xlab <- paste0(x_axis, " (", round(variance[as.numeric(substr(x_axis, 3, 3))], 1), "%)")
  ylab <- paste0(y_axis, " (", round(variance[as.numeric(substr(y_axis, 3, 3))], 1), "%)")
  
  ggplot(df, aes_string(x = x_axis, y = y_axis, color = "group")) +
    geom_point(size = 3) +
    labs(x = xlab, y = ylab, title = paste("PCoA :", x_axis, "vs", y_axis)) +
    theme_minimal() +
    theme(legend.title = element_blank())
}

# 5. Tracer plusieurs couples d’axes
plots <- list(
  plot_pcoa_axes(df_pcoa, "PC1", "PC2", variance_explained, groupes_fusion),
  plot_pcoa_axes(df_pcoa, "PC1", "PC3", variance_explained, groupes_fusion),
  plot_pcoa_axes(df_pcoa, "PC2", "PC3", variance_explained, groupes_fusion),
  plot_pcoa_axes(df_pcoa, "PC3", "PC4", variance_explained, groupes_fusion)
)

# Combine les 4 plots en une grille (2x2)
combined_plot2 <- ggarrange(plotlist = plots,
                            ncol = 2, nrow = 2,
                            common.legend = TRUE,
                            legend = "bottom")

# Affiche le plot combiné
print(combined_plot2)

# Sauvegarde en PDF
ggsave("PcOA_ENERGY_KEGG-Hell_4axes.pdf",
       plot = combined_plot2,
       width = 12,
       height = 10)



#===================================================================================
#===================================================================================
#PCOA ALL KEGG
#===================================================================================
#===================================================================================

# 1. Lire le fichier Excel
donnees <- read_excel("PCOA-all-KEGG-Hell.xlsx")

# 2. Extraire les noms des MAG (première colonne) et des échantillons (colonnes suivantes)
noms_mag <- donnees[[1]]  # Noms des MAG (première colonne)
noms_echantillons <- colnames(donnees)[-1]  # Noms des échantillons (colonnes 2 à n)

# 3. Créer la matrice d'abondances (lignes = MAG, colonnes = échantillons)
matrice_abondances <- as.matrix(donnees[, -1])  # Toutes les colonnes sauf la première
rownames(matrice_abondances) <- noms_mag  # Noms des MAG en lignes
colnames(matrice_abondances) <- noms_echantillons  # Noms des échantillons en colonnes

# Vérifier la structure
print(dim(matrice_abondances))
print(colnames(matrice_abondances))

# 4. Calculer la distance de Bray-Curtis (transposer pour que les échantillons soient en lignes)
distance <- vegan::vegdist(t(matrice_abondances), method = "bray")

# 5. Réaliser la PCoA
pcoa <- ape::pcoa(distance)

# 6. Extraire les coordonnées et assigner les noms des échantillons
coords <- pcoa$vectors
#rownames(coords) <- noms_echantillons  # Assigner les noms des échantillons aux coordonnées
rownames(coords) <- colnames(matrice_abondances)  # Assigner les noms des échantillons aux coordonnées


# 7. Créer le data frame pour ggplot2
df_pcoa <- data.frame(
  PC1 = coords[, 1],
  PC2 = coords[, 2],
  Echantillon = rownames(coords)
)

# Extraire le nom du site (ex: "CAP" de "CAP-Rep1")
df_pcoa$Site <- sapply(strsplit(df_pcoa$Echantillon, "-"), function(x) x[1])

# Extraire le nom du site (ex: "CAP-Rep" de "CAP-Rep1")
#df_pcoa$Site <- gsub("-Rep[0-9]", "", df_pcoa$Site)

# Vérifier les noms des sites
print(unique(df_pcoa$Site))
print(df_pcoa$Echantillon)

# Extraire les valeurs propres
eigenvalues <- pcoa$values$Eigenvalues

# Afficher les valeurs propres
print(eigenvalues)

# Calculer la variance totale
total_variance <- sum(eigenvalues)

# Calculer les pourcentages de variance pour chaque axe
percentage_variance <- (eigenvalues / total_variance) * 100

# Afficher les pourcentages de variance pour les deux premiers axes
percentage_PC1 <- round(percentage_variance[1], 1)
percentage_PC2 <- round(percentage_variance[2], 1)

# Afficher les résultats
cat("PC1 explique", percentage_PC1, "% de la variance totale.\n")
cat("PC2 explique", percentage_PC2, "% de la variance totale.\n")

# Définir les symboles pour chaque site
symbols <- c("CAP" = 15, "NTE" = 16, "LL" = 17, "Y3" = 18)  # 15 = carré, 16 = cercle, 17 = triangle, 18 = losange

# Créer le graphique
p <- ggplot(df_pcoa, aes(x = PC1, y = PC2, color = Site, shape = Site)) +
  geom_point(color = "black",size = 2) +
  geom_text_repel(aes(label = Echantillon), size = 2, box.padding = 0.5, color = "black") +  # Affiche les noms des échantillons
  ggforce::geom_mark_ellipse(aes(group = Site), color = "black", linewidth = 0.2) +  # Ellipses de confiance en noir
  scale_shape_manual(values = symbols) +
  scale_color_grey(start = 0.2, end = 0.8) +  # Nuances de gris pour les couleurs
  labs(
    #title = "Analyse en Coordonnées Principales (PCoA)",
    x = paste0("PC1 (", percentage_PC1, "%)"),
    y = paste0("PC2 (", percentage_PC2, "%)"),
    shape = "Site",
    #color = "Site"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size =6),   # Taille des valeurs sur l'axe x
    axis.text.y = element_text(size = 6),   # Taille des valeurs sur l'axe y
    axis.title.x = element_text(size = 8), # Taille de la légende de l'axe x
    axis.title.y = element_text(size = 8)  # Taille de la légende de l'axe y
  )

# Ajouter le "A" en dehors du cadre
p_with_letterA <- ggdraw(p) +
  draw_label("A", x = 0.05, y = 0.97, size = 12, fontface = "bold")

# Afficher le graphique
print(p_with_letterA)


# Sauvegarder le graphique
ggsave("pcoa_ALL_KEGG-Hell.pdf", p_with_letterA, width = 8, height = 6, dpi = 300)


# Afficher la structure de l'objet pcoa
str(pcoa)

# Ouvrir une connexion vers un fichier texte
file_path <- "resultats_PCoA-ALL-KEGG-Hell.txt"
sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
cat("RESULTATS PcOA ALL KEGG-Hell ")
cat("\n\n% Variance totale sur les 2 axes PC1 et PC2\n\n")
cat("PC1 explique", percentage_PC1, "% de la variance totale.\n")
cat("PC2 explique", percentage_PC2, "% de la variance totale.\n")
# Fermer la connexion pour revenir à la console
sink()  # Arrête la redirection




##Test statistique sur la PCOA

#Anosim

# Créer un vecteur de groupes (ex: "CAP20", "NTE20", etc.)
# Exemple : Si vos échantillons sont nommés "CAP20-A", "CAP20-B", etc.
groupes <- gsub("-.*", "", rownames(coords))  # Extrait "CAP20" de "CAP20-A"

# Exécuter l'ANOSIM
anosim_result <- anosim(distance, groupes, permutations = 999)
print(anosim_result)

#PERMANOVA

# Créer un data frame avec les groupes
groupes_df <- data.frame(Site = groupes)

# Exécuter la PERMANOVA
adonis_result <- adonis2(distance ~ Site, data = groupes_df, permutations = 999)
print(adonis_result)

# Ouvrir une connexion vers un fichier texte
file_path <- "resultats_PCoA-ALL-KEGG-Hell.txt"
sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
cat("\n\nTest satistique sur la PCoA\n\n ")
print(anosim_result)
print(adonis_result)
# Fermer la connexion pour revenir à la console
sink()  # Arrête la redirection


##Analyse de corrélation entre les abondances des MAG et les coordonnées de la PCoA.

# Extraire les coordonnées de la PCoA (déjà calculées)
coords <- pcoa$vectors
rownames(coords) <- colnames(matrice_abondances)  # Noms des échantillons

# Transposer la matrice d'abondances pour avoir les MAG en colonnes
#matrice_transposee <- t(matrice_abondances)


# Calculer les corrélations de Spearman entre les abondances des MAG et les axes de la PCoA
cor_pcoa_mag <- cor(t(matrice_abondances), coords, method = "spearman")

# Afficher les corrélations pour PC1
print(cor_pcoa_mag[, 1])

# Afficher les corrélations pour PC2
print(cor_pcoa_mag[, 2])

# Extraire les MAG les plus corrélés avec PC1 et PC2, avec leurs valeurs
top_mag_pc1_with_cor <- head(sort(abs(cor_pcoa_mag[, 1]), decreasing = TRUE), 5)
top_mag_pc2_with_cor <- head(sort(abs(cor_pcoa_mag[, 2]), decreasing = TRUE), 5)

# Afficher les résultats
print(top_mag_pc1_with_cor)
print(top_mag_pc2_with_cor)

# Ouvrir une connexion vers un fichier texte
file_path <- "resultats_PCoA-ALL-KEGG-Hell.txt"
sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
cat("\n\nTest de correlation Abondances/axes PcOA\n\n ")
print(top_mag_pc1_with_cor)
print(top_mag_pc2_with_cor)
# Fermer la connexion pour revenir à la console
sink()  # Arrête la redirection

### 5 métabolismes les plus abondants par site
# Extraire les noms des sites à partir des noms de colonnes
sites <- sub("-Rep[0-9]+$", "", colnames(matrice_abondances))

# Transformer la matrice en format long
df_long <- as.data.frame(matrice_abondances) %>%
  rownames_to_column(var = "Metabolisme") %>%  # Les noms de lignes deviennent la colonne MAG
  pivot_longer(
    cols = -Metabolisme,
    names_to = "Site_Replicat",
    values_to = "Abondance"
  ) %>%
  separate(Site_Replicat, into = c("Site", "Replicat"), sep = "-Rep", remove = FALSE)

# Vérifier que df_long est correctement créé
head(df_long)

# Calculer l'abondance moyenne par métabolisme et par site
abondance_moyenne <- df_long %>%
  group_by(Metabolisme, Site) %>%
  summarise(Abondance_Moyenne = mean(Abondance), .groups = "drop")

# Identifier les 5 métabolismes les plus abondants pour chaque site
top5_metabolismes_par_site <- abondance_moyenne %>%
  group_by(Site) %>%
  top_n(5, Abondance_Moyenne) %>%
  arrange(Site, desc(Abondance_Moyenne))

# Afficher les résultats
print(top5_metabolismes_par_site)

# Ouvrir une connexion vers un fichier texte
file_path <- "resultats_PCoA-ALL-KEGG-Hell.txt"
sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
cat("\n\n5 KEGG les plus abondants par site\n\n ")
print(top5_metabolismes_par_site)
# Fermer la connexion pour revenir à la console
sink()  # Arrête la redirection


# Fonction pour tester si les 5 métabolismes les plus abondants sont significativement plus abondants que les autres
test_significativite_metabolismes <- function(df_long, site) {
  # Filtrer les données pour le site spécifique
  df_site <- df_long %>% filter(Site == site)
  
  # Calculer l'abondance moyenne par métabolisme
  abondance_moyenne <- df_site %>%
    group_by(Metabolisme) %>%
    summarise(Abondance_Moyenne = mean(Abondance), .groups = "drop")
  
  # Identifier les 5 métabolismes les plus abondants
  top5 <- abondance_moyenne %>%
    top_n(5, Abondance_Moyenne) %>%
    pull(Metabolisme)
  
  # Créer une colonne pour indiquer si le métabolisme est dans le top 5
  df_site <- df_site %>%
    mutate(Top5 = ifelse(Metabolisme %in% top5, "Top5", "Autres"))
  
  # Test de normalité (Shapiro-Wilk)
  normalite <- df_site %>%
    group_by(Top5) %>%
    shapiro_test(Abondance)
  
  # Test d'homogénéité des variances (Levene)
  homogenite <- levene_test(df_site, Abondance ~ Top5)
  
  # Choisir le test approprié en fonction de la normalité et de l'homogénéité des variances
  if (all(normalite$p > 0.05) & homogenite$p > 0.05) {
    # Si les données sont normales et les variances homogènes, utiliser ANOVA + Tukey
    test_result <- df_site %>%
      anova_test(Abondance ~ Top5) %>%
      get_anova_table()
    
    posthoc_result <- df_site %>%
      tukey_hsd(Abondance ~ Top5)
    
    cat("\nTest ANOVA pour le site", site, ":\n")
    print(test_result)
    
    cat("\nTest post-hoc de Tukey pour le site", site, ":\n")
    print(posthoc_result)
    
    # Ouvrir une connexion vers un fichier texte
    file_path <- "resultats_PCoA-ALL-KEGG-Hell.txt"
    sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
    cat("\n\nTest de significativité abandnces ANOVA+Tukey\n\n ")
    cat("\nTest ANOVA pour le site", site, ":\n")
    print(test_result)
    cat("\nTest post-hoc de Tukey pour le site", site, ":\n")
    print(posthoc_result)
    # Fermer la connexion pour revenir à la console
    sink()  # Arrête la redirection
    
  } else {
    # Si les données ne sont pas normales ou les variances ne sont pas homogènes, utiliser Kruskal-Wallis + Dunn
    test_result <- df_site %>%
      kruskal_test(Abondance ~ Top5)
    
    posthoc_result <- df_site %>%
      dunn_test(Abondance ~ Top5, p.adjust.method = "bonferroni")
    
    cat("\nTest de Kruskal-Wallis pour le site", site, ":\n")
    print(test_result)
    
    cat("\nTest post-hoc de Dunn pour le site", site, ":\n")
    print(posthoc_result)
    
    # Ouvrir une connexion vers un fichier texte
    file_path <- "resultats_PCoA-ALL-KEGG-Hell.txt"
    sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
    cat("\n\nTest de significativité abandnces Kruskal Wallis+Dunn\n\n ")
    cat("\nTest de Kruskal-Wallis pour le site", site, ":\n")
    print(test_result)
    cat("\nTest post-hoc de Dunn pour le site", site, ":\n")
    print(posthoc_result)
    # Fermer la connexion pour revenir à la console
    sink()  # Arrête la redirection
    
  }
}

# Appliquer le test pour chaque site
sites_uniques <- unique(df_long$Site)
for (site in sites_uniques) {
  test_significativite_metabolismes(df_long, site)
}

#visualisation boxplot

# Fonction pour retirer les outliers
remove_outliers <- function(data, variable) {
  Q1 <- quantile(data[[variable]], 0.25, na.rm = TRUE)
  Q3 <- quantile(data[[variable]], 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  data %>%
    filter(.data[[variable]] >= lower_bound, .data[[variable]] <= upper_bound)
}

# Préparer les données pour chaque site
sites <- unique(df_long$Site)
plots <- list()

for (site in sites) {
  df_site <- df_long %>%
    filter(Site == site)
  
  abondance_moyenne_site <- df_site %>%
    group_by(Metabolisme) %>%
    summarise(Abondance_Moyenne = mean(Abondance), .groups = "drop")
  
  top5_site <- abondance_moyenne_site %>%
    top_n(5, Abondance_Moyenne) %>%
    pull(Metabolisme)
  
  # Ajouter la colonne Top5
  df_site <- df_site %>%
    mutate(Top5 = ifelse(Metabolisme %in% top5_site, "Top5", "Autres"))
  
  # Retirer les outliers
  df_site_no_outliers <- remove_outliers(df_site, "Abondance")
  
  # Calculer la valeur p pour le test de Kruskal-Wallis
  p_value <- df_site_no_outliers %>%
    kruskal_test(Abondance ~ Top5) %>%
    .$p
  
  p_text <- ifelse(p_value < 0.0001, "p < 0.0001", paste0("p = ", round(p_value, 4)))
  
  p <- ggplot(df_site_no_outliers, aes(x = Top5, y = Abondance, fill = Top5)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_grey(start = 0.7, end = 0.3) +
    labs(title = paste("Site", site),
         x = "Group",
         y = "Abundance") +
    annotate("text", x = 1.5, y = Inf, label = p_text, vjust = 2, hjust = 0.5) +
    theme_minimal() +
    theme(legend.position = "none",  # Retire la légende supérieure
          plot.title = element_text(hjust = 0.5))
  
  plots[[site]] <- p
}

# Combiner les graphiques en une seule feuille
combined_plot <- ggarrange(plotlist = plots, ncol = 2, nrow = 2, common.legend = FALSE)

# Afficher le graphique combiné
print(combined_plot)

# Sauvegarder le graphique
ggsave(filename = "abondance_top5_ALL_KEGG-Hell.pdf", plot = combined_plot, device = "pdf", width = 8, height = 6, dpi = 300)


### Analyse EdgR


# -------- Entrées --------
# matrice_abondances : matrice des comptages (lignes = features, colonnes = échantillons)
# groupes : vecteur facteur ou caractère des groupes des échantillons

# Exemple d'import :
# matrice_abondances <- read.csv("data_counts.csv", row.names=1)
# groupes <- factor(c("CAP", "CAP", "NTE", "NTE", "LL", "LL", "Y3", "Y3"))  

# 2. Fonction analyse edgeR avec correction Bonferroni

analyze_edgeR_pair_bonferroni <- function(counts_matrix, groupes_vector, groupe1, groupe2) {
  
  cat("Analyse edgeR Bonferroni pour :", groupe1, "vs", groupe2, "\n")
  
  samples_keep <- which(groupes_vector %in% c(groupe1, groupe2))
  counts_sub <- counts_matrix[, samples_keep]
  group_sub <- factor(groupes_vector[samples_keep], levels = c(groupe1, groupe2))
  
  dge <- DGEList(counts = counts_sub, group = group_sub)
  
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  
  if(nrow(dge) == 0) {
    warning("Aucun gène/pathway après filtrage pour cette comparaison.")
    return(NULL)
  }
  
  dge <- calcNormFactors(dge)
  design <- model.matrix(~ group_sub)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef=2)
  res <- topTags(qlf, n=Inf)$table
  
  res$Bonferroni <- p.adjust(res$PValue, method = "bonferroni")
  
  n_sig_bonf <- sum(res$Bonferroni < 0.05, na.rm = TRUE)
  cat("Nombre d'espèces/pathways significatifs (Bonferroni < 0.05) :", n_sig_bonf, "\n\n")
  
  # Ouvrir une connexion vers un fichier texte
  file_path <- "resultats_PCoA-ALL-KEGG-Hell.txt"
  sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
  cat("\n\nAbondances différentielle KEGG EdgeR\n\n ")
  cat("Nombre d'espèces/pathways significatifs (Bonferroni < 0.05) :", n_sig_bonf, "\n\n")
  # Fermer la connexion pour revenir à la console
  sink()  # Arrête la redirection
  
  return(res)
}

# 3. Analyse toutes paires de groupes

group_levels <- unique(groupes)
pairs <- combn(group_levels, 2, simplify = FALSE)
res_list_bonf <- list()

for(pair in pairs) {
  res_pair <- analyze_edgeR_pair_bonferroni(matrice_abondances, groupes, pair[1], pair[2])
  pair_name <- paste(pair, collapse = "_vs_")
  res_list_bonf[[pair_name]] <- res_pair
}

# 4. Affichage console

for(name in names(res_list_bonf)) {
  cat("\n=== Résultats pour :", name, "===\n")
  res_df <- res_list_bonf[[name]]
  if(is.null(res_df)) {
    cat("Aucun résultat disponible.\n")
  } else {
    print(head(res_df, n=20))
  }
}

# 5. Export dans un fichier Excel, une feuille par comparaison

wb <- createWorkbook()

for(name in names(res_list_bonf)) {
  res_df <- res_list_bonf[[name]]
  if(!is.null(res_df)) {
    # Remettre les noms de ligne dans une colonne
    res_df_out <- cbind(Feature = rownames(res_df), res_df)
    addWorksheet(wb, sheetName = substr(name, 1, 31))  
    # Limite Excel: nom de feuille max 31 chars, on tronque si besoin
    writeData(wb, sheet = substr(name, 1, 31), x = res_df_out)
  }
}

# Enregistrer le fichier Excel
output_filename <- "edgeR_Bonferroni_results_ALL-KEGG-Hell.xlsx"
saveWorkbook(wb, file = output_filename, overwrite = TRUE)

cat("\nTous les résultats ont été exportés dans :", output_filename, "\n")

### Tests statistiques sur PCOA et groupes

#Vérifier les groupes
table(groupes)

#distance en matrice
distance_matrix <- as.matrix(distance)  # ton objet distance mis en matrice

# Relancer la pairwise PERMANOVA
pairwise_results_Bonferroni <- pairwise.adonis(distance_matrix, groupes, p.adjust.m = "bonferroni", perm = 999)

cat("\navec Bonferroni:\n")
print(pairwise_results_Bonferroni)

# Ouvrir une connexion vers un fichier texte
file_path <- "resultats_PCoA-ALL-KEGG-Hell.txt"
sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
cat("\n\nTest de significativité différence groupes PcOA\n\n ")
cat("avec Bonferroni:\n")
print(pairwise_results_Bonferroni)
sink()  # Arrête la redirection

#autre ajustement FDR
pairwise_results_FDR <- pairwise.adonis(distance_matrix, groupes, p.adjust.m = "fdr", perm = 999)

cat("\navec FDR:\n")
print(pairwise_results_FDR)

# Ouvrir une connexion vers un fichier texte
file_path <- "resultats_PCoA-ALL-KEGG-Hell.txt"
sink(file_path, append = TRUE)  # Redirige toute la sortie vers le fichier
cat("\n\nTest de significativité différence groupes PcOA\n\n ")
cat("\navec FDR:\n")
print(pairwise_results_FDR)
# Fermer la connexion pour revenir à la console
sink()  # Arrête la redirection


## visaulisation de plus d'axes de la PcoA

# 1. Calcul PCoA avec 4 axes
pcoa_result <- cmdscale(distance_matrix, k = 4, eig = TRUE)

# 2. Calcul % variance expliquée pour chaque axe
eig_vals <- pcoa_result$eig
variance_explained <- eig_vals / sum(eig_vals) * 100

# 3. Prépare dataframe pour ggplot
df_pcoa <- as.data.frame(pcoa_result$points)
colnames(df_pcoa) <- paste0("PC", 1:ncol(df_pcoa))
df_pcoa$group <- groupes  # Remplace par ton facteur de groupes

# 4. Fonction pour tracer deux axes donnés
plot_pcoa_axes <- function(df, x_axis, y_axis, variance, groups){
  xlab <- paste0(x_axis, " (", round(variance[as.numeric(substr(x_axis, 3, 3))], 1), "%)")
  ylab <- paste0(y_axis, " (", round(variance[as.numeric(substr(y_axis, 3, 3))], 1), "%)")
  
  ggplot(df, aes_string(x = x_axis, y = y_axis, color = "group")) +
    geom_point(size = 3) +
    labs(x = xlab, y = ylab, title = paste("PCoA :", x_axis, "vs", y_axis)) +
    theme_minimal() +
    theme(legend.title = element_blank())
}

# 5. Tracer plusieurs couples d’axes
plots <- list(
  plot_pcoa_axes(df_pcoa, "PC1", "PC2", variance_explained, groupes),
  plot_pcoa_axes(df_pcoa, "PC1", "PC3", variance_explained, groupes),
  plot_pcoa_axes(df_pcoa, "PC2", "PC3", variance_explained, groupes),
  plot_pcoa_axes(df_pcoa, "PC3", "PC4", variance_explained, groupes)
)

# Combine les 4 plots en une grille (2x2)
combined_plot2 <- ggarrange(plotlist = plots,
                            ncol = 2, nrow = 2,
                            common.legend = TRUE,
                            legend = "bottom")

# Affiche le plot combiné
print(combined_plot2)

# Sauvegarde en PDF
ggsave("PcOA_ALL_KEGG-Hell_4axes.pdf",
       plot = combined_plot2,
       width = 12,
       height = 10)


#=============================
#combined PLOT Figure 1
#============================
# Combiner les deux graphiques
combined_plot <- p_with_letterA + p_with_letterB

# Afficher le résultat
print(combined_plot)

# Sauvegarder en TIFF
ggsave("Figure1-Hell.TIFF", plot = combined_plot, device = "tiff", dpi = 300, width = 16, height = 8, units = "cm")

ggsave("Figure1-Hell.png", plot = combined_plot, device = "png", dpi = 300, width = 16, height = 8, units = "cm")



