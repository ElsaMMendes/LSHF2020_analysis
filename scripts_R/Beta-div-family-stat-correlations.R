# Charger les bibliothèques nécessaires
library(dplyr)
library(vegan)
library(ggplot2)
library(ggrepel)
library(readxl)
library(tibble)
library(tidyverse)

# 1. Lire le fichier Excel contenant les données
data <- read_excel("Diversite-alpha-beta-MAGs.xlsx", sheet = "per replicat")

# 2. Préparer les données d'abondance
for (col in c("rep_1", "rep_2", "rep_3")) {
  data[[col]] <- gsub(",", ".", as.character(data[[col]]))
  data[[col]] <- suppressWarnings(as.numeric(data[[col]]))
  data[[col]][is.na(data[[col]])] <- 0
}

# 3. Créer une matrice de comptage large (échantillons en colonnes, MAGs en lignes)
count_data_long <- data %>%
  select(site, MAG, starts_with("rep_")) %>%
  pivot_longer(cols = starts_with("rep_"),
               names_to = "replicate",
               values_to = "abundance") %>%
  mutate(sample_id = paste(site, replicate, sep = "_")) %>%
  select(MAG, sample_id, abundance)

# 4. Ajouter les informations taxonomiques
taxonomic_data <- data %>%
  select(MAG, family) %>%
  drop_na(family) %>%
  distinct(MAG, .keep_all = TRUE)

# 5. Agrégation des abondances par family
family_abundance <- count_data_long %>%
  left_join(taxonomic_data, by = "MAG") %>%
  group_by(sample_id, family) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0)

# 6. Créer une matrice d'abondance avec les family en lignes et les échantillons en colonnes
family_abundance_matrix <- family_abundance %>%
  column_to_rownames(var = "family") %>%
  as.matrix()

# 7. Normalisation
family_abundance_norm <- family_abundance_matrix / colSums(family_abundance_matrix)

# 8. Calcul des distances de Bray-Curtis entre les échantillons
distance_matrix_family <- vegdist(t(family_abundance_norm), method = "bray")

# 9. Analyse PCoA
pcoa_result_family <- cmdscale(distance_matrix_family, k = 2, eig = TRUE) 
#x.ret = TRUE)
pcoa_df_family <- as.data.frame(pcoa_result_family$points)
colnames(pcoa_df_family) <- c("PC1", "PC2")

# Ajouter les noms des échantillons et les sites
pcoa_df_family$sample_id <- rownames(pcoa_result_family$points)
pcoa_df_family$Site <- sapply(strsplit(pcoa_df_family$sample_id, "_"), function(x) x[1])

# 10. Calcul des pourcentages de variance expliquée
eigen_values <- pcoa_result_family$eig[1:2]
total_variance <- sum(eigen_values)
percentage_PC1 <- round(eigen_values[1] / total_variance * 100, 1)
percentage_PC2 <- round(eigen_values[2] / total_variance * 100, 1)

cat("PC1 explique", percentage_PC1, "% de la variance totale.\n")
cat("PC2 explique", percentage_PC2, "% de la variance totale.\n")

# Ouvrir une connexion vers un fichier texte
file_path <- "family_correlations_beta_div.txt"
sink(file_path)  # Redirige toute la sortie vers le fichier

cat("PC1 explique", percentage_PC1, "% de la variance totale.\n")
cat("PC2 explique", percentage_PC2, "% de la variance totale.\n")

# Fermer la connexion pour revenir à la console
sink()  # Arrête la redirection



# Définir les symboles pour chaque site
symbols <- c("CAP" = 15, "NTE" = 16, "LL" = 17, "Y3" = 18)  # 15 = carré, 16 = cercle, 17 = triangle, 18 = losange


# Créer le graphique
p <- ggplot(pcoa_df_family, aes(x = PC1, y = PC2, shape = Site)) +
  geom_point(size = 3, color = "black") +  # Points en noir
  scale_shape_manual(values = symbols) +  # Utiliser les symboles définis
  geom_text_repel(aes(label = sample_id), size = 3, box.padding = 0.5, color = "black") +  # Affiche les noms des échantillons
  ggforce::geom_mark_ellipse(aes(group = Site), level = 0.95, color = "black", linewidth = 0.5) +  # Ellipses de confiance en noir
  #labs(
    #title = "PCoA (Bray-Curtis) - family Level",
   # x = paste0("PC1"),
   # y = paste0("PC2"),
    #shape = "Site"
 # ) +
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
  axis.text.x = element_text(size =10),   # Taille des valeurs sur l'axe x
axis.text.y = element_text(size = 10),   # Taille des valeurs sur l'axe y
axis.title.x = element_text(size = 12), # Taille de la légende de l'axe x
axis.title.y = element_text(size = 12))  # Taille de la légende de l'axe y# Supprime la légende

# Afficher le graphique
print(p)

# Sauvegarder le graphique
ggsave("pcoa_plot_family_symbols.pdf", p, width = 8, height = 6, dpi = 300)
ggsave("pcoa_plot_family_symbols.png", p, width = 8, height = 6, dpi = 300)

# 15. Tests statistiques
adonis_result_family <- adonis2(distance_matrix_family ~ Site, data = data.frame(Site = pcoa_df_family$Site))
cat("\nRésultats PERMANOVA (family Level):\n")
print(adonis_result_family)
write.csv(as.data.frame(t(as.matrix(adonis_result_family))), "permanova_results_beta_family.csv")

anosim_result_family <- anosim(distance_matrix_family, pcoa_df_family$Site)
cat("\nRésultats ANOSIM (family Level):\n")
print(anosim_result_family)
write.csv(data.frame(Statistic = anosim_result_family$R,
                     p.value = anosim_result_family$p.value),
          "anosim_results_beta_family.csv")

# 16. Dendrogramme de clustering
hclust_result_family <- hclust(distance_matrix_family, method = "ward.D2")
plot(hclust_result_family, cex = 0.6, hang = -1, main = "Clustering hiérarchique (Bray-Curtis) - family Level")
dendextend::color_branches(hclust_result_family, group = pcoa_df_family$Site)
dev.copy(pdf, "dendrogram_family.pdf", width = 10, height = 8)
dev.off()



# 11. Préparation des données pour les corrélations
family_abundance_norm_t <- t(family_abundance_norm)
rownames(family_abundance_norm_t) <- colnames(family_abundance_norm)

# Vérifiez que les noms des lignes de family_abundance_norm_t et les sample_id de pcoa_df_family sont les mêmes
all(rownames(family_abundance_norm_t) == pcoa_df_family$sample_id)

# Si les noms ne sont pas dans le même ordre, réorganisez-les
family_abundance_norm_t <- family_abundance_norm_t[match(pcoa_df_family$sample_id, rownames(family_abundance_norm_t)), ]

# Vérifiez les valeurs manquantes dans family_abundance_norm_t
any(is.na(family_abundance_norm_t))

# Remplacez les valeurs NA par 0
family_abundance_norm_t[is.na(family_abundance_norm_t)] <- 0

# Vérifiez les colonnes constantes
constant_columns <- apply(family_abundance_norm_t, 2, function(x) var(x) == 0)
table(constant_columns)

# Si des colonnes sont constantes, retirez-les
if (any(constant_columns)) {
  family_abundance_norm_t <- family_abundance_norm_t[, !constant_columns, drop = FALSE]
}

# Vérifiez les dimensions après nettoyage
dim(family_abundance_norm_t)

# 12. Calcul des corrélations et des p-values pour PC1 et PC2
PC1_correlations <- numeric(ncol(family_abundance_norm_t))
PC2_correlations <- numeric(ncol(family_abundance_norm_t))
PC1_pvalues <- numeric(ncol(family_abundance_norm_t))
PC2_pvalues <- numeric(ncol(family_abundance_norm_t))

names(PC1_correlations) <- colnames(family_abundance_norm_t)
names(PC2_correlations) <- colnames(family_abundance_norm_t)
names(PC1_pvalues) <- colnames(family_abundance_norm_t)
names(PC2_pvalues) <- colnames(family_abundance_norm_t)

# Calcul des corrélations et p-values pour PC1
for (i in seq_along(colnames(family_abundance_norm_t))) {
  family <- colnames(family_abundance_norm_t)[i]
  test_result <- try(cor.test(as.numeric(family_abundance_norm_t[, family]), as.numeric(pcoa_df_family$PC1), method = "spearman"), silent = TRUE)
  if (!inherits(test_result, "try-error")) {
    PC1_correlations[i] <- test_result$estimate
    PC1_pvalues[i] <- test_result$p.value
  } else {
    PC1_correlations[i] <- NA
    PC1_pvalues[i] <- NA
  }
}

# Calcul des corrélations et p-values pour PC2
for (i in seq_along(colnames(family_abundance_norm_t))) {
  family <- colnames(family_abundance_norm_t)[i]
  test_result <- try(cor.test(as.numeric(family_abundance_norm_t[, family]), as.numeric(pcoa_df_family$PC2), method = "spearman"), silent = TRUE)
  if (!inherits(test_result, "try-error")) {
    PC2_correlations[i] <- test_result$estimate
    PC2_pvalues[i] <- test_result$p.value
  } else {
    PC2_correlations[i] <- NA
    PC2_pvalues[i] <- NA
  }
}

# 13. Créer un data frame avec les corrélations et les p-values
significant_family <- data.frame(
  family = names(PC1_correlations),
  PC1_correlation = PC1_correlations,
  PC1_pvalue = PC1_pvalues,
  PC2_correlation = PC2_correlations,
  PC2_pvalue = PC2_pvalues,
  stringsAsFactors = FALSE
) %>%
  mutate(
    abs_PC1 = abs(PC1_correlation),
    abs_PC2 = abs(PC2_correlation),
    significant_PC1 = PC1_pvalue < 0.05,
    significant_PC2 = PC2_pvalue < 0.05
  ) %>%
  filter(significant_PC1 | significant_PC2)  # Filtrer les lignes où au moins une corrélation est significative

# Ouvrir une connexion vers un fichier texte
file_path <- "resultats_family_beta_div-top20.txt"
sink(file_path)  # Redirige toute la sortie vers le fichier

cat("PC1 explique", percentage_PC1, "% de la variance totale.\n")
cat("PC2 explique", percentage_PC2, "% de la variance totale.\n")

cat("\nRésultats PERMANOVA (family Level):\n")
print(adonis_result_family)
cat("\nRésultats ANOSIM (family Level):\n")
print(anosim_result_family)

# 14. Afficher les familys les plus corrélés avec PC1 et PC2 (significatifs uniquement)
cat("\nfamilys significativement corrélés avec PC1 (top 20) :\n")
print(
  significant_family %>%
    filter(significant_PC1) %>%
    select(family, PC1_correlation, PC1_pvalue) %>%
    arrange(desc(abs(PC1_correlation))) %>%
    head(20)
)

cat("\nfamilys significativement corrélés avec PC2 (top 20) :\n")
print(
  significant_family %>%
    filter(significant_PC2) %>%
    select(family, PC2_correlation, PC2_pvalue) %>%
    arrange(desc(abs(PC2_correlation))) %>%
    head(20)
)

# Fermer la connexion pour revenir à la console
sink()  # Arrête la redirection
