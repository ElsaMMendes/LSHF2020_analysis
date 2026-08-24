#analyse de clustering hiérarchique 

setwd("~/Documents/Documents Alain/ARTICLE/Metagenomic LSHF avril 2024/Papier Mendes et al mars 2025/tableau a faire")

# Installer les packages (une seule fois)
# install.packages(c("readxl", "vegan", "dendextend", "ggplot2", "reshape2", "factoextra"))
install.packages("pvclust")
library(readxl)
library(vegan)
library(dendextend)
library(ggplot2)
library(reshape2)
library(factoextra)
library(pvclust)
install.packages("ade4")
library(ade4)
library(phangorn)
library(compositions)  # Pour la normalisation CLR
library(ape)

# 1. Lire les données Excel (entrez le chemin correct du fichier)
data_raw <- read_excel("PCOA-MAG-Function-gMethoxi-RelAb2.xlsx")

# Vérifier la structure des données
head(data_raw)

# Supposons :
# - les colonnes sont les réplicats (nom colonnes : "CAP-Rep1", "CAP-Rep2", ...),
# - la première colonne est le nom des fonctions

# 2. Extraire la matrice d’abondance (enlever la colonne des noms)
# Adapter ici si la colonne des fonctions a un autre nom ou position
fun_names <- data_raw[[1]]        # Première colonne = noms fonctions
data_abund <- as.data.frame(data_raw[,-1])
rownames(data_abund) <- fun_names

# 3. Transposer la matrice pour avoir samples (réplicats) en LIGNES et fonctions en COLONNES
data_abund_t <- t(data_abund)


# Vérifier noms d’échantillons/colonnes
sample_names <- rownames(data_abund_t)
print(sample_names)

# 4. Extraire l’information site depuis les noms d’échantillons

# Exemple : "CAP-Rep1" -> site = "CAP"
library(stringr)
site <- str_extract(sample_names, "^[A-Z0-9]+") # extraction du début jusqu’au premier "-"
print(site)


# 5. Calculer la matrice de distance Bray-Curtis entre échantillons
dist_bc <- vegdist(data_abund_t, method = "bray")


# 6. Clustering hiérarchique (méthode ward.D2)
hc <- hclust(dist_bc, method = "ward.D2")
dist_bc <- as.matrix(dist_bc)
boot <- pvclust(dist_bc, method.hclust = "ward.D2", nboot = 1000)

plot(boot, print.pv = "au", print.num = FALSE)

# convertir cluster en objet dendrogramme
dend <- as.dendrogram(hc)

# associer une couleur au site
site_colors <- c(CAP = "#E41A1C", NTE = "#377EB8", LL = "#4DAF4A", Y3 = "#FF7F00")
colors_dend <- site_colors[site]

# appliquer couleurs sur les labels
labels_colors(dend) <- colors_dend

# tracer le dendrogramme
plot(dend, main = "Clustering hiérarchique des réplicats par site")
legend("topright", legend = names(site_colors), fill = site_colors, border = NA, bty = "n")

#=======================================================
# 9. Heatmap des fonctions les plus variables

library(pheatmap)
library(stringr)

# Définition top50 fonctions variables
var_fun <- apply(data_abund, 1, var)
top50 <- names(sort(var_fun, decreasing = TRUE))[1:50]

# Extraction sous-matrice et conversion en matrice
data_mat <- as.matrix(data_abund[top50, ])

# Conserver uniquement les lignes sans NA
data_mat <- data_mat[complete.cases(data_mat), ]

# Clustering sur lignes et colonnes
dist_rows <- vegdist(data_mat, method = "bray")
hc_rows <- hclust(dist_rows, method = "ward.D2")

#set.seed(123)
#pv_col <- pvclust(dist_bc, method.hclust = "ward.D2", nboot = 1000)

dist_cols <- vegdist(t(data_mat), method = "bray")
hc_cols <- hclust(dist_cols, method = "ward.D2")

###############################################
# Préparer annotation colonnes
samples <- colnames(data_mat)
site <- str_extract(samples, "^[A-Z0-9]+")
names(site) <- samples
annotation_col <- data.frame(Site = site)
rownames(annotation_col) <- samples

site_colors <- c(CAP = "gray20", NTE = "gray40", LL = "white", Y3 = "gray80")

# Plot heatmap
p <- pheatmap(data_mat,
         cluster_rows = hc_rows,
         cluster_cols = hc_cols,
         annotation_col = annotation_col,
         annotation_colors = list(Site = site_colors),
         fontsize_row = 8,
         fontsize_col = 8,
         cellheight = 10,   # hauteur des cellules
         cellwidth = 10,    # largeur des cellules
         main = "",
         annotation_legend = FALSE
)
print(p)

ggsave("heatmap_RelAb2-function.pdf", p, width = 9, height = 9, dpi = 300)
ggsave("heatmap-RelAbA2-function.png", p, width = 9, height = 9, dpi = 300, bg = "white")


# Ouvrir un périphérique graphique PNG
#png(
#  filename = "heatmap_function-gMethoxi-relAb.png",  # Nom du fichier de sortie
#  width = 22,                    # Largeur en pouces
#  height = 22,                    # Hauteur en pouces
#  units = "cm",                  # Unité : pouces
#  res = 300                      # Résolution (DPI)
#)

# Créer et afficher le heatmap
#pheatmap(
#  data_mat,
#  cluster_rows = hc_rows,
# cluster_cols = boot$hclust,
#  annotation_col = annotation_col,
#  annotation_colors = list(Site = site_colors),
#  fontsize_row = 8,
#  fontsize_col = 8,
#  cellheight = 10,
# cellwidth = 10,
# main = "",
# annotation_legend = FALSE
#)

# Fermer le périphérique graphique pour sauvegarder le fichier
#dev.off()

##############
#test statistique sur les fonctions les plus variables
###############
# ===================================================
# 10. Tests statistiques fonction par fonction (4 sites, 3 réplicats)
# ===================================================

library(tidyverse)
library(rstatix)
library(dplyr)

# Garder uniquement la sous‑matrice des 50 fonctions retenues par la heatmap
top50_mat <- as.matrix(data_abund[top50, ])          # 50 lignes, colonnes = échantillons
top50_mat <- top50_mat[complete.cases(top50_mat), ]  # enlever fonctions avec NA

# Convertir en long format (one row per fonction & sample)
top50_df <- as.data.frame(top50_mat) %>%
  tibble::rownames_to_column("Function") %>%
  pivot_longer(-Function, names_to = "Sample", values_to = "Abundance") %>%
  mutate(
    Site = str_extract(Sample, "^[A-Z0-9]+"),
    Site = factor(Site, levels = c("CAP", "NTE", "Y3", "LL"))
  )

# 10.1 Kruskal‑Wallis par fonction
kw_results <- top50_df %>%
  group_by(Function) %>%
  kruskal_test(Abundance ~ Site) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(p_adj)

# 10.2 Identifier fonctions significatives (p_adj < 0,05)
sig_functions <- kw_results %>%
  filter(p_adj < 0.05) %>%
  pull(Function)

# 10.3 Post‑hoc Dunn pour les fonctions significatives
dunn_results <- top50_df %>%
  filter(Function %in% sig_functions) %>%
  group_by(Function) %>%
  dunn_test(Abundance ~ Site, p.adjust.method = "BH") %>%
  ungroup()

# 10.4 Export des résultats
write.csv(kw_results, "kruskal_wallis_heatmap_top50-function-RelAb2.csv", row.names = FALSE)
write.csv(dunn_results, "dunn_posthoc_heatmap_top50-function-RelAb2.csv", row.names = FALSE)

# 10.5 Tableau de sortie joignable à la heatmap
functions_summary <- kw_results %>%
  select(Function, p, p_adj, n) %>%
  mutate(
    significant = ifelse(p_adj < 0.05, "yes", "no"),
    logarithm_p_adj = -log10(p_adj)
  )

write.csv(functions_summary, "functions_heatmap_stat_summary-RelAb2.csv", row.names = FALSE)

# 10.6 Affichage rapide
cat("Nombre de fonctions testées:", nrow(kw_results), "\n")
cat("Nombre de fonctions significatives (p_adj < 0.05) :", nrow(functions_summary %>% filter(significant == "yes")), "\n")


##############################################################
#heatmap toutes fonctions avec etoile sur les significatives
##############################################################
pval_vec <- functions_summary$p_adj
names(pval_vec) <- functions_summary$Function

lab_row <- rownames(data_mat)  

lab_row_star <- sapply(lab_row, function(f) {
  p <- pval_vec[f]
  if (is.na(p)) {
    f
  } else if (p < 0.001) {
    paste0(f, " ***")
  } else if (p < 0.01) {
    paste0(f, " **")
  } else if (p < 0.05) {
    paste0(f, " *")
  } else {
    f
  }
})

p2 <- pheatmap(data_mat,
              cluster_rows = hc_rows,
              cluster_cols = hc_cols,
              labels_row = lab_row_star,
              annotation_col = annotation_col,
              annotation_colors = list(Site = site_colors),
              fontsize_row = 8,
              fontsize_col = 8,
              cellheight = 10,   # hauteur des cellules
              cellwidth = 10,    # largeur des cellules
              main = "",
              annotation_legend = FALSE
)
print(p2)

ggsave("heatmap-RelAb2-function-Signif_star.pdf", p2, width = 9, height = 9, dpi = 300)
ggsave("heatmap-RelAb2-function-Signif_star.png", p2, width = 8, height = 9, dpi = 300, bg = "white")

# ================================================
# 11. HEATMAP des fonctions SIGNIFICATIVES (norm déjà faite)
# ================================================

library(tidyverse)
library(pheatmap)

# 11.1 Récupérer les fonctions significatives du test précédent
sig_functions_vec <- functions_summary %>%
  filter(significant == "yes") %>%
  pull(Function)

if (length(sig_functions_vec) == 0) {
 cat("Aucune fonction significative (p_adj < 0.05) → pas de heatmap spécifique.\n")
} else {
  # 11.2 Garder la matrice relAb (déjà normalisée par colonne = échantillons)
  # supposons que ta matrice relAb s'appelle `data_relAb` (créée avant)
  # par exemple par : data_relAb <- data_abund / colSums(data_abund)
  
  sig_mat_rel <- data_abund[sig_functions_vec, ]
  sig_mat_rel <- sig_mat_rel[complete.cases(sig_mat_rel), ]
  
  # Clustering sur lignes
  dist_rows <- vegdist(sig_mat_rel, method = "bray")
  hc_rows <- hclust(dist_rows, method = "ward.D2")
  
  
  # 11.3 Normalisation fonction par fonction (z-score) pour la heatmap, pas pour Bray‑Curtis
  #sig_mat_z <- t(scale(t(sig_mat_rel), center = TRUE, scale = TRUE))
  
  # 11.4 Préparation annotation Site
  samples <- colnames(sig_mat_rel)
  site <- str_extract(samples, "^[A-Z0-9]+")
  annotation_col <- data.frame(Site = factor(site))
  rownames(annotation_col) <- samples
  
  site_colors <- list(
    Site = c(
      CAP = "gray20",
      NTE = "gray40",
      LL  = "white",
      Y3  = "gray80"
    )
  )
  
  # 11.5 Clustering stable (sans z‑score pour la distance)
  # Bray‑Curtis sur les colonnes via la matrice CLR-modif (toujours > 0)
  #dist_cols <- vegdist(t(sig_mat_rel), method = "bray")
  #hc_cols   <- hclust(dist_cols, method = "ward.D2")
  
  # 11.6 Sauvegarde finale
  write.csv(sig_mat_rel, "heatmap_significants_function-Hell.csv", quote = FALSE)
 #write.csv(as.data.frame(sig_mat_z), "heatmap_significants_CLR_zscore.csv", quote = FALSE)
  
  p3 <- pheatmap::pheatmap(
    sig_mat_rel,
    cluster_rows = hc_rows,        # clustering par défaut (euclidien) sur z‑score
    cluster_cols = hc_cols,     # clustering par Bray‑Curtis sur relAb (> 0)
    annotation_col = annotation_col,
    annotation_colors = site_colors,
    fontsize_row = 10,
    fontsize_col = 8,
    cellheight = 10,
    cellwidth = 10,
    main = "",
    annotation_legend = FALSE
  )
  
  print(p3)
  
  ggsave("heatmap-RelAb2-function-significant.pdf", p3, width = 10, height = 10, dpi = 300)
  
  # 11.7 Heatmap (sans z‑score en ligne, clustering colonnes OK)
  #png("heatmap_functions_significatives_function-RelAbl.png", width = 1800, height = 2000, res = 200)
 # pheatmap::pheatmap(
 #   sig_mat_rel,
 #   cluster_rows = hc_rows,        # clustering par défaut (euclidien) sur z‑score
 #  cluster_cols = hc_cols,     # clustering par Bray‑Curtis sur relAb (> 0)
 #   annotation_col = annotation_col,
 #   annotation_colors = site_colors,
#    fontsize_row = 10,
 #   fontsize_col = 8,
 #   cellheight = 10,
 #   cellwidth = 10,
 #   main = "",
 #   annotation_legend = FALSE
 # )
 # dev.off()
}

