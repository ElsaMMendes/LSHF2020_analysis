#SET directory
setwd("~/Documents/Documents Alain/ARTICLE/Metagenomic LSHF avril 2024/Papier Mendes et al mars 2025/tableau a faire")


# ================================
# 1. PACKAGES
# ================================
library(readxl)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(ggrepel)
library(tibble)
library(grid)  # pour les flèches


#Fichiers normalisés pour db-RDA
env <- read_excel("chimie-envir-SeaWonly-zscore.xlsx")

# 1. Lire le fichier Excel contenant les données
fun <- read_excel("Diversite-alpha-beta-MAGs.xlsx", sheet = "per replicat")

# 2. Préparer les données d'abondance
for (col in c("rep_1", "rep_2", "rep_3")) {
  fun[[col]] <- gsub(",", ".", as.character(fun[[col]]))
  fun[[col]] <- suppressWarnings(as.numeric(fun[[col]]))
  fun[[col]][is.na(fun[[col]])] <- 0
}

# 3. Créer une matrice de comptage large (échantillons en colonnes, MAGs en lignes)
count_fun_long <- fun %>%
  select(site, MAG, starts_with("rep_")) %>%
  pivot_longer(cols = starts_with("rep_"),
               names_to = "replicate",
               values_to = "abundance") %>%
  mutate(sample_id = paste(site, replicate, sep = "_")) %>%
  select(MAG, sample_id, abundance)

# 4. Ajouter les informations taxonomiques
taxonomic_fun <- fun %>%
  select(MAG, family) %>%
  drop_na(family) %>%
  distinct(MAG, .keep_all = TRUE)

# 5. Agrégation des abondances par family
family_abundance <- count_fun_long %>%
  left_join(taxonomic_fun, by = "MAG") %>%
  group_by(sample_id, family) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0)

# 6. Créer une matrice d'abondance avec les family en lignes et les échantillons en colonnes
family_abundance_matrix <- family_abundance %>%
  column_to_rownames(var = "family") %>%
  as.matrix()

# 7. Normalisation
fun_mat <- family_abundance_matrix / colSums(family_abundance_matrix)

fun_mat <- t(fun_mat)

# ================================
# 3. NETTOYAGE DES NOMS
# ================================
#clean_names <- function(x) {
#  x %>%
#    gsub("_", "-", .) %>%
#    gsub("rep", "Rep", ., ignore.case = TRUE) %>%
#    gsub("CAPp", "CAP", .) %>%
#    trimws()
#}

#env$Site <- clean_names(env$Site)
#colnames(fun)[-1] <- clean_names(colnames(fun)[-1])

# ================================
# 4. MATRICE DES FONCTIONS
# ================================
fun_mat <- family_abundance %>%
  column_to_rownames(var = colnames(family_abundance)[1]) %>%
  t() %>%
  as.data.frame()

fun_mat$Site <- rownames(fun_mat)

# ================================
# 5. FUSION DES DONNÉES
# ================================
data_merged <- inner_join(env, fun_mat, by = "Site")

# ================================
# 6. VARIABLES ENVIRONNEMENTALES
# ⚠️ ADAPTE ICI SI BESOIN
# ================================
env_vars <- data_merged %>%
  select(
#    pH_diffus,
#    Eh_diffus,
#    Fe_diffus,
#    H2S_diffus,
    pH_seawater,
    Eh_seawater,
    Fe_seawater,
    Fe2O3,
  )
# variable qualitative substratum (facteur)
substratum_factor <- factor(data_merged$substratum)

# garder lignes complètes
complete_rows <- complete.cases(env_vars)& !is.na(substratum_factor)
env_vars <- env_vars[complete_rows, ]
substratum_factor <- substratum_factor[complete_rows]


# ================================
# 7. MATRICE DES FONCTIONS
# ================================
fun_vars <- data_merged %>%
  select(colnames(fun_mat)[colnames(fun_mat) != "Site"])

fun_vars <- fun_vars[complete_rows, ]

# ================================
# 8. DISTANCE + DB-RDA
# ================================
dist_fun <- vegdist(fun_vars[complete_rows, ], method = "bray")

# On ajoute substratum_factor comme variable explicative qualitative dans le modèle,
# via la formule en R avec + substratum_factor (par nom variable)
dbrda_model <- capscale(dist_fun ~ . + substratum_factor, data = as.data.frame(env_vars))


# ================================
# 9. TESTS STAT
# ================================
print(anova(dbrda_model))
print(anova(dbrda_model, by = "terms"))

# ================================
# 10. SCORES POUR GGPLOT
# ================================
# Sites
# Supposons que complete_rows est déjà défini (lignes complètes)
scores_sites <- scores(dbrda_model, display = "sites") %>%
  as.data.frame()

scores_sites$Site <- fun_mat$Site[complete_rows]


# Groupe (CAP, NTE, LL, Y3)
scores_sites$Site_group <- gsub("_rep_[0-9]+", "", scores_sites$Site)
print(unique(scores_sites$Site_group))


# Variables environnementales UNIQUEMENT
scores_env <- scores(dbrda_model, display = "bp") %>%
  as.data.frame()

scores_env$Variable <- rownames(scores_env)

# sécurité anti-bug
scores_env <- scores_env %>%
  filter(Variable %in% colnames(env_vars))

# ================================
# 11. SCALING DES FLÈCHES
# ================================
scores_env <- scores_env %>%
  mutate(across(where(is.numeric), ~ . * 0.5))


# ================================
# 12. VARIANCE EXPLIQUÉE
# ================================
eig_vals <- eigenvals(dbrda_model)
var_explained <- eig_vals / sum(eig_vals)

xlab <- paste0("db-RDA1 (", round(var_explained[1] * 100, 1), "%)")
ylab <- paste0("db-RDA2 (", round(var_explained[2] * 100, 1), "%)")

# ================================
# Représentation de substratum (centroïdes)
# ================================

# Calcul des centroïdes des sites par substratum
centroids <- aggregate(cbind(CAP1 = scores_sites$CAP1, CAP2 = scores_sites$CAP2),
                       by = list(substratum = substratum_factor),
                       FUN = mean)


# ================================
# 13. PLOT GGPLOT
# ================================

# Définir les symboles pour chaque site
symbols <- c("CAP" = 15, "NTE" = 16, "Y3" = 18, "LL" = 17)  # 15 = carré, 16 = cercle, 17 = triangle, 18 = losange

# Créer le graphique ggplot avec les top 5 variables discriminantes

dbRDA_p1 <- ggplot() +
  geom_point(data = scores_sites, aes(x = CAP1, y = CAP2, shape = Site_group), size = 3, color = "black") +  # Points en noir
  geom_segment(data = scores_env,
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length = unit(0.1, "cm")),
               color = "black",  # Flèches en gris
               linewidth = 0.2) +
  geom_text(data = scores_env,
                  aes(x = CAP1, y = CAP2, label = Variable),
                  color = "black",  # Texte en gris
                  size = 4,
            nudge_x = -0.10) +
  geom_point(data = centroids,
             aes(x = CAP1, y = CAP2),
             color = "black",
             size = 4, shape = 25, fill = "gray80") +  # losanges pour centroïdes substratum
  geom_text_repel(data = centroids,
                  aes(x = CAP1, y = CAP2, label = substratum),
                  color = "black",
                #fontface = "bold",
                 size = 4,
                nudge_y = 0.10,
                nudge_x =-0.10) +
  scale_shape_manual(values = symbols) +  # Utiliser les symboles définis
  scale_color_manual(values = rep("gray50", length(unique(scores_sites$Site_group)))) +  # Échelle de couleur pour les ellipses
  ggforce::geom_mark_ellipse(data = scores_sites, aes(x = CAP1, y = CAP2, group = Site_group), color = "gray50", linewidth = 0.5) +  # Ellipses en gris
  theme_minimal(base_size = 7) +
  labs(title = "",
       x = xlab,
       y = ylab,
       shape = "Site") +
  #theme(plot.title = element_text(hjust = 0.5, face = "bold"),
       # legend.position = "none")

theme(
  legend.position = "none",
  panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
  axis.text.x = element_text(size =10),   # Taille des valeurs sur l'axe x
  axis.text.y = element_text(size = 10),   # Taille des valeurs sur l'axe y
  axis.title.x = element_text(size = 12), # Taille de la légende de l'axe x
  axis.title.y = element_text(size = 12))  # Taille de la légende de l'axe y# Supprime la légende


# Ajouter le "A" en dehors du cadre
p_with_letterA <- ggdraw(dbRDA_p1) +
  draw_label("A", x = 0.05, y = 0.97, size = 12, fontface = "bold")


#afficher plot
print(p_with_letterA)

# Sauvegarder le graphique
ggsave("db-RDA-SeaWonly-substratum-family.pdf", p_with_letterA, width = 8, height = 6, dpi = 300)
ggsave("db-RDA-SeaWonly-substratum-family.png", p_with_letterA, width = 8, height = 6, dpi = 300)

##################DIFFUS

#Fichiers normalisés pour db-RDA
env <- read_excel("chimie-envir-ter-difonly-zscore-dbRDA.xlsx")

# 1. Lire le fichier Excel contenant les données
fun <- read_excel("Diversite-alpha-beta-MAGs-withoutLL.xlsx", sheet = "per replicat")

# 2. Préparer les données d'abondance
for (col in c("rep_1", "rep_2", "rep_3")) {
  fun[[col]] <- gsub(",", ".", as.character(fun[[col]]))
  fun[[col]] <- suppressWarnings(as.numeric(fun[[col]]))
  fun[[col]][is.na(fun[[col]])] <- 0
}

# 3. Créer une matrice de comptage large (échantillons en colonnes, MAGs en lignes)
count_fun_long <- fun %>%
  select(site, MAG, starts_with("rep_")) %>%
  pivot_longer(cols = starts_with("rep_"),
               names_to = "replicate",
               values_to = "abundance") %>%
  mutate(sample_id = paste(site, replicate, sep = "_")) %>%
  select(MAG, sample_id, abundance)

# 4. Ajouter les informations taxonomiques
taxonomic_fun <- fun %>%
  select(MAG, family) %>%
  drop_na(family) %>%
  distinct(MAG, .keep_all = TRUE)

# 5. Agrégation des abondances par family
family_abundance <- count_fun_long %>%
  left_join(taxonomic_fun, by = "MAG") %>%
  group_by(sample_id, family) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0)

# 6. Créer une matrice d'abondance avec les family en lignes et les échantillons en colonnes
family_abundance_matrix <- family_abundance %>%
  column_to_rownames(var = "family") %>%
  as.matrix()

# 7. Normalisation
fun_mat <- family_abundance_matrix / colSums(family_abundance_matrix)

fun_mat <- t(fun_mat)

# ================================
# 3. NETTOYAGE DES NOMS
# ================================
#clean_names <- function(x) {
#  x %>%
#    gsub("_", "-", .) %>%
#    gsub("rep", "Rep", ., ignore.case = TRUE) %>%
#    gsub("CAPp", "CAP", .) %>%
#    trimws()
#}

#env$Site <- clean_names(env$Site)
#colnames(fun)[-1] <- clean_names(colnames(fun)[-1])

# ================================
# 4. MATRICE DES FONCTIONS
# ================================
fun_mat <- family_abundance %>%
  column_to_rownames(var = colnames(family_abundance)[1]) %>%
  t() %>%
  as.data.frame()

fun_mat$Site <- rownames(fun_mat)

# ================================
# 5. FUSION DES DONNÉES
# ================================
data_merged <- inner_join(env, fun_mat, by = "Site")

# ================================
# 6. VARIABLES ENVIRONNEMENTALES
# ⚠️ ADAPTE ICI SI BESOIN
# ================================
env_vars <- data_merged %>%
  select(
    pH_diffus,
   Eh_diffus,
    Fe_diffus,
    #    H2S_diffus,
    # pH_seawater,
    # Eh_seawater,
    # Fe_seawater,
    Fe2O3,
  )
# variable qualitative substratum (facteur)
substratum_factor <- factor(data_merged$substratum)

# garder lignes complètes
complete_rows <- complete.cases(env_vars)& !is.na(substratum_factor)
env_vars <- env_vars[complete_rows, ]
substratum_factor <- substratum_factor[complete_rows]


# ================================
# 7. MATRICE DES FONCTIONS
# ================================
fun_vars <- data_merged %>%
  select(colnames(fun_mat)[colnames(fun_mat) != "Site"])

fun_vars <- fun_vars[complete_rows, ]

# ================================
# 8. DISTANCE + DB-RDA
# ================================
dist_fun <- vegdist(fun_vars[complete_rows, ], method = "bray")

# On ajoute substratum_factor comme variable explicative qualitative dans le modèle,
# via la formule en R avec + substratum_factor (par nom variable)
dbrda_model <- capscale(dist_fun ~ . + substratum_factor, data = as.data.frame(env_vars))


# ================================
# 9. TESTS STAT
# ================================
print(anova(dbrda_model))
print(anova(dbrda_model, by = "terms"))

# ================================
# 10. SCORES POUR GGPLOT
# ================================
# Sites
# Supposons que complete_rows est déjà défini (lignes complètes)
scores_sites <- scores(dbrda_model, display = "sites") %>%
  as.data.frame()

scores_sites$Site <- fun_mat$Site[complete_rows]


# Groupe (CAP, NTE, LL, Y3)
scores_sites$Site_group <- gsub("_rep_[0-9]+", "", scores_sites$Site)
print(unique(scores_sites$Site_group))


# Variables environnementales UNIQUEMENT
scores_env <- scores(dbrda_model, display = "bp") %>%
  as.data.frame()

scores_env$Variable <- rownames(scores_env)

# sécurité anti-bug
scores_env <- scores_env %>%
  filter(Variable %in% colnames(env_vars))

# ================================
# 11. SCALING DES FLÈCHES
# ================================
scores_env <- scores_env %>%
  mutate(across(where(is.numeric), ~ . * 0.5))


# ================================
# 12. VARIANCE EXPLIQUÉE
# ================================
eig_vals <- eigenvals(dbrda_model)
var_explained <- eig_vals / sum(eig_vals)

xlab <- paste0("db-RDA1 (", round(var_explained[1] * 100, 1), "%)")
ylab <- paste0("db-RDA2 (", round(var_explained[2] * 100, 1), "%)")

# ================================
# Représentation de substratum (centroïdes)
# ================================

# Calcul des centroïdes des sites par substratum
centroids <- aggregate(cbind(CAP1 = scores_sites$CAP1, CAP2 = scores_sites$CAP2),
                       by = list(substratum = substratum_factor),
                       FUN = mean)


# ================================
# 13. PLOT GGPLOT
# ================================

# Définir les symboles pour chaque site
symbols <- c("CAP" = 15, "NTE" = 16, "Y3" = 18, "LL" = 17)  # 15 = carré, 16 = cercle, 17 = triangle, 18 = losange

# Créer le graphique ggplot avec les top 5 variables discriminantes

dbRDA_p2 <- ggplot() +
  geom_point(data = scores_sites, aes(x = CAP1, y = CAP2, shape = Site_group), size = 3, color = "black") +  # Points en noir
  geom_segment(data = scores_env,
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length = unit(0.1, "cm")),
               color = "black",  # Flèches en gris
               linewidth = 0.4) +
  geom_text_repel(data = scores_env,
                  aes(x = CAP1, y = CAP2, label = Variable),
                  color = "black",  # Texte en gris
                  size = 4) +
  geom_point(data = centroids,
             aes(x = CAP1, y = CAP2), 
             color = "black",
             size = 4, shape = 25, fill = "gray80") +  # losanges pour centroïdes substratum
  geom_text_repel(data = centroids,
                  aes(x = CAP1, y = CAP2, label = substratum), 
                  color = "black",
                  #fontface = "",
                  size = 4,
                  nudge_y = 0.10, 
                  nudge_x = 0.10) +
  scale_shape_manual(values = symbols) +  # Utiliser les symboles définis
  scale_color_manual(values = rep("gray50", length(unique(scores_sites$Site_group)))) +  # Échelle de couleur pour les ellipses
  ggforce::geom_mark_ellipse(data = scores_sites, aes(x = CAP1, y = CAP2, group = Site_group), color = "gray50", linewidth = 0.5) +  # Ellipses en gris
  theme_minimal(base_size = 7) +
  labs(title = "",
       x = xlab,
       y = ylab,
       shape = "Site") +
  #theme(plot.title = element_text(hjust = 0.5, face = "bold"),
  # legend.position = "none")
  
  theme(
    legend.position = "none",
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    axis.text.x = element_text(size =10),   # Taille des valeurs sur l'axe x
    axis.text.y = element_text(size = 10),   # Taille des valeurs sur l'axe y
    axis.title.x = element_text(size = 12), # Taille de la légende de l'axe x
    axis.title.y = element_text(size = 12))  # Taille de la légende de l'axe y# Supprime la légende

# Ajouter le "A" en dehors du cadre
p_with_letterB <- ggdraw(dbRDA_p2) +
  draw_label("B", x = 0.05, y = 0.97, size = 12, fontface = "bold")


#afficher plot
print(p_with_letterB)


#afficher plot
print(p_with_letterB)

# Sauvegarder le graphique
ggsave("db-RDA-Difonly-substratum-family.pdf", p_with_letterB, width = 8, height = 6, dpi = 300)
ggsave("db-RDA-Difonly-substratum-family.png", p_with_letterB, width = 8, height = 6, dpi = 300)


#=============================
#combined PLOT Figure 1
#============================
# Combiner les deux graphiques
combined_plot <- p_with_letterA / p_with_letterB

# Afficher le résultat
print(combined_plot)

# Sauvegarder en TIFF
ggsave("Figure4.TIFF", plot = combined_plot, device = "tiff", dpi = 300, width = 20, height = 29, units = "cm")

ggsave("Figure4.png", plot = combined_plot, device = "png", dpi = 300, width = 20, height = 29, units = "cm")



