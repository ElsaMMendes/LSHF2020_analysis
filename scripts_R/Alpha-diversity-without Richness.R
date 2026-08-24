setwd("~/Documents/Documents Alain/ARTICLE/Metagenomic LSHF avril 2024/Papier Mendes et al mars 2025/tableau a faire")
install.packages("pgirmess")
install.packages("PMCMRplus")
install.packages("ggpubr")

# Charger les bibliothèques nécessaires
library(readxl)
library(dplyr)
library(ggplot2)
library(vegan)
library(tidyr)
library(patchwork)
library(rstatix)
library(car)
library(pgirmess)
library(PMCMRplus)
library(FSA)
library(broom)
library(ggpubr)
library(readr)
library(cowplot)

# 1. Lire le fichier Excel
# Remplacez "votre_fichier.xlsx" par le chemin vers votre fichier
data <- read_excel("Diversite-alpha-beta-MAGs.xlsx" , sheet = "per replicat")

# Afficher les noms des colonnes pour vérifier
print(colnames(data))

# 2. Préparer les données
abundance_columns <- c("rep_1", "rep_2", "rep_3")  # Remplacez par les noms exacts de vos colonnes

# Vérifier si les colonnes existent dans le jeu de données
if(!all(abundance_columns %in% colnames(data))) {
  stop("Une ou plusieurs colonnes d'abondance ne sont pas trouvées. Vérifiez les noms des colonnes.")
}

# Convertir les colonnes d'abondance en numérique
data <- data %>%
  mutate(across(all_of(abundance_columns), as.numeric))

# 3. Calculer les indices de diversité alpha par réplicat et par site
calculate_alpha_diversity <- function(abundance_vector) {
  if (length(abundance_vector) == 0 || all(is.na(abundance_vector))) {
    return(data.frame(shannon = NA, simpson = NA))
  }
  proportions <- abundance_vector / sum(abundance_vector, na.rm = TRUE)
  if (any(is.nan(proportions))) {
    return(data.frame(shannon = NA, simpson = NA))
  }
  shannon <- diversity(proportions, index = "shannon")
  simpson <- diversity(proportions, index = "simpson")
  #richness <- length(abundance_vector)
  
  return(data.frame(shannon = shannon, simpson = simpson))
}

# Préparer les données pour chaque réplicat
diversity_by_replicate <- data %>%
  pivot_longer(
    cols = all_of(abundance_columns),
    names_to = "replicate",
    values_to = "abundance"
  ) %>%
  group_by(site, replicate) %>%
  summarise(
    diversity = list(calculate_alpha_diversity(abundance))
  ) %>%
  unnest(cols = c(diversity)) %>%
  mutate(replicate = factor(replicate, levels = abundance_columns))

# 4. Calculer les paramètres des boxplots pour chaque indice et par site
boxplot_params <- diversity_by_replicate %>%
  group_by(site, replicate) %>%
  summarise(
    # Pour Shannon
    median_shannon = median(shannon, na.rm = TRUE),
    q1_shannon = quantile(shannon, 0.25, na.rm = TRUE),
    q3_shannon = quantile(shannon, 0.75, na.rm = TRUE),
    min_shannon = min(shannon, na.rm = TRUE),
    max_shannon = max(shannon, na.rm = TRUE),
    iqr_shannon = IQR(shannon, na.rm = TRUE),
    lower_whisker_shannon = q1_shannon - 1.5 * iqr_shannon,
    upper_whisker_shannon = q3_shannon + 1.5 * iqr_shannon,
    
    # Pour Simpson
    median_simpson = median(simpson, na.rm = TRUE),
    q1_simpson = quantile(simpson, 0.25, na.rm = TRUE),
    q3_simpson = quantile(simpson, 0.75, na.rm = TRUE),
    min_simpson = min(simpson, na.rm = TRUE),
    max_simpson = max(simpson, na.rm = TRUE),
    iqr_simpson = IQR(simpson, na.rm = TRUE),
    lower_whisker_simpson = q1_simpson - 1.5 * iqr_simpson,
    upper_whisker_simpson = q3_simpson + 1.5 * iqr_simpson,
    
    # Pour Richness
   # median_richness = median(richness, na.rm = TRUE),
   # q1_richness = quantile(richness, 0.25, na.rm = TRUE),
   # q3_richness = quantile(richness, 0.75, na.rm = TRUE),
   # min_richness = min(richness, na.rm = TRUE),
   # max_richness = max(richness, na.rm = TRUE),
   # iqr_richness = IQR(richness, na.rm = TRUE),
   # lower_whisker_richness = q1_richness - 1.5 * iqr_richness,
   # upper_whisker_richness = q3_richness + 1.5 * iqr_richness
  )

# Sauvegarder les paramètres des boxplots dans un fichier CSV
write_csv(boxplot_params, "boxplot_diversity_indices_parameters_per_replicate-without Richness.csv")

#Extraire les indice de boxplot_params
colonnes_souhaitees <- c("site", "replicate", "median_shannon", "median_simpson")
diversity_indices_per_replicate <- boxplot_params[, colonnes_souhaitees, drop = FALSE]
write.csv(diversity_indices_per_replicate, "diversity_indices_per_replicat-without Richness.csv", row.names = FALSE)

### Visualisation. avec indicateurs stat

# 1. Lire le fichier CSV contenant les indices de diversité
diversity_data <- read_csv("diversity_indices_per_replicate.csv")

# Vérifier que les colonnes existent
required_columns <- c("site", "shannon", "simpson")
if (!all(required_columns %in% colnames(diversity_data))) {
  stop("Certaines colonnes nécessaires sont manquantes dans le fichier CSV. Vérifiez les noms des colonnes.")
}

# Ouvrir une connexion vers un fichier texte
file_path <- "resultats_alpha_div.txt"
sink(file_path)  # Redirige toute la sortie vers le fichier


# Fonction pour effectuer les analyses statistiques pour un indice donné
analyze_index <- function(data, index_name) {
  cat(paste("\nAnalyse pour l'indice", index_name, "\n"))
  
  # Test de normalité
  normality_test <- shapiro.test(data[[index_name]])
  cat("Test de normalité (Shapiro-Wilk):\n")
  print(normality_test)
  
  # Test d'homogénéité des variances
  levene_test_result <- try(car::leveneTest(data[[index_name]], data$site), silent = TRUE)
  if (inherits(levene_test_result, "try-error")) {
    cat("\nLevene's test could not be computed properly.\n")
    levene_test_result <- list(`Pr(>F)` = c(NA))
  } else {
    cat("\nTest d'homogénéité des variances (Levene):\n")
    print(levene_test_result)
  }
  
  # Décider du test à utiliser
  if (normality_test$p.value > 0.05 & (!is.na(levene_test_result$`Pr(>F)`[1]) & levene_test_result$`Pr(>F)`[1] > 0.05)) {
    cat("\nLes hypothèses de normalité et d'homogénéité des variances sont respectées. Utilisation de l'ANOVA.\n")
    
    # ANOVA
    aov_result <- aov(as.formula(paste0(index_name, " ~ site")), data = data)
    cat("\nRésultats de l'ANOVA:\n")
    print(summary(aov_result))
    
    # Test post-hoc de Tukey
    tukey_result <- TukeyHSD(aov_result)
    cat("\nRésultats du test post-hoc de Tukey:\n")
    print(tukey_result)
    
    # Convertir les résultats de Tukey en data frame
    tukey_df <- as.data.frame(tukey_result$site)
    tukey_df$Comparison <- rownames(tukey_df)
    rownames(tukey_df) <- NULL
    
    return(list(method = "ANOVA", test_result = aov_result, posthoc = tukey_df))
  } else {
    cat("\nLes hypothèses ne sont pas respectées. Utilisation du test de Kruskal-Wallis.\n")
    
    # Test de Kruskal-Wallis
    kruskal_result <- kruskal.test(as.formula(paste0(index_name, " ~ site")), data = data)
    cat("\nRésultats du test de Kruskal-Wallis:\n")
    print(kruskal_result)
    
    # Test post-hoc de Dunn
    dunn_result <- FSA::dunnTest(x = data[[index_name]], g = data$site, method = "bonferroni")
    cat("\nRésultats du test post-hoc de Dunn:\n")
    print(dunn_result)
    
    # Convertir les résultats de Dunn en data frame
    dunn_df <- as.data.frame(dunn_result$res)
    
    return(list(method = "Kruskal-Wallis", test_result = kruskal_result, posthoc = dunn_df))
  }
}

# Fermer la connexion pour revenir à la console
sink()  # Arrête la redirection


# Fonction pour créer un boxplot en niveaux de gris avec annotations de significativité
create_grayscale_boxplot <- function(data, index_name, posthoc_results, test_type) {
  # Définir le nom complet de l'indice pour le titre
  index_full_name <- switch(index_name,
                            "shannon" = "Shannon",
                            "simpson" = "Simpson")
  
  p <- ggplot(data, aes(x = site, y = .data[[index_name]], fill = site)) +
    geom_boxplot() +
    labs(x = "Site", y = "Value", title = "") +  # Titre = nom de l'indice, ordonnée = "Value"
    scale_fill_grey(start = 0.75, end = 0.75) +  # Échelle de gris
    theme_minimal(base_family = "") +
    theme(
      legend.position = "none",
      axis.title = element_text(face = "bold"),
      panel.background = element_rect(fill = "white", color = "grey50"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      #plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),  # Titre centré et en gras
        axis.text.x = element_text(size =10),   # Taille des valeurs sur l'axe x
        axis.text.y = element_text(size = 10),   # Taille des valeurs sur l'axe y
        axis.title.x = element_text(size = 12), # Taille de la légende de l'axe x
        axis.title.y = element_text(size = 12)  # Taille de la légende de l'axe y
    )
  
  # Ajouter les annotations de significativité
  if (test_type == "ANOVA" && nrow(posthoc_results) > 0 && "p adj" %in% colnames(posthoc_results)) {
    sig_pairs <- posthoc_results[posthoc_results$`p adj` < 0.05, ]
    if (nrow(sig_pairs) > 0) {
      for (i in 1:nrow(sig_pairs)) {
        pair <- strsplit(as.character(sig_pairs$Comparison[i]), "-")[[1]]
        if (length(pair) == 2) {
          p_value <- sig_pairs$`p adj`[i]
          p <- p + geom_signif(
            comparisons = list(c(pair[1], pair[2])),
            annotations = paste0("p = ", format(p_value, digits = 3)),
            y_position = max(data[[index_name]]) + 0.1 * i,
            tip_length = 0.01,
            color = "black",
            linecolor = "black"
          )
        }
      }
    }
  } else if (test_type == "Kruskal-Wallis" && nrow(posthoc_results) > 0 && "P.adj" %in% colnames(posthoc_results)) {
    sig_pairs <- posthoc_results[posthoc_results$P.adj < 0.05, ]
    if (nrow(sig_pairs) > 0) {
      for (i in 1:nrow(sig_pairs)) {
        pair <- strsplit(as.character(sig_pairs$Comparison[i]), " - ")[[1]]
        if (length(pair) == 2) {
          p_value <- sig_pairs$P.adj[i]
          p <- p + geom_signif(
            comparisons = list(c(pair[1], pair[2])),
            annotations = paste0("P.adj = ", format(p_value, digits = 3)),
            y_position = max(data[[index_name]]) + 0.01 * i,
            tip_length = 0.01,
            color = "black",
            linecolor = "black"
          )
        }
      }
    }
  }
  
  
  return(p)
}


# Effectuer les analyses pour chaque indice
indices <- c("shannon", "simpson")
results <- list()
plots <- list()

for (index in indices) {
  cat(paste("\n\nAnalyse pour l'indice", index, "\n"))
  results[[index]] <- analyze_index(diversity_data, index)
  
  # Créer le boxplot en niveaux de gris avec annotations
  plots[[index]] <- create_grayscale_boxplot(
    diversity_data, index,
    results[[index]]$posthoc,
    results[[index]]$method
  )
  
  # Sauvegarder les résultats
  if (results[[index]]$method == "ANOVA") {
    write_csv(broom::tidy(results[[index]]$test_result), paste0(index, "_anova_results-bis.csv"))
    write_csv(results[[index]]$posthoc, paste0(index, "_posthoc_results-bis.csv"))
  } else {
    write_csv(broom::tidy(results[[index]]$test_result), paste0(index, "_kruskal_results-bis.csv"))
    write_csv(results[[index]]$posthoc, paste0(index, "_dunn_results-bis.csv"))
  }
}

# Combiner les 2 boxplots en un seul graphique (côte à côte) SANS TITRE GLOBAL
combined_plot <- gridExtra::grid.arrange(
  grobs = plots,
  ncol = 2  # 3 colonnes pour les 3 graphiques côte à côte
)

# Ajouter lettres "A" et "B" en dehors du cadre
p_with_letters <- ggdraw(combined_plot) +
  draw_label("A", x = 0.05, y = 0.97, size = 12, fontface = "bold") +
  draw_label("B", x = 0.57, y = 0.97, size = 12, fontface = "bold")

print(p_with_letters)

# Sauvegarder le graphique combiné
ggsave("Alpha_div_stat-without Richness.pdf", p_with_letters, width = 7, height = 6, dpi = 300)
ggsave("Alpha_div_stat-without Richness.png", p_with_letters, width = 7, height = 6, dpi = 300)
