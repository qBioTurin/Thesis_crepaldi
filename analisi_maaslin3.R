
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.23")

#options("install.lock"=FALSE)
#BiocManager::install("biobakery/maaslin3")
#BiocManager::install("remotes")

###########MAASLIN3 X VARIAZIONE NEL TEMPO SOLO MAV
source("Libraries.R")
main_el= "tp"
confunders = c("sex", "bmi", "age_at_t0")
output_folder= "C:/Users/gialu/Desktop/tesi_magistrale/Analisi_pts_completi_app/Risultati_Maaslin/Risultati_MAV_app/risultati_maaslin_eukaryota"

#QUI CAMBIA IL PERCORSO DELL'RDS DI INTERESSE
rds = readRDS("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\SUPERVISED_DEC\\Eukaryota_Supervised_decontam.rds")
otu <- rds@otu_table
taxa <- rds@tax_table
metadata_df <- as(rds@sam_data, "data.frame")
unique(metadata_df$tp)


taxa <- as.data.frame(taxa)
taxa$Genus_Species <- paste(taxa$Genus, taxa$Species, sep = "_")
Tnew <- taxa %>% select(Genus_Species)
Tnew$taxid <- rownames(Tnew)
if (confunders %>% length() == 0) {
  formulaInUse <- main_el
} else {
  formulaInUse <- paste(main_el, paste(confunders, collapse = "+"), sep = "+")
}
otu <- as.data.frame(otu)
otu$taxid <- rownames(otu)
# Merge the two data frames by taxid
merged_df <- merge(otu, Tnew, by = "taxid")
rownames(merged_df) <- merged_df$Genus_Species
merged_df <- merged_df %>% select(-c(taxid, Genus_Species))
merged_df <- t(merged_df)

metadata_df$tp <- factor(
  as.character(metadata_df$tp),
  levels = c("T0", "T3", "T12", "T15", "T24", "T36"),
  ordered = FALSE
)
unique(metadata_df$tp)
metadata_df$sex = as.factor(metadata_df$sex)
levels(metadata_df$sex) = c("M", "F")
metadata_df$lesioni =as.factor(metadata_df$lesioni)
levels(metadata_df$lesioni) = c("Yes", "No")
metadata_df$relapse =as.factor(metadata_df$relapse)
levels(metadata_df$relapse) = c("Yes", "No")

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

coef_plot_vars <- c("tp T3", "tp T12", "tp T15", "tp T24", "tp T36")

discr_els_formatted <- unlist(lapply(confunders, function(d) {
  if (is.factor(metadata_df[[d]])) {
    paste(d, levels(metadata_df[[d]])) # tutti i livelli
  } else {
    d # continuous
  }
}))

##for stability of results across runs
formulaInUse_fin = paste0("~ ",formulaInUse, "+(1|original)")
set.seed(2204)
median_comparison_abundance_in = FALSE
fx_effects = NULL
random_effects =  NULL
group_effects = NULL
ordered_effects = NULL
strata_effects = NULL
Sys.setenv(OMP_NUM_THREADS=1, MKL_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1)
fit_out <- maaslin3(
  input_data = merged_df,
  input_metadata = metadata_df,
  output = output_folder,
  formula = formulaInUse_fin,
  fixed_effects = fx_effects,
  plot_summary_plot = TRUE,
  small_random_effects = TRUE,
  reference = c('tp,T0'),
  standardize= TRUE,#o continuous guardare se standardizzare o meno le variabili continue
  coef_plot_vars = coef_plot_vars,
  random_effects = random_effects,
  group_effects = group_effects,
  ordered_effects = ordered_effects,
  strata_effects = strata_effects,
  normalization = "TSS",
  heatmap_vars = c("tp T3", "tp T12", "tp T15", "tp T24", "tp T36", "sex F", "bmi", "age_at_t0"),
  transform = "LOG",
  augment = TRUE,
  max_significance = 0.1,
  median_comparison_abundance = median_comparison_abundance_in,
  median_comparison_prevalence = FALSE,
  max_pngs = 30,
  summary_plot_first_n = 200,
  # pval_filter = 0.05,
  cores = 1
)




########## MAASLIN3 PER HEALTHY VS MALATI SOLO GROUP ##########

main_el = "group"
confunders = c("sex", "bmi", "age_at_t0")
output_folder = "C:/Users/gialu/Desktop/tesi_magistrale/Analisi_pts_completi_app/Risultati_Maaslin/Risultati_app_healty_plus_MAV/Risultati_maaslin_complete_eukaryota"

# QUI CAMBIA IL PERCORSO DELL'RDS DI INTERESSE
rds = readRDS("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_app_healty_plus_MAV\\shared_results\\SUPERVISED_DEC\\Eukaryota_Supervised_decontam.rds")

otu <- rds@otu_table
taxa <- rds@tax_table
metadata_df <- as(rds@sam_data, "data.frame")
unique(metadata_df$group)

# tassonomia
taxa <- as.data.frame(taxa)
taxa$Genus_Species <- paste(taxa$Genus, taxa$Species, sep = "_")
Tnew <- taxa %>%
  select(Genus_Species)
Tnew$taxid <- rownames(Tnew)

# OTU
otu <- as.data.frame(otu)
otu$taxid <- rownames(otu)

# merge OTU + taxa
merged_df <- merge(otu, Tnew, by = "taxid")
rownames(merged_df) <- merged_df$Genus_Species
merged_df <- merged_df %>%
  select(-c(taxid, Genus_Species))

# trasposizione (righe = campioni, colonne = taxa)
merged_df <- t(merged_df)

# GROUP come factor
metadata_df$group <- factor(
  as.character(metadata_df$group),
  levels = c("Healthy", "Affected"),
  ordered = FALSE
)
unique(metadata_df$group)

# sex factor
metadata_df$sex = as.factor(metadata_df$sex)
levels(metadata_df$sex) = c("M", "F")

# formula
if (confunders %>% length() == 0) {
  formulaInUse <- main_el
} else {
  formulaInUse <- paste(
    main_el,
    paste(confunders, collapse = "+"),
    sep = "+"
  )
}

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

coef_plot_vars <- c(
  "group Affected"
)

discr_els_formatted <- unlist(
  lapply(confunders, function(d) {
    if (is.factor(metadata_df[[d]])) {
      paste(d, levels(metadata_df[[d]]))
    } else {
      d
    }
  })
)

#formula finale (nessun random effect: qui non c'e' ripetizione di soggetto
formulaInUse_fin = paste0("~ ", formulaInUse)

set.seed(2204)
median_comparison_abundance_in = FALSE
fx_effects = NULL
random_effects = NULL
group_effects = NULL
ordered_effects = NULL
strata_effects = NULL

Sys.setenv(
  OMP_NUM_THREADS=1,
  MKL_NUM_THREADS=1,
  OPENBLAS_NUM_THREADS=1
)

fit_out <- maaslin3(
  input_data = merged_df,
  input_metadata = metadata_df,
  output = output_folder,
  
  formula = formulaInUse_fin,
  
  fixed_effects = fx_effects,
  
  plot_summary_plot = TRUE,
  reference = c("group,Healthy"),
  
  standardize = TRUE,
  
  coef_plot_vars = coef_plot_vars,
  
  random_effects = random_effects,
  group_effects = group_effects,
  ordered_effects = ordered_effects,
  strata_effects = strata_effects,
  
  normalization = "TSS",
  
  heatmap_vars = c(
    "group Affected",
    "sex F",
    "bmi",
    "age_at_t0"
  ),
  
  transform = "LOG",
  
  augment = TRUE,
  
  max_significance = 0.1,
  
  median_comparison_abundance = median_comparison_abundance_in,
  
  median_comparison_prevalence = FALSE,
  
  max_pngs = 30,
  
  summary_plot_first_n = 200,
  
  cores = 1
)





