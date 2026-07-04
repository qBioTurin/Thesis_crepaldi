
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.23")

#options("install.lock"=FALSE)
#BiocManager::install("biobakery/maaslin3")
#BiocManager::install("remotes")


source("Libraries.R")
#QUI CAMBIA IL MAIN EL E I CONFUNDERS SE MI SERVE
# Create data for maaslin3
main_el= "tp"
confunders = c("sex", "bmi", "age_at_t0")
output_folder= "C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\risultati_maaslin_eukaryota"

#QUI CAMBIA IL PERCORSO DELL'RDS DI INTERESSE
rds = readRDS("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\shared_results_kraken_20260704_172753\\shared_results\\SUPERVISED_DEC\\Eukaryota_Supervised_decontam.rds")
otu <- rds@otu_table
taxa <- rds@tax_table
metadata_df <- as(rds@sam_data, "data.frame")

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

metadata_df$tp =as.factor(metadata_df$tp)
levels(metadata_df$tp) = c("T0", "T3", "T12", "T15", "T24", "36")
metadata_df$sex = as.factor(metadata_df$sex)
levels(metadata_df$sex) = c("M", "F")
metadata_df$lesioni =as.factor(metadata_df$lesioni)
levels(metadata_df$lesioni) = c("Yes", "No")
metadata_df$relapse =as.factor(metadata_df$relapse)
levels(metadata_df$relapse) = c("Yes", "No")

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}
# el_plot <- levels(metadata_df[[main_el]])[1]
el_plot2 <- levels(metadata_df[[main_el]])[2]
# main_el_plot <- paste(main_el, el_plot)
second_el_plot <- paste(main_el, el_plot2)
coef_plot_vars <- c(second_el_plot)

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
random_effects = NULL
group_effects = NULL
ordered_effects = NULL
strata_effects = NULL
heatmap_vars = NULL
Sys.setenv(OMP_NUM_THREADS=1, MKL_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1)
fit_out <- maaslin3(
  input_data = merged_df,
  input_metadata = metadata_df,
  output = output_folder,
  formula = formulaInUse_fin,
  fixed_effects = fx_effects,
  plot_summary_plot = TRUE,
  coef_plot_vars = coef_plot_vars,
  heatmap_vars = heatmap_vars,
  random_effects = random_effects,
  group_effects = group_effects,
  ordered_effects = ordered_effects,
  strata_effects = strata_effects,
  normalization = "TSS",
  transform = "LOG",
  augment = TRUE,
  standardize = TRUE,
  max_significance = 0.1,
  median_comparison_abundance = median_comparison_abundance_in,
  median_comparison_prevalence = FALSE,
  max_pngs = 30,
  summary_plot_first_n = 200,
  # pval_filter = 0.05,
  cores = 1
)

