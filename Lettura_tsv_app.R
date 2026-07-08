source("Libraries.R")

######risultati con pazienti totali#####
#batteri
all_definitive_bacteria = read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\risultati_maaslin_bacteria\\all_results.tsv")

significant_bacteria_prevalence <- all_definitive_bacteria %>%
  filter(model == "prevalence", metadata == "tp",
    pval_individual < 0.05 | qval_individual < 0.05)

significant_bacteria_abundance <- all_definitive_bacteria %>%
  filter(model == "abundance", metadata == "tp",
    pval_individual < 0.05 | qval_individual < 0.05)




#eukarioti
all_definitive_eukaryota = read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\risultati_maaslin_eukaryota\\all_results.tsv")

significant_eukaryota_prevalence <- all_definitive_eukaryota %>%
  filter(model == "prevalence", metadata == "tp",
         pval_individual < 0.05 | qval_individual < 0.05)

significant_eukaryota_abundance <- all_definitive_eukaryota %>%
  filter(model == "abundance", metadata == "tp",
         pval_individual < 0.05 | qval_individual < 0.05)

numero_taxa_eukariota = readRDS("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\shared_results\\SUPERVISED_DEC\\Eukaryota_Supervised_decontam.rds")
