source("Libraries.R")

######risultati con dati solo 7 pazienti####
#all.result.tsv del 23/06
all_7_pts <- read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\time_MAV\\MAASLIN3\\tp_plus_sex_age_at_t0_bmi_2026-06-24\\all_results.tsv")
#sigificant.result del 23/06
signifcant_7_pts <- read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\time_MAV\\MAASLIN3\\tp_plus_sex_age_at_t0_bmi_2026-06-24\\significant_results.tsv")

#FILTRO 
significant_result_tp <- dati %>%
  filter(metadata == "tp") %>%
  filter(pval_individual < 0.05 | qval_individual < 0.05)



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




