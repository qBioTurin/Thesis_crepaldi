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


######risultati con pazienti totali
#batteri
all_definitive_bacteria = read.delim("")
significant_bacteria <- all_definitive_bacteria %>%
  filter(metadata == "tp") %>%
  filter(pval_individual < 0.05 | qval_individual < 0.05)


#archea
all_definitive_archea = read.delim("")
significant_archea <- all_definitive_archea %>%
  filter(metadata == "tp") %>%
  filter(pval_individual < 0.05 | qval_individual < 0.05)

#eukarioti
all_definitive_eukariota = read.delim("")
significant_eukariota <- all_definitive_eukariota %>%
  filter(metadata == "tp") %>%
  filter(pval_individual < 0.05 | qval_individual < 0.05)



