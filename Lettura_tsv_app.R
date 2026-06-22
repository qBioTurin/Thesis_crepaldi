source("Libraries.R")

#all.result.tsv del 05/06
dati <- read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\time_MAV\\MAASLIN3\\tp_plus_sex_age_at_t0_bmi_2026-06-05\\all_results.tsv")

#sigificant.result del 05/06
signifcant_dati <- read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\time_MAV\\MAASLIN3\\tp_plus_sex_age_at_t0_bmi_2026-06-05\\significant_results.tsv")

str(dati)

dati_filtrati <- dati %>%
  filter(metadata == "tp") %>%
  filter(pval_joint < 0.05 | qval_joint < 0.05)
