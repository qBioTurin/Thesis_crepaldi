source("Libraries.R")

#all.result.tsv del 23/06
dati <- read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\time_MAV\\MAASLIN3\\tp_plus_sex_age_at_t0_bmi_2026-06-23\\all_results.tsv")

#sigificant.result del 23/06
signifcant_dati <- read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\time_MAV\\MAASLIN3\\tp_plus_sex_age_at_t0_bmi_2026-06-23\\significant_results.tsv")

str(dati)

#ho utilizzato | perchè altrimenti non ci sono corrispondenze con & 
significant_result_tp <- dati %>%
  filter(metadata == "tp") %>%
  filter(pval_joint < 0.05 | qval_joint < 0.05)

write.csv(significant_result_tp, "C:/Users/gialu/Desktop/tesi_magistrale/Git.hub/Thesis_crepaldi/Input/significant_result_tp.csv", row.names = FALSE)
