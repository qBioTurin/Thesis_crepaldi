source("Libraries.R")

######risultati con pazienti totali#####
#batteri
all_definitive_bacteria = read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\risultati_maaslin_bacteria\\all_results.tsv")

significant_bacteria_prevalence_pvalue <- all_definitive_bacteria %>%
  filter(model == "prevalence", metadata == "tp",
    pval_individual < 0.05)

significant_bacteria_prevalence_qvalue <- all_definitive_bacteria %>%
  filter(model == "prevalence", metadata == "tp", qval_individual < 0.1)

significant_bacteria_abundance_pvalue <- all_definitive_bacteria %>%
  filter(model == "abundance", metadata == "tp",
    pval_individual < 0.05)

significant_bacteria_abundance_qvalue <- all_definitive_bacteria %>%
  filter(model == "abundance", metadata == "tp",
         qval_individual < 0.1)

#eukarioti
all_definitive_eukaryota = read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\risultati_maaslin_eukaryota\\all_results.tsv")

significant_eukaryota_prevalence_pvalue <- all_definitive_eukaryota %>%
  filter(model == "prevalence", metadata == "tp",
         pval_individual < 0.05)

significant_eukaryota_prevalence_qvalue <- all_definitive_eukaryota %>%
  filter(model == "prevalence", metadata == "tp",
         qval_individual < 0.1)

significant_eukaryota_abundance_pvalue <- all_definitive_eukaryota %>%
  filter(model == "abundance", metadata == "tp",
         pval_individual < 0.05)

significant_eukaryota_abundance_qvalue <- all_definitive_eukaryota %>%
  filter(model == "abundance", metadata == "tp",
         qval_individual < 0.1)

#Archaea
all_definitive_archaea = read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\risultati_maaslin_archea\\all_results.tsv")

significant_archaea_prevalence_pvalue <- all_definitive_archaea %>%
  filter(model == "prevalence", metadata == "tp",
         pval_individual < 0.05)

significant_archaea_prevalence_qvalue <- all_definitive_archaea %>%
  filter(model == "prevalence", metadata == "tp",
         qval_individual < 0.1)

significant_archaea_abundance_pvalue <- all_definitive_archaea %>%
  filter(model == "abundance", metadata == "tp",
         pval_individual < 0.05)

significant_archaea_abundance_qvalue <- all_definitive_archaea %>%
  filter(model == "abundance", metadata == "tp",
         qval_individual < 0.1)

#lettura numero di taxa / taxa table
numero_taxa_bacteria = readRDS("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\SUPERVISED_DEC\\Bacteria_Supervised_decontam.rds")
numero_taxa_archaea = readRDS("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\SUPERVISED_DEC\\Archaea_Supervised_decontam.rds")
numero_taxa_eukariota = readRDS("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\SUPERVISED_DEC\\Eukaryota_Supervised_decontam.rds")


##grafici summary plot per ogni matrice:

#batteri prevalence pvalue
significant_bacteria_prevalence_pvalue$neg_log10_p <- -log10(significant_bacteria_prevalence_pvalue$pval_individual)
max_abs_coef <- max(abs(significant_bacteria_prevalence_pvalue$coef), na.rm = TRUE)
ggplot(significant_bacteria_prevalence_pvalue,
       aes(y = reorder(feature, -pval_individual))) +
  geom_segment(aes(x = 0,
                   xend = coef,
                   yend = reorder(feature, -pval_individual)),
               color = "grey70",
               linewidth = 0.8) +
  geom_point(aes(x = coef,
                 color = neg_log10_p),
             size = 4) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "black") +
  scale_color_gradient(low = "lightblue",high = "darkblue",
    name = expression(-log[10](p-value))) +
  scale_x_continuous(limits = c(-max_abs_coef, max_abs_coef)) +
  labs(
    x = "Beta coefficient",
    y = "Species",
    title = "Significant Bacteria prevalence (pvalue)"
  ) +
  
  theme_classic()

#batteri prevalence qvalue
significant_bacteria_prevalence_qvalue$neg_log10_q <- -log10(significant_bacteria_prevalence_qvalue$qval_individual)
max_abs_coef <- max(abs(significant_bacteria_prevalence_qvalue$coef), na.rm = TRUE)

ggplot(significant_bacteria_prevalence_qvalue,
       aes(y = reorder(feature, -qval_individual))) +
  geom_segment(aes(x = 0,
                   xend = coef,
                   yend = reorder(feature, -qval_individual)),
               color = "grey70",
               linewidth = 0.8) +
  geom_point(aes(x = coef,
                 color = neg_log10_q),
             size = 4) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "black") +
  scale_color_gradient(low = "lightblue", high = "darkblue",
                       name = expression(-log[10](q-value))) +
  scale_x_continuous(limits = c(-max_abs_coef, max_abs_coef)) +
  labs(
    x = "Beta coefficient",
    y = "Species",
    title = "Significant Bacteria prevalence (qvalue)"
  ) +
  theme_classic()


#batteri abundance pvalue


#batteri abundance qvalue


#eukaryota prevalence pvalue


#eukaryota prevalence qvalue NO!!!
#eukaryota abundance qvalue NO!!!

#eukaryota abundance pvalue


#archaea prevalence pvalue e qvalue NO!!!
#archaea abundance qvalue NO!!!
#archae abundance pvalue


###############---------------------------------------------###########
###analisi matrici healthy vs affected

#batteri
healty_affected_bacteria = read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_app_healty_plus_MAV\\Risultati_maaslin_complete_bacteria\\all_results.tsv")

significant_bacteria_prevalence_pvalue <- healty_affected_bacteria %>%
  filter(model == "prevalence", metadata == "group",
         pval_individual < 0.05)

significant_bacteria_prevalence_qvalue <- healty_affected_bacteria %>%
  filter(model == "prevalence", metadata == "group", qval_individual < 0.1)

significant_bacteria_abundance_pvalue <- healty_affected_bacteria %>%
  filter(model == "abundance", metadata == "group",
         pval_individual < 0.05)

significant_bacteria_abundance_qvalue <- healty_affected_bacteria %>%
  filter(model == "abundance", metadata == "group",
         qval_individual < 0.1)

#eukarioti
healty_affected_eukaryota = read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_app_healty_plus_MAV\\Risultati_maaslin_complete_eukaryota\\all_results.tsv")

significant_eukaryota_prevalence_pvalue <- healty_affected_eukaryota %>%
  filter(model == "prevalence", metadata == "group",
         pval_individual < 0.05)

significant_eukaryota_prevalence_qvalue <- healty_affected_eukaryota %>%
  filter(model == "prevalence", metadata == "group",
         qval_individual < 0.1)

significant_eukaryota_abundance_pvalue <- healty_affected_eukaryota %>%
  filter(model == "abundance", metadata == "group",
         pval_individual < 0.05)

significant_eukaryota_abundance_qvalue <- healty_affected_eukaryota %>%
  filter(model == "abundance", metadata == "group",
         qval_individual < 0.1)

#Archaea
healty_affected_archaea = read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_app_healty_plus_MAV\\Risultati_maaslin_complete_archaea\\all_results.tsv")

significant_archaea_prevalence_pvalue <- healty_affected_archaea %>%
  filter(model == "prevalence", metadata == "group",
         pval_individual < 0.05)

significant_archaea_prevalence_qvalue <- healty_affected_archaea %>%
  filter(model == "prevalence", metadata == "group",
         qval_individual < 0.1)

significant_archaea_abundance_pvalue <- healty_affected_archaea %>%
  filter(model == "abundance", metadata == "group",
         pval_individual < 0.05)

significant_archaea_abundance_qvalue <- healty_affected_archaea %>%
  filter(model == "abundance", metadata == "group",
         qval_individual < 0.1)

