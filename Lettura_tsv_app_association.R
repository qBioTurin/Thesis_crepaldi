source("Libraries.R")

############################################
######risultati con pazienti totali#####
############################################

#batteri
all_definitive_bacteria = read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\risultati_maaslin_bacteria\\all_results.tsv")

significant_bacteria_prevalence_pvalue <- all_definitive_bacteria %>%
  filter(model == "prevalence", metadata == "tp",
         pval_individual < 0.05)

#write.csv(significant_bacteria_prevalence_pvalue, "C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\matrici_prev_abund\\Significant_Bacteria_prevalence_pvalue.csv", row.names = FALSE)


significant_bacteria_prevalence_qvalue <- all_definitive_bacteria %>%
  filter(model == "prevalence", metadata == "tp", qval_individual < 0.1)

#write.csv(significant_bacteria_prevalence_qvalue, "C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\matrici_prev_abund\\Significant_Bacteria_prevalence_qvalue.csv", row.names = FALSE)


significant_bacteria_abundance_pvalue <- all_definitive_bacteria %>%
  filter(model == "abundance", metadata == "tp",
         pval_individual < 0.05)

#write.csv(significant_bacteria_abundance_pvalue, "C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\matrici_prev_abund\\Significant_Bacteria_abundance_pvalue.csv", row.names = FALSE)


significant_bacteria_abundance_qvalue <- all_definitive_bacteria %>%
  filter(model == "abundance", metadata == "tp",
         qval_individual < 0.1)

#write.csv(significant_bacteria_abundance_qvalue, "C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\matrici_prev_abund\\Significant_bacteria_abundance_qvalue.csv", row.names = FALSE)

significant_bacteria_abundance_qvalue$feature

#####
#ASSOCIATION PLOTS
#####
rds = readRDS("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\SUPERVISED_DEC\\Bacteria_Supervised_decontam.rds")
metadata_df <- as(rds@sam_data, "data.frame")
unique(metadata_df$tp)

metadata_df$tp <- factor(
  as.character(metadata_df$tp),
  levels = c("T0", "T3", "T12", "T15", "T24", "T36"))
unique(metadata_df$tp)


sig_species <- significant_bacteria_abundance_qvalue %>%
  filter(qval_individual < 0.1) %>%
  pull(feature) %>%
  unique()

maaslin_transformed <- read.delim(
  "C:/Users/gialu/Desktop/tesi_magistrale/Analisi_pts_completi_app/Risultati_Maaslin/Risultati_MAV_app/risultati_maaslin_bacteria/features/data_transformed.tsv",
  row.names = 1,
  check.names = FALSE
)

#alistipes onderdonkii
sp <- "Alistipes_onderdonkii"
if(sp %in% colnames(maaslin_transformed)){
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  print(
    ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
      geom_boxplot() +
      geom_jitter(width = 0.15, size = 2) +
      theme_bw() +
      labs(title = "Alistipes onderdonkii", 
           subtitle = "Beta coefficient: 3.803 | qvalue: 0.044 | Value: T24",
           tag = "Significant Abundance",
           x = "Timepoint", y = "MaAsLin3 transformed abundance: Log2(relative abundance)") ) +
    theme(
      plot.tag = element_markdown(hjust = 1, size = 12),
      plot.tag.position = c(0.95, 0.98))
 } else {
  message("Specie non trovata: ", sp)
}

#eghertella lenta
sp <- "Eggerthella_lenta"
if(sp %in% colnames(maaslin_transformed)){
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  print(
    ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
      geom_boxplot() +
      geom_jitter(width = 0.15, size = 2) +
      theme_bw() +
      labs(title = "Eggerthella lenta", 
           subtitle = "Beta coefficient: -2.158 | qvalue: 0.044 | Value: T3-T12",
           tag = "Significant Abundance",
           x = "Timepoint", y = "MaAsLin3 transformed abundance: Log2(relative abundance)") ) +
    theme(
      plot.tag = element_markdown(hjust = 1, size = 12),
      plot.tag.position = c(0.95, 0.98))
} else {
  message("Specie non trovata: ", sp)
}

#Parabacteroides_chongii
sp <- "Parabacteroides_chongii"
if(sp %in% colnames(maaslin_transformed)){
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  print(
    ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
      geom_boxplot() +
      geom_jitter(width = 0.15, size = 2) +
      theme_bw() +
      labs(title = "Parabacteroides chongii", 
           subtitle = "Beta coefficient: 3.035 | qvalue: 0.044| Value: T24",
           tag = "Significant Abundance",
           x = "Timepoint", y = "MaAsLin3 transformed abundance: Log2(relative abundance)") ) +
    theme(
      plot.tag = element_markdown(hjust = 1, size = 12),
      plot.tag.position = c(0.95, 0.98))
} else {
  message("Specie non trovata: ", sp)
}

#Collinsella_sp. i05-0019-G5
sp <- "Collinsella_sp. i05-0019-G5"
if(sp %in% colnames(maaslin_transformed)){
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  print(
    ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
      geom_boxplot() +
      geom_jitter(width = 0.15, size = 2) +
      theme_bw() +
      labs(title = "Collinsella_sp. i05-0019-G5", 
           subtitle = "Beta coefficient: -1.799 | qvalue: 0.052 | Value: T3",
           tag = "Significant Abundance",
           x = "Timepoint", y = "MaAsLin3 transformed abundance: Log2(relative abundance)") ) +
    theme(
      plot.tag = element_markdown(hjust = 1, size = 12),
      plot.tag.position = c(0.95, 0.98))
} else {
  message("Specie non trovata: ", sp)
}

#Gordonibacter_pamelaeae
sp <- "Gordonibacter_pamelaeae"
if(sp %in% colnames(maaslin_transformed)){
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  print(
    ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
      geom_boxplot() +
      geom_jitter(width = 0.15, size = 2) +
      theme_bw() +
      labs(title = "Gordonibacter pamelaeae", 
           subtitle = "Beta coefficient: -2.508 | qvalue: 0.052 | Value: T3 & T12",
           tag = "Significant Abundance",
           x = "Timepoint", y = "MaAsLin3 transformed abundance: Log2(relative abundance)") ) +
    theme(
      plot.tag = element_markdown(hjust = 1, size = 12),
      plot.tag.position = c(0.95, 0.98))
} else {
  message("Specie non trovata: ", sp)
}

#Segatella_copri
sp <- "Segatella_copri"
if(sp %in% colnames(maaslin_transformed)){
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  print(
    ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
      geom_boxplot() +
      geom_jitter(width = 0.15, size = 2) +
      theme_bw() +
      labs(title = "Segatella copri", 
           subtitle = "Beta coefficient: 1.892 | qvalue: 0.052 | Value: T3",
           tag = "Significant Abundance",
           x = "Timepoint", y = "MaAsLin3 transformed abundance: Log2(relative abundance)") ) +
    theme(
      plot.tag = element_markdown(hjust = 1, size = 12),
      plot.tag.position = c(0.95, 0.98))
} else {
  message("Specie non trovata: ", sp)
}

#Coprococcus_catus
sp <- "Coprococcus_catus"
if(sp %in% colnames(maaslin_transformed)){
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  print(
    ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
      geom_boxplot() +
      geom_jitter(width = 0.15, size = 2) +
      theme_bw() +
      labs(title = "Coprococcus catus", 
           subtitle = "Beta coefficient: 2.363 | qvalue: 0.072 | Value: T15",
           tag = "Significant Abundance",
           x = "Timepoint", y = "MaAsLin3 transformed abundance: Log2(relative abundance)") ) +
    theme(
      plot.tag = element_markdown(hjust = 1, size = 12),
      plot.tag.position = c(0.95, 0.98))
} else {
  message("Specie non trovata: ", sp)
}

#Bifidobacterium_angulatum
sp <- "Bifidobacterium_angulatum"
if(sp %in% colnames(maaslin_transformed)){
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  print(
    ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
      geom_boxplot() +
      geom_jitter(width = 0.15, size = 2) +
      theme_bw() +
      labs(title = "Bifidobacterium angulatum", 
           subtitle = "Beta coefficient: -4.123 | qvalue: 0.099 | Value: T15",
           tag = "Significant Abundance",
           x = "Timepoint", y = "MaAsLin3 transformed abundance: Log2(relative abundance)") ) +
    theme(
      plot.tag = element_markdown(hjust = 1, size = 12),
      plot.tag.position = c(0.95, 0.98))
} else {
  message("Specie non trovata: ", sp)
}




#########
#eukarioti
##########
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
###########
#Archaea
##########
all_definitive_archaea = read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\risultati_maaslin_archaea\\all_results.tsv")

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


##########################################################################################################################
#SUMMARY PLOT
####################################
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
significant_bacteria_abundance_pvalue$neg_log10_p <- -log10(significant_bacteria_abundance_pvalue$pval_individual)
max_abs_coef <- max(abs(significant_bacteria_abundance_pvalue$coef), na.rm = TRUE)
ggplot(significant_bacteria_abundance_pvalue,
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
    title = "Significant Bacteria abundance (pvalue)"
  ) +
  
  theme_classic()


#batteri abundance qvalue
significant_bacteria_abundance_qvalue$neg_log10_q <- -log10(significant_bacteria_abundance_qvalue$qval_individual)
max_abs_coef <- max(abs(significant_bacteria_abundance_qvalue$coef), na.rm = TRUE)

ggplot(significant_bacteria_abundance_qvalue,
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
    title = "Significant Bacteria abundance (qvalue)"
  ) +
  theme_classic()

#eukaryota prevalence pvalue
significant_eukaryota_prevalence_pvalue$neg_log10_p <- -log10(significant_eukaryota_prevalence_pvalue$pval_individual)
max_abs_coef <- max(abs(significant_eukaryota_prevalence_pvalue$coef), na.rm = TRUE)
ggplot(significant_eukaryota_prevalence_pvalue,
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
    title = "Significant Eukaryota prevalence (pvalue)"
  ) +
  
  theme_classic()


#eukaryota prevalence qvalue NO!!!
#eukaryota abundance qvalue NO!!!

#eukaryota abundance pvalue
significant_eukaryota_abundance_pvalue$neg_log10_p <- -log10(significant_eukaryota_abundance_pvalue$pval_individual)
max_abs_coef <- max(abs(significant_eukaryota_abundance_pvalue$coef), na.rm = TRUE)
ggplot(significant_eukaryota_abundance_pvalue,
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
    title = "Significant Eukaryota abundance (pvalue)"
  ) +
  
  theme_classic()

#archaea prevalence pvalue e qvalue NO!!!
#archaea abundance qvalue NO!!!
#archaea abundance pvalue
significant_archaea_abundance_pvalue$neg_log10_p <- -log10(significant_archaea_abundance_pvalue$pval_individual)
max_abs_coef <- max(abs(significant_archaea_abundance_pvalue$coef), na.rm = TRUE)
ggplot(significant_archaea_abundance_pvalue,
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
    title = "Significant Archaea abundance (pvalue)"
  ) +
  
  theme_classic()

###############---------###############################################------------------------------------###########
###analisi matrici healthy vs affected
##########################################################################--------------------------------------------
#batteri
healty_affected_bacteria = read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_app_healty_plus_MAV\\Risultati_maaslin_complete_bacteria\\all_results.tsv")

significant_bacteria_prevalence_pvalue_h <- healty_affected_bacteria %>%
  filter(model == "prevalence", metadata == "group",
         pval_individual < 0.05)

#write.csv(significant_bacteria_prevalence_pvalue_h ,"C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\matrici_prev_abund\\analisi sani vs malati\\significant_bacteria_prevalence_pvalue_h.csv")

significant_bacteria_prevalence_qvalue_h <- healty_affected_bacteria %>%
  filter(model == "prevalence", metadata == "group", qval_individual < 0.1)

#write.csv(significant_bacteria_prevalence_qvalue_h ,"C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\matrici_prev_abund\\analisi sani vs malati\\significant_bacteria_prevalence_qvalue_h.csv")


significant_bacteria_abundance_pvalue_h <- healty_affected_bacteria %>%
  filter(model == "abundance", metadata == "group",
         pval_individual < 0.05)

#write.csv(significant_bacteria_abundance_pvalue_h ,"C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\matrici_prev_abund\\analisi sani vs malati\\significant_bacteria_abundance_pvalue_h.csv")

significant_bacteria_abundance_qvalue_h <- healty_affected_bacteria %>%
  filter(model == "abundance", metadata == "group",
         qval_individual < 0.1)

#write.csv(significant_bacteria_abundance_qvalue_h ,"C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\matrici_prev_abund\\analisi sani vs malati\\significant_bacteria_abundance_qvalue_h.csv")


#eukarioti
healty_affected_eukaryota = read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_app_healty_plus_MAV\\Risultati_maaslin_complete_eukaryota\\all_results.tsv")

significant_eukaryota_prevalence_pvalue_h <- healty_affected_eukaryota %>%
  filter(model == "prevalence", metadata == "group",
         pval_individual < 0.05)

significant_eukaryota_prevalence_qvalue_h <- healty_affected_eukaryota %>%
  filter(model == "prevalence", metadata == "group",
         qval_individual < 0.1)

significant_eukaryota_abundance_pvalue_h <- healty_affected_eukaryota %>%
  filter(model == "abundance", metadata == "group",
         pval_individual < 0.05)

significant_eukaryota_abundance_qvalue_h <- healty_affected_eukaryota %>%
  filter(model == "abundance", metadata == "group",
         qval_individual < 0.1)

#Archaea
healty_affected_archaea = read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_app_healty_plus_MAV\\Risultati_maaslin_complete_archaea\\all_results.tsv")

significant_archaea_prevalence_pvalue_h <- healty_affected_archaea %>%
  filter(model == "prevalence", metadata == "group",
         pval_individual < 0.05)

significant_archaea_prevalence_qvalue_h <- healty_affected_archaea %>%
  filter(model == "prevalence", metadata == "group",
         qval_individual < 0.1)

significant_archaea_abundance_pvalue_h <- healty_affected_archaea %>%
  filter(model == "abundance", metadata == "group",
         pval_individual < 0.05)

significant_archaea_abundance_qvalue_h <- healty_affected_archaea %>%
  filter(model == "abundance", metadata == "group",
         qval_individual < 0.1)
##grafici summary plot per ogni matrice:

#batteri prevalence pvalue
significant_bacteria_prevalence_pvalue_h$neg_log10_p <- -log10(significant_bacteria_prevalence_pvalue_h$pval_individual)
max_abs_coef <- max(abs(significant_bacteria_prevalence_pvalue_h$coef), na.rm = TRUE)
ggplot(significant_bacteria_prevalence_pvalue_h,
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

#batteri prevalence qvalue #######NON VA PERCHE' 0


#batteri abundance pvalue
significant_bacteria_abundance_pvalue_h$neg_log10_p <- -log10(significant_bacteria_abundance_pvalue_h$pval_individual)
max_abs_coef <- max(abs(significant_bacteria_abundance_pvalue_h$coef), na.rm = TRUE)
ggplot(significant_bacteria_abundance_pvalue_h,
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
    title = "Significant Bacteria abundance (pvalue)"
  ) +
  
  theme_classic()


####batteri abundance qvalue NO!!!!!!!


#eukaryota prevalence pvalue NO!!!!!

#eukaryota prevalence qvalue NO!!!


#eukaryota abundance qvalue 
significant_eukaryota_abundance_qvalue_h$neg_log10_q <- -log10(significant_eukaryota_abundance_qvalue_h$qval_individual)
max_abs_coef <- max(abs(significant_eukaryota_abundance_qvalue_h$coef), na.rm = TRUE)

ggplot(significant_eukaryota_abundance_qvalue_h,
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
    title = "Significant Bacteria abundance (qvalue)"
  ) +
  theme_classic()



#eukaryota abundance pvalue
significant_eukaryota_abundance_pvalue_h$neg_log10_p <- -log10(significant_eukaryota_abundance_pvalue_h$pval_individual)
max_abs_coef <- max(abs(significant_eukaryota_abundance_pvalue_h$coef), na.rm = TRUE)
ggplot(significant_eukaryota_abundance_pvalue_h,
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
    title = "Significant Eukaryota abundance (pvalue)"
  ) +
  
  theme_classic()

#archaea prevalence pvalue e qvalue NO!!!
#archaea abundance qvalue NO!!!
#archaea abundance pvalue NO!!!!








