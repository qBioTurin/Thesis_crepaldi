#creo grafico association composito

#batteri
all_definitive_bacteria = read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\risultati_maaslin_bacteria\\all_results.tsv")

significant_bacteria_prevalence_pvalue <- all_definitive_bacteria %>%
  filter(model == "prevalence", metadata == "tp",
         pval_individual <= 0.05)

#write.csv(significant_bacteria_prevalence_pvalue, "C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\matrici_prev_abund\\Significant_Bacteria_prevalence_pvalue.csv", row.names = FALSE)


significant_bacteria_prevalence_qvalue <- all_definitive_bacteria %>%
  filter(model == "prevalence", metadata == "tp", qval_individual<=0.1)

#write.csv(significant_bacteria_prevalence_qvalue, "C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\matrici_prev_abund\\Significant_Bacteria_prevalence_qvalue.csv", row.names = FALSE)


significant_bacteria_abundance_pvalue <- all_definitive_bacteria %>%
  filter(model == "abundance", metadata == "tp",
         pval_individual<=0.05)

#write.csv(significant_bacteria_abundance_pvalue, "C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\matrici_prev_abund\\Significant_Bacteria_abundance_pvalue.csv", row.names = FALSE)


significant_bacteria_abundance_qvalue <- all_definitive_bacteria %>%
  filter(model == "abundance", metadata == "tp",
         qval_individual<=0.1)

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
association_alistipes <- if(sp %in% colnames(maaslin_transformed)){
  
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  
  ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
    geom_boxplot() +
    geom_jitter(width = 0.15, size = 2) +
    theme_bw() +
    labs(
      title = "Alistipes onderdonkii", 
      subtitle = "β-coefficient: 3.803 | qvalue: 0.044 | Value: T24",
      x = "Timepoint", 
      y = "Log2(relative abundance)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
  
}


#eghertella lenta
sp <- "Eggerthella_lenta"
associations_eggerthella=if(sp %in% colnames(maaslin_transformed)){
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
 
    ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
      geom_boxplot() +
      geom_jitter(width = 0.15, size = 2) +
      theme_bw() +
      labs(title = "Eggerthella lenta", 
           subtitle = "Values: T3 -> β-coefficient: -1.723 | qvalue: 0.084
           T12 -> β-coefficient: -2.158 | qvalue: 0.044",
           x = "Timepoint", y = "Log2(relative abundance)")  +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
}

#Parabacteroides_chongii
sp <- "Parabacteroides_chongii"
association_parabacteroides = if(sp %in% colnames(maaslin_transformed)){
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
 
    ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
      geom_boxplot() +
      geom_jitter(width = 0.15, size = 2) +
      theme_bw() +
      labs(title = "Parabacteroides chongii", 
           subtitle = "β-coefficient: 3.035 | qvalue: 0.044| Value: T24",
          
           x = "Timepoint", y = "Log2(relative abundance)")  +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
} 

#Collinsella_sp. i05-0019-G5
sp <- "Collinsella_sp. i05-0019-G5"
association_collinsella = if(sp %in% colnames(maaslin_transformed)){
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
 
    ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
      geom_boxplot() +
      geom_jitter(width = 0.15, size = 2) + coord_cartesian(ylim = c(-15, NA))+
      theme_bw() +
      labs(title = "Collinsella_sp. i05-0019-G5", 
           subtitle = "β-coefficient: -1.799 | qvalue: 0.052 | Value: T3",
           
           x = "Timepoint", y = "Log2(relative abundance)")  +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
} 

#Gordonibacter_pamelaeae
sp <- "Gordonibacter_pamelaeae"
assocaition_gordonibacter = if(sp %in% colnames(maaslin_transformed)){
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
 
    ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
      geom_boxplot() +
      geom_jitter(width = 0.15, size = 2) +
      theme_bw() +
      labs(title = "Gordonibacter pamelaeae", 
           subtitle = "Values: T3 -> β-coefficient: -2.397 | qvalue: 0.052 
           T12 -> β-coefficient: -2.508 | qvalue: 0.052" ,
          
           x = "Timepoint", y = "Log2(relative abundance)")  +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
}

#Segatella_copri
sp <- "Segatella_copri"
assocaition_segatella = if(sp %in% colnames(maaslin_transformed)){
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  
    ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
      geom_boxplot() +
      geom_jitter(width = 0.15, size = 2) +
      theme_bw() +
      labs(title = "Segatella copri", 
           subtitle = "β-coefficient: 1.892 | qvalue: 0.052 | Value: T3",

           x = "Timepoint", y = "Log2(relative abundance)")  +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
} 

#Coprococcus_catus
sp <- "Coprococcus_catus"
assocaition_coprpcoccus = if(sp %in% colnames(maaslin_transformed)){
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  
    ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
      geom_boxplot() +
      geom_jitter(width = 0.15, size = 2) +
      theme_bw() +
      labs(title = "Coprococcus catus", 
           subtitle = "β-coefficient: 2.363 | qvalue: 0.072 | Value: T15",
           
           x = "Timepoint", y = "Log2(relative abundance)")  +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
} 

#Bifidobacterium_angulatum
sp <- "Bifidobacterium_angulatum"
assocaition_bifidobacterium = if(sp %in% colnames(maaslin_transformed)){
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  
    ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
      geom_boxplot() +
      geom_jitter(width = 0.15, size = 2) + coord_cartesian(ylim = c(-21, NA))+
      theme_bw() +
      labs(title = "Bifidobacterium angulatum", 
           subtitle = "β-coefficient: -4.123 | qvalue: 0.099 | Value: T15",
           
           x = "Timepoint", y = "Log2(relative abundance)")  +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
} 


# Estraggo la legenda da un grafico
legend <- get_legend(
  association_alistipes +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
)

# Tolgo la legenda da tutti i grafici
association_alistipes <- association_alistipes +
  theme(legend.position = "none")

associations_eggerthella <- associations_eggerthella +
  theme(legend.position = "none")

association_parabacteroides <- association_parabacteroides +
  theme(legend.position = "none")

association_collinsella <- association_collinsella +
  theme(legend.position = "none")

assocaition_gordonibacter <- assocaition_gordonibacter +
  theme(legend.position = "none")

assocaition_segatella <- assocaition_segatella +
  theme(legend.position = "none")

assocaition_coprpcoccus <- assocaition_coprpcoccus +
  theme(legend.position = "none")

assocaition_bifidobacterium <- assocaition_bifidobacterium +
  theme(legend.position = "none")


# Creo la griglia 2 colonne × 4 righe
association_grid <- plot_grid(
  association_collinsella,assocaition_segatella,
    associations_eggerthella,assocaition_gordonibacter,
  assocaition_bifidobacterium,assocaition_coprpcoccus,
  association_alistipes,association_parabacteroides,
  ncol = 2,
  nrow = 4,
  align = "hv"
)


# Aggiungo la legenda unica a destra
association_grid_legend <- plot_grid(
  association_grid,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.12) #cosi non lo mette fuori dalla pagina non modificare più
)


# Creo il titolo generale
title <- ggdraw() +
  draw_label(
    "Bacterial species statistically significant for abundance",
    fontface = "bold",
    size = 18,
    hjust = 0.5)


# Metto titolo + grafici + legenda insieme
association_composito_batteri <- plot_grid(
  title,
  association_grid_legend,
  ncol = 1,
  rel_heights = c(0.06, 1))

ggsave(
  "association_composito_batteri_2.jpeg",
  association_composito_batteri,
  width = 14,
  height = 15,
  dpi = 300)





################archaea association
rds = readRDS("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\SUPERVISED_DEC\\Archaea_Supervised_decontam.rds")
metadata_df <- as(rds@sam_data, "data.frame")
unique(metadata_df$tp)

metadata_df$tp <- factor(
  as.character(metadata_df$tp),
  levels = c("T0", "T3", "T12", "T15", "T24", "T36"))
unique(metadata_df$tp)


sig_species <- significant_archaea_abundance_qvalue %>%
  filter(qval_individual < 0.1) %>%
  pull(feature) %>%
  unique()

maaslin_transformed <- read.delim(
  "C:/Users/gialu/Desktop/tesi_magistrale/Analisi_pts_completi_app/Risultati_Maaslin/Risultati_MAV_app/risultati_maaslin_archaea/features/data_transformed.tsv",
  row.names = 1,
  check.names = FALSE
)

#methanobrevibacter
sp <- "Methanobrevibacter_sp. TLL-48-HuF1"
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
      labs(title = "Methanobrevibacter_sp. TLL-48-HuF1",
           subtitle = "Value: T36 | β-coefficient: -0.873 | qvalue: 0.059",
           tag = "Significant Abundance",
           x = "Timepoint", y = "Log2(relative abundance)") ) +
    theme(
      plot.title = element_text(face = "bold", size = 14), plot.tag = element_markdown(hjust = 1, size = 12),
      plot.tag.position = c(0.95, 0.98))
    
} 



################eukaryota association
rds = readRDS("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\SUPERVISED_DEC\\Eukaryota_Supervised_decontam.rds")
metadata_df <- as(rds@sam_data, "data.frame")
unique(metadata_df$tp)

metadata_df$tp <- factor(
  as.character(metadata_df$tp),
  levels = c("T0", "T3", "T12", "T15", "T24", "T36"))
unique(metadata_df$tp)


sig_species <- significant_eukaryota_abundance_qvalue %>%
  filter(qval_individual < 0.1) %>%
  pull(feature) %>%
  unique()

maaslin_transformed <- read.delim(
  "C:/Users/gialu/Desktop/tesi_magistrale/Analisi_pts_completi_app/Risultati_Maaslin/Risultati_MAV_app/risultati_maaslin_eukaryota/features/data_transformed.tsv",
  row.names = 1,
  check.names = FALSE
)

#saccaromices
sp <- "Saccharomyces_cerevisiae"
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
      labs(title = "Saccharomyces cerevisiae",
           subtitle = "Values: T12 -> β-coefficient: -0.187 | qvalue: 0.010  &  T24 -> β-coefficient: -0.332 | qvalue: 0.002",
           tag = "Significant Abundance",
           x = "Timepoint", y = "Log2(relative abundance)") ) +
    theme(
      plot.title = element_text(face = "bold", size = 14), plot.tag = element_markdown(hjust = 1, size = 12),
      plot.tag.position = c(0.95, 0.98))
} 



################################
#ASSOCIATIONS SANI MALATI DA AFFIANCARE ALLE SIGNIFICATIVE
#####
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
#ASSOCIATION PLOTS
#####
rds = readRDS("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_app_healty_plus_MAV\\shared_results\\SUPERVISED_DEC\\Bacteria_Supervised_decontam.rds")
metadata_df <- as(rds@sam_data, "data.frame")
unique(metadata_df$tp)

metadata_df$group <- factor(
  as.character(metadata_df$group),
  levels = c("Healthy", "Affected"))
unique(metadata_df$group)


sig_species <- significant_bacteria_abundance_pvalue_h %>%
  filter(pval_individual < 0.05) %>%
  pull(feature) %>%
  unique()

maaslin_transformed <- read.delim(
  "C:/Users/gialu/Desktop/tesi_magistrale/Analisi_pts_completi_app/Risultati_Maaslin/Risultati_app_healty_plus_MAV/Risultati_maaslin_complete_bacteria/features/data_transformed.tsv",
  row.names = 1,
  check.names = FALSE
)

#####
sp <- "Bifidobacterium_angulatum"
if(sp %in% colnames(maaslin_transformed)){
  plot_df <- data.frame(
    group = metadata_df$group,
    abundance = maaslin_transformed[, sp]
  )
  print(
    ggplot(
      plot_df %>% filter(group == "Healthy"),
      aes(x = group, y = abundance, fill = group)) +
      geom_boxplot() +
      geom_jitter(width = 0.15, size = 2) + scale_fill_manual(
        values = c(
          "Healthy" = "#FFF68F",
          "Affected" = "#7FFFD4"
        ),
        labels = c(
          "Healthy" = "HD",
          "Affected" = "MS"
        )
      ) + coord_cartesian(ylim = c(-21, NA))+
      scale_x_discrete(labels = c("Healthy" = "HD"))+
      theme_bw() +
      labs(title = "Bifidobacterium angualtum", 
           subtitle = ".",
           tag = "Significant Abundance",
           x = "Timepoint", y = "MaAsLin3 transformed abundance: Log2(relative abundance)") ) +
    theme(
      plot.tag = element_markdown(hjust = 1, size = 12),
      plot.tag.position = c(0.95, 0.98))
} 


####### 
sp <- "Collinsella_sp. i05-0019-G5"
if(sp %in% colnames(maaslin_transformed)){
  plot_df <- data.frame(
    group = metadata_df$group,
    abundance = maaslin_transformed[, sp]
  )
  print(
    ggplot(
      plot_df %>% filter(group == "Healthy"),
      aes(x = group, y = abundance, fill = group)) +
      geom_boxplot() +
      geom_jitter(width = 0.15, size = 2) + scale_fill_manual(values = c("Healthy" = "#FFF68F"
                                                                         )) + coord_cartesian(ylim = c(-15
                                                                                                      , -6))+
      scale_x_discrete(labels = c("Healthy" = "HD"))+
      theme_bw() +
      labs(title = "Collinsella sp.i05-0019-G5", 
           subtitle = "Collinsella sp.i05-0019-G5",
           tag = "Significant Abundance",
           x = "Timepoint", y = "MaAsLin3 transformed abundance: Log2(relative abundance)") ) +
    theme(
      plot.tag = element_markdown(hjust = 1, size = 12),
      plot.tag.position = c(0.95, 0.98))
} 



#definitivo per impaginazione tesi
############################################################
# ASSOCIATION PLOTS COMPOSITI + HEALTHY DONORS
############################################################

all_definitive_bacteria <- read.delim(
  "C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\risultati_maaslin_bacteria\\all_results.tsv"
)

significant_bacteria_prevalence_pvalue <- all_definitive_bacteria %>%
  filter(
    model == "prevalence",
    metadata == "tp",
    pval_individual <= 0.05)

significant_bacteria_prevalence_qvalue <- all_definitive_bacteria %>%
  filter(
    model == "prevalence",
    metadata == "tp",
    qval_individual <= 0.1)

significant_bacteria_abundance_pvalue <- all_definitive_bacteria %>%
  filter(
    model == "abundance",
    metadata == "tp",
    pval_individual <= 0.05)

significant_bacteria_abundance_qvalue <- all_definitive_bacteria %>%
  filter(
    model == "abundance",
    metadata == "tp",
    qval_individual <= 0.1)

rds <- readRDS(
  "C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\SUPERVISED_DEC\\Bacteria_Supervised_decontam.rds")

metadata_df <- as(rds@sam_data, "data.frame")

metadata_df$tp <- factor(
  as.character(metadata_df$tp),
  levels = c("T0", "T3", "T12", "T15", "T24", "T36"))


maaslin_transformed <- read.delim(
  "C:/Users/gialu/Desktop/tesi_magistrale/Analisi_pts_completi_app/Risultati_Maaslin/Risultati_MAV_app/risultati_maaslin_bacteria/features/data_transformed.tsv",
  row.names = 1,
  check.names = FALSE)

#alistipes
sp <- "Alistipes_onderdonkii"

association_alistipes <- if (sp %in% colnames(maaslin_transformed)) {
  
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  
  ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
    geom_boxplot() +
    geom_jitter(width = 0.15, size = 2) +
    theme_bw() +
    labs(
      title = "Alistipes onderdonkii",
      subtitle = "β-coefficient: 3.803 | qvalue: 0.044 | Value: T24",
      x = "Timepoint",
      y = "Log2(relative abundance)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
  
}


# EGGERTHELLA
sp <- "Eggerthella_lenta"

associations_eggerthella <- if (sp %in% colnames(maaslin_transformed)) {
  
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  
  ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
    geom_boxplot() +
    geom_jitter(width = 0.15, size = 2) +
    theme_bw() +
    labs(
      title = "Eggerthella lenta",
      subtitle = "Values: T3 -> β-coefficient: -1.723 | qvalue: 0.084\nT12 -> β-coefficient: -2.158 | qvalue: 0.044",
      x = "Timepoint",
      y = "Log2(relative abundance)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
  
} 

#parabacteroides
sp <- "Parabacteroides_chongii"

association_parabacteroides <- if (sp %in% colnames(maaslin_transformed)) {
  
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  
  ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
    geom_boxplot() +
    geom_jitter(width = 0.15, size = 2) +
    theme_bw() +
    labs(
      title = "Parabacteroides chongii",
      subtitle = "β-coefficient: 3.035 | qvalue: 0.044 | Value: T24",
      x = "Timepoint",
      y = "Log2(relative abundance)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
  
} 

# COLLINSella
sp <- "Collinsella_sp. i05-0019-G5"

association_collinsella <- if (sp %in% colnames(maaslin_transformed)) {
  
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  
  ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
    geom_boxplot() +
    geom_jitter(width = 0.15, size = 2) +
    coord_cartesian(ylim = c(-15, -6)) +
    theme_bw() +
    labs(
      title = "Collinsella sp. i05-0019-G5",
      subtitle = "β-coefficient: -1.799 | qvalue: 0.052 | Value: T3",
      x = "Timepoint",
      y = "Log2(relative abundance)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
  
} 


#GORDONIBACTER
sp <- "Gordonibacter_pamelaeae"

association_gordonibacter <- if (sp %in% colnames(maaslin_transformed)) {
  
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  
  ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
    geom_boxplot() +
    geom_jitter(width = 0.15, size = 2) +
    theme_bw() +
    labs(
      title = "Gordonibacter pamelaeae",
      subtitle = "Values: T3 -> β-coefficient: -2.397 | qvalue: 0.052\nT12 -> β-coefficient: -2.508 | qvalue: 0.052",
      x = "Timepoint",
      y = "Log2(relative abundance)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
  
} 


#SEGATELLA
sp <- "Segatella_copri"

association_segatella <- if (sp %in% colnames(maaslin_transformed)) {
  
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  
  ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
    geom_boxplot() +
    geom_jitter(width = 0.15, size = 2) +
    theme_bw() +
    labs(
      title = "Segatella copri",
      subtitle = "β-coefficient: 1.892 | qvalue: 0.052 | Value: T3",
      x = "Timepoint",
      y = "Log2(relative abundance)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
  
} 



# COPROCOCCUS
sp <- "Coprococcus_catus"

association_coprococcus <- if (sp %in% colnames(maaslin_transformed)) {
  
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  
  ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
    geom_boxplot() +
    geom_jitter(width = 0.15, size = 2) +
    theme_bw() +
    labs(
      title = "Coprococcus catus",
      subtitle = "β-coefficient: 2.363 | qvalue: 0.072 | Value: T15",
      x = "Timepoint",
      y = "Log2(relative abundance)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
  
} 

# BIFIDOBACTERIUM
sp <- "Bifidobacterium_angulatum"

association_bifidobacterium <- if (sp %in% colnames(maaslin_transformed)) {
  
  plot_df <- data.frame(
    tp = metadata_df$tp,
    abundance = maaslin_transformed[, sp]
  )
  
  ggplot(plot_df, aes(x = tp, y = abundance, fill = tp)) +
    geom_boxplot() +
    geom_jitter(width = 0.15, size = 2) +
    coord_cartesian(ylim = c(-21, NA)) +
    theme_bw() +
    labs(
      title = "Bifidobacterium angulatum",
      subtitle = "β-coefficient: -4.123 | qvalue: 0.099 | Value: T15",
      x = "Timepoint",
      y = "Log2(relative abundance)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
  
} else {
  message("Specie non trovata: ", sp)
}



#HEALTHY DONORS

healthy_affected_bacteria <- read.delim(
  "C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_app_healty_plus_MAV\\Risultati_maaslin_complete_bacteria\\all_results.tsv")

rds_HD <- readRDS(
  "C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_app_healty_plus_MAV\\shared_results\\SUPERVISED_DEC\\Bacteria_Supervised_decontam.rds")

metadata_HD <- as(rds_HD@sam_data, "data.frame")

metadata_HD$group <- factor(
  as.character(metadata_HD$group),
  levels = c("Healthy", "Affected"))

maaslin_transformed_HD <- read.delim(
  "C:/Users/gialu/Desktop/tesi_magistrale/Analisi_pts_completi_app/Risultati_Maaslin/Risultati_app_healty_plus_MAV/Risultati_maaslin_complete_bacteria/features/data_transformed.tsv",
  row.names = 1,
  check.names = FALSE)



# hD: COLLINSella

sp <- "Collinsella_sp. i05-0019-G5"

association_HD_collinsella <- if (sp %in% colnames(maaslin_transformed_HD)) {
  
  plot_df_HD <- data.frame(
    group = metadata_HD$group,
    abundance = maaslin_transformed_HD[, sp]
  )
  
  ggplot(
    plot_df_HD %>% filter(group == "Healthy"),
    aes(x = group, y = abundance, fill = group)
  ) +
    geom_boxplot() +
    geom_jitter(width = 0.15, size = 2) +
    coord_cartesian(ylim = c(-15, -6)) +
    scale_x_discrete(labels = c("Healthy" = "HD")) +
    scale_fill_manual(
      values = c("Healthy" = "#FFF68F")
    ) +
    theme_bw() +
    labs(
      x = NULL,
      y = "Log2(relative abundance)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
  
} 



# HD: BIFIDOBACTERIUM

sp <- "Bifidobacterium_angulatum"

association_HD_bifidobacterium <- if (sp %in% colnames(maaslin_transformed_HD)) {
  
  plot_df_HD <- data.frame(
    group = metadata_HD$group,
    abundance = maaslin_transformed_HD[, sp]
  )
  
  ggplot(
    plot_df_HD %>% filter(group == "Healthy"),
    aes(x = group, y = abundance, fill = group)
  ) +
    geom_boxplot() +
    geom_jitter(width = 0.15, size = 2) +
    coord_cartesian(ylim = c(-21, -9)) +
    scale_x_discrete(labels = c("Healthy" = "HD")) +
    scale_fill_manual(
      values = c("Healthy" = "#FFF68F")
    ) +
    theme_bw() +
    labs(
      x = NULL,
      y = "Log2(relative abundance)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
  
} else {
  message("Specie non trovata: ", sp)
}



# tolgo legende singoli grafici
association_alistipes <- association_alistipes +
  theme(legend.position = "none")

associations_eggerthella <- associations_eggerthella +
  theme(legend.position = "none")

association_parabacteroides <- association_parabacteroides +
  theme(legend.position = "none")

association_collinsella <- association_collinsella +
  theme(legend.position = "none")

association_gordonibacter <- association_gordonibacter +
  theme(legend.position = "none")

association_segatella <- association_segatella +
  theme(legend.position = "none")

association_coprococcus <- association_coprococcus +
  theme(legend.position = "none")

association_bifidobacterium <- association_bifidobacterium +
  theme(legend.position = "none")

association_HD_collinsella <- association_HD_collinsella +
  theme(legend.position = "none")

association_HD_bifidobacterium <- association_HD_bifidobacterium +
  theme(legend.position = "none")



# 16.GRIGLIA ORIGINALE 
association_grid_original <- plot_grid(
  association_collinsella,
  association_segatella,
  
  associations_eggerthella,
  association_gordonibacter,
  
  association_bifidobacterium,
  association_coprococcus,
  
  association_alistipes,
  association_parabacteroides,
  
  ncol = 2,
  nrow = 4,
  align = "hv")

# Gli HD rimangono più piccoli dei grafici principali

association_HD_collinsella_small <- ggdraw() +
  draw_plot(
    association_HD_collinsella,
    x = 0.10,
    y = 0.07,
    width = 0.85,
    height = 0.75)

association_HD_bifidobacterium_small <- ggdraw() +
  draw_plot(
    association_HD_bifidobacterium,
    x = 0.10,
    y = 0.07,
    width = 0.80,
    height = 0.75)


#hd colonna

# Righe:
# 1 = HD Collinsella
# 2 = vuoto
# 3 = HD Bifidobacterium
# 4 = vuoto

HD_column <- plot_grid(
  association_HD_collinsella_small,
  ggdraw(),
  association_HD_bifidobacterium_small,
  ggdraw(),
  ncol = 1,
  nrow = 4,
  align = "v"
)

#unione

association_complete_grid <- plot_grid(
  HD_column,
  association_grid_original,
  ncol = 2,
  rel_widths = c(0.28, 2),
  align = "h")


# LEGENDA UNICA
timepoint_colors <- setNames(
  scales::hue_pal()(6),
  c("T0", "T3", "T12", "T15", "T24", "T36")
)

legend_df <- data.frame(
  x = 1:7,
  y = 1,
  category = factor(
    c(
      "T0",
      "T3",
      "T12",
      "T15",
      "T24",
      "T36",
      "HD"),
    levels = c(
      "T0",
      "T3",
      "T12",
      "T15",
      "T24",
      "T36",
      "HD")))

legend_plot <- ggplot(
  legend_df,
  aes(x = x, y = y, fill = category)) +
  geom_point(
    shape = 22,
    size = 5) +
  scale_fill_manual(
    name = NULL,
    values = c(
      timepoint_colors,
      "HD" = "#FFF68F")) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 11))

legend <- get_legend(legend_plot)

association_grid_legend <- plot_grid(
  association_complete_grid,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.12)
)



# titolo figura
title <- ggdraw() +
  draw_label(
    "Bacterial species statistically significant for abundance",
    fontface = "bold",
    size = 18,
    hjust = 0.5)



# COMPOSITO FINALE
association_composito_batteri <- plot_grid(
  title,
  association_grid_legend,
  ncol = 1,
  rel_heights = c(0.06, 1))

ggsave(
  "association_composito_batteri_HD.jpeg",
  association_composito_batteri,
  width = 18,
  height = 15,
  dpi = 300)




#########################################
####################################
#tabella generale anche con sani e pvalue
#####################

cleaned_metadata <- read.csv(
  "Input/metadata_cleaned_completo.csv",
  sep = ";"
)

cleaned_metadata_SANI <- cleaned_metadata %>%
  filter(grepl("^(HV|AMCTRL|MBCTRL|NBCTRL|RRCTRL)", ID))

names(cleaned_metadata_SANI)[names(cleaned_metadata_SANI) == "data_nascita"] <- "Data_di_nascita"
names(cleaned_metadata_SANI)[names(cleaned_metadata_SANI) == "age_at_T0"] <- "age_at_t0"

cols_keep <- c(
  "ID",
  "Data_di_nascita",
  "naive",
  "sex",
  "bmi",
  "age_at_t0")

cleaned_metadata_SANI_reduced <- cleaned_metadata_SANI %>%
  select(any_of(cols_keep))

bmi_fun <- function(x) {
  case_when(
    x < 18.5 ~ "Underweight (<18.5)",
    x >= 18.5 & x < 25 ~ "Normal weight (18.5-24.9)",
    x >= 25 ~ "Overweight (≥25)",
    TRUE ~ NA_character_)}

medas_fun <- function(x) {
  case_when(
    x >= 50 ~ "High adherence (≥50)",
    x < 50  ~ "Low adherence (<50)",
    TRUE ~ NA_character_)}

merged <- read.csv2("C:/Users/gialu/Desktop/merged.csv", stringsAsFactors = FALSE)

merged[merged$ID == "MAV32", c("bmi", "smoking_habit", "physical_activity")] <- list(33.20, "yes", "sedentary")
merged[merged$ID == "TEC25", c("bmi", "smoking_habit", "physical_activity")] <- list(21.88, "no", "sedentary")
merged[merged$ID == "MAV31", c("bmi", "smoking_habit", "physical_activity")] <- list(27.43, "no", "sedentary")
merged[merged$ID == "TEC26", c("bmi", "smoking_habit", "physical_activity")] <- list(22.46, "yes", "active")
merged[merged$ID == "TEC27", c("bmi", "smoking_habit", "physical_activity")] <- list(25.17, "no", "partially_active")
merged[merged$ID == "TEC32", c("physical_activity")] <- list(NA)
merged[merged$ID == "MAV15", c("bmi")] <- list(32)

merged <- merged %>%
  mutate(
    gadolinium_contrast = case_when(
      gadolinium_contrast %in% c("Yes", "1") ~ 1,
      gadolinium_contrast %in% c("No", "0")  ~ 0,
      TRUE ~ NA_real_
    ),
    bmi = as.numeric(bmi),
    bmi_cat = bmi_fun(bmi),
    MEDAS = as.numeric(MEDAS),
    medas_cat = medas_fun(MEDAS))

cleaned_metadata_SANI_reduced <- cleaned_metadata_SANI_reduced %>%
  mutate(bmi_cat = bmi_fun(bmi))


p_age <- wilcox.test(
  age_at_T0 ~ cohort,
  data = bind_rows(
    merged %>%
      select(age_at_T0) %>%
      mutate(cohort = "pwMS"),
    cleaned_metadata_SANI_reduced %>%
      select(age_at_t0) %>%
      rename(age_at_T0 = age_at_t0) %>%
      mutate(cohort = "HD")))$p.value

p_sex <- fisher.test(
  table(
    c(as.character(merged$sex),
      as.character(cleaned_metadata_SANI_reduced$sex)),
    c(rep("pwMS", nrow(merged)),
      rep("HD", nrow(cleaned_metadata_SANI_reduced)))))$p.value

p_bmi <- fisher.test(
  table(
    c(as.character(merged$bmi_cat),
      as.character(cleaned_metadata_SANI_reduced$bmi_cat)),
    c(rep("pwMS", nrow(merged)),
      rep("HD", nrow(cleaned_metadata_SANI_reduced)))))$p.value


n_ctrl <- nrow(cleaned_metadata_SANI_reduced)

ctrl_age <- sprintf(
  "%.0f ± %.0f",
  mean(cleaned_metadata_SANI_reduced$age_at_t0, na.rm = TRUE),
  sd(cleaned_metadata_SANI_reduced$age_at_t0, na.rm = TRUE))

ctrl_F_txt <- sprintf(
  "%d (%.0f%%)",
  sum(cleaned_metadata_SANI_reduced$sex == "F", na.rm = TRUE),
  100 * mean(cleaned_metadata_SANI_reduced$sex == "F", na.rm = TRUE))

ctrl_M_txt <- sprintf(
  "%d (%.0f%%)",
  sum(cleaned_metadata_SANI_reduced$sex == "M", na.rm = TRUE),
  100 * mean(cleaned_metadata_SANI_reduced$sex == "M", na.rm = TRUE))

ctrl_bmi <- cleaned_metadata_SANI_reduced %>%
  count(bmi_cat) %>%
  mutate(
    txt = sprintf("%d (%.1f%%)", n, 100 * n / sum(n))  )

bmi_lookup <- setNames(
  ctrl_bmi$txt,
  ctrl_bmi$bmi_cat)


tab2 <- merged %>%
  select(
    GROUP,
    age_at_T0,
    sex,
    bmi_cat,
    naive,
    gc_treatment,
    smoking_habit,
    physical_activity,
    medas_cat,
    lesion_burden,
    gadolinium_contrast,
    bone_marrow_lesions,
    subtentorial_lesions) %>%
  tbl_summary(
    by = GROUP,
    missing = "no",
    type = list(
      gadolinium_contrast ~ "dichotomous",
      gc_treatment ~ "dichotomous"),
    value = list(
      gadolinium_contrast ~ 1,
      gc_treatment ~ "positive"),
    statistic = list(
      age_at_T0 ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"),
    label = list(
      age_at_T0 ~ "Age at baseline (years)",
      sex ~ "Sex",
      bmi_cat ~ "BMI category",
      naive ~ "Treatment-naïve",
      gc_treatment ~ "Glucocorticoid at T0",
      smoking_habit ~ "Smoking habit",
      physical_activity ~ "Physical activity",
      medas_cat ~ "MEDAS adherence",
      lesion_burden ~ "Lesion burden",
      gadolinium_contrast ~ "Gadolinium contrast",
      bone_marrow_lesions ~ "Spinal cord lesions",
      subtentorial_lesions ~ "Subtentorial lesions")) %>%
  add_overall() %>%
  add_n(statistic = "{N_miss}") %>%
  add_p(
    test = list(
      age_at_T0 ~ "wilcox.test",
      all_categorical() ~ "fisher.test"),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  modify_table_body(
    ~ .x %>%
      mutate(
        HD = case_when(
          variable == "age_at_T0" &
            row_type == "label" ~ ctrl_age,
          
          variable == "sex" &
            label == "F" ~ ctrl_F_txt,
          
          variable == "sex" &
            label == "M" ~ ctrl_M_txt,
          
          variable == "bmi_cat" &
            row_type == "level" ~ bmi_lookup[label],
                    TRUE ~ "" ))) %>%
  modify_header(
    HD ~ paste0("**HD**<br>N = ", n_ctrl),
    n ~ "**Missing**" ) %>%
  bold_labels()

tab2

# Conversione da gtsummary a gt
tab2_gt <- as_gt(tab2)

# Imposta il percorso di Edge
Sys.setenv(  CHROMOTE_CHROME = "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe")

# Esporta la tabella in PNG
gt::gtsave(
  data = tab2_gt,
  filename = "tabella_anche_con_sani.png",
  zoom = 2)

####con edss
#################aggiungendo EDSS
tab2 <- merged %>%
  select(
    GROUP,
    age_at_T0,
    sex,
    bmi_cat,
    naive,
    gc_treatment,
    EDSS_T0,
    smoking_habit,
    physical_activity,
    medas_cat,
    lesion_burden,
    gadolinium_contrast,
    bone_marrow_lesions,
    subtentorial_lesions
  ) %>%
  tbl_summary(
    by = GROUP,
    missing = "no",
    type = list(
      gadolinium_contrast ~ "dichotomous",
      gc_treatment ~ "dichotomous"
    ),
    value = list(
      gadolinium_contrast ~ 1,
      gc_treatment ~ "positive"
    ),
    statistic = list(
      age_at_T0 ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    # Uniforma i decimali a 1 cifra percentuale per tutte le colonne
    digits = list(
      all_categorical() ~ c(0, 1),
      all_continuous() ~ c(0, 0)
    ),
    label = list(
      age_at_T0 ~ "Age at baseline (years)",
      sex ~ "Sex",
      bmi_cat ~ "BMI category",
      naive ~ "Treatment-naïve",
      gc_treatment ~ "Glucocorticoid at T0",
      EDSS_T0 ~ "EDSS",
      smoking_habit ~ "Smoking habit",
      physical_activity ~ "Physical activity",
      medas_cat ~ "MEDAS adherence",
      lesion_burden ~ "Lesion burden",
      gadolinium_contrast ~ "Gadolinium contrast",
      bone_marrow_lesions ~ "Spinal cord lesions",
      subtentorial_lesions ~ "Subtentorial lesions"
    )
  ) %>%
  add_overall() %>%
  add_n(statistic = "{N_miss}") %>%
  add_p(
    test = list(
      age_at_T0 ~ "wilcox.test",
      all_categorical() ~ "fisher.test"
    ),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) %>%
  modify_table_body(
    ~ .x %>%
      mutate(
        HD = case_when(
          variable == "age_at_T0" &
            row_type == "label" ~ ctrl_age,
          
          variable == "sex" &
            label == "F" ~ ctrl_F_txt,
          
          variable == "sex" &
            label == "M" ~ ctrl_M_txt,
          
          variable == "bmi_cat" &
            row_type == "level" ~ bmi_lookup[label],
          
          TRUE ~ ""
        )
      )
  ) %>%
  modify_header(
    HD ~ paste0("**HD**<br>N = ", n_ctrl),
    n ~ "**Missing**"
  ) %>%
  bold_labels()
tab2
# Conversione da gtsummary a gt
tab2_gt <- as_gt(tab2)

# Imposta il percorso di Edge
Sys.setenv( CHROMOTE_CHROME = "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe")

# Esporta la tabella in PNG
gt::gtsave(
  data = tab2_gt,
  filename = "tabella_anche_con_sani_percentuali_aggiornate.png",
  zoom = 2)



#######################
######unione dei due summary plot sani vs malati per impaginare con un unica legenda
#batteri abundance pval

healty_affected_bacteria = read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_app_healty_plus_MAV\\Risultati_maaslin_complete_bacteria\\all_results.tsv")

significant_bacteria_prevalence_pvalue_h <- healty_affected_bacteria %>%
  filter(model == "prevalence", metadata == "group",
         pval_individual < 0.05)

significant_bacteria_abundance_pvalue_h <- healty_affected_bacteria %>%
  filter(model == "abundance", metadata == "group",
         pval_individual < 0.05)



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
    title = "Bacterial species statistically significant for abundance (p-value)") +
  theme_classic()

#batteri prevalence pval
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
    title = "Bacterial species statistically significant for prevalence (p-value)") +
  theme_classic()

##########################################
alty_affected_bacteria = read.delim("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_app_healty_plus_MAV\\Risultati_maaslin_complete_bacteria\\all_results.tsv")

significant_bacteria_prevalence_pvalue_h <- healty_affected_bacteria %>%
  filter(model == "prevalence", metadata == "group",
         pval_individual < 0.05)

significant_bacteria_abundance_pvalue_h <- healty_affected_bacteria %>%
  filter(model == "abundance", metadata == "group",
         pval_individual < 0.05)

significant_bacteria_abundance_pvalue_h$neg_log10_p <- 
  -log10(significant_bacteria_abundance_pvalue_h$pval_individual)

max_abs_coef <- max(
  abs(significant_bacteria_abundance_pvalue_h$coef), 
  na.rm = TRUE
)

#abundance
p_abundance <- ggplot(
  significant_bacteria_abundance_pvalue_h,
  aes(y = reorder(feature, -pval_individual))
) +
  geom_segment(
    aes(
      x = 0,
      xend = coef,
      yend = reorder(feature, -pval_individual)
    ),
    color = "grey70",
    linewidth = 0.8
  ) +
  geom_point(
    aes(
      x = coef,
      color = neg_log10_p
    ),
    size = 4
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  scale_color_gradient(
    low = "lightblue",
    high = "darkblue",
    name = expression(-log[10](p-value))
  ) +
  scale_x_continuous(
    limits = c(-max_abs_coef, max_abs_coef)
  ) +
  labs(
    x = "Beta coefficient",
    y = "Species",
    title = "Bacterial species statistically significant for abundance (p-value)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 13)
  )

#prevalence
significant_bacteria_prevalence_pvalue_h$neg_log10_p <- 
  -log10(significant_bacteria_prevalence_pvalue_h$pval_individual)

max_abs_coef <- max(
  abs(significant_bacteria_prevalence_pvalue_h$coef), 
  na.rm = TRUE
)

p_prevalence <- ggplot(
  significant_bacteria_prevalence_pvalue_h,
  aes(y = reorder(feature, -pval_individual))
) +
  geom_segment(
    aes(
      x = 0,
      xend = coef,
      yend = reorder(feature, -pval_individual)
    ),
    color = "grey70",
    linewidth = 0.8
  ) +
  geom_point(
    aes(
      x = coef,
      color = neg_log10_p
    ),
    size = 4
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  scale_color_gradient(
    low = "lightblue",
    high = "darkblue",
    name = expression(-log[10](p-value))
  ) +
  scale_x_continuous(
    limits = c(-max_abs_coef, max_abs_coef)
  ) +
  labs(
    x = "Beta coefficient",
    y = "Species",
    title = "Bacterial species statistically significant for prevalence (p-value)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 13),
    legend.position = "none"
  )


#combino

combined_plot <- p_abundance / p_prevalence +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.justification = "top")

# A e B
combined_plot <- combined_plot +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 14))

combined_plot

ggsave(
  "bacterial_abundance_prevalence_pvalue.png",
  plot = combined_plot,
  width = 10,
  height = 14,
  units = "in",
  dpi = 300)


#################survival con scritte più grandi
merged <- read.csv2("C:/Users/gialu/Desktop/merged.csv", stringsAsFactors = FALSE)
merged[merged$ID == "MAV32", c("bmi", "smoking_habit", "physical_activity")] <- list(33.20, "yes", "sedentary")
merged[merged$ID == "TEC25", c("bmi", "smoking_habit", "physical_activity")] <- list(21.88, "no", "sedentary")
merged[merged$ID == "MAV31", c("bmi", "smoking_habit", "physical_activity")] <- list(27.43, "no", "sedentary")
merged[merged$ID == "TEC26", c("bmi", "smoking_habit", "physical_activity")] <- list(22.46, "yes", "active")
merged[merged$ID == "TEC27", c("bmi", "smoking_habit", "physical_activity")] <- list(25.17, "no", "partially_active")
merged[merged$ID == "TEC32", c( "physical_activity")] <- list(NA)
merged[merged$ID == "MAV15", c( "bmi")] <- list(32)

merged <- merged %>%
  rowwise() %>%
  mutate(
    EDSS_last = {
      x <- c(EDSS_T36, EDSS_T30, EDSS_T24,
             EDSS_T18, EDSS_T12, EDSS_T6, EDSS_T0)
      x[which(!is.na(x))[1]]
    }
  ) %>%
  ungroup()

merged <- merged %>%
  mutate(
    EDSS_worsening = case_when(
      EDSS_T0 == 0 ~ (EDSS_last - EDSS_T0) >= 1.5,
      EDSS_T0 >= 1 & EDSS_T0 <= 5 ~ (EDSS_last - EDSS_T0) >= 1,
      EDSS_T0 > 5 ~ (EDSS_last - EDSS_T0) >= 0.5,
      TRUE ~ FALSE
    )
  )

merged <- merged %>%
  mutate(
    event_6_C  = EDSS_worsening | relapse_T6  == 1 | lesioni_T6  == 1,
    event_12_C = EDSS_worsening | relapse_T12 == 1 | lesioni_T12 == 1,
    event_18_C = EDSS_worsening | relapse_T18 == 1 | lesioni_T18 == 1,
    event_24_C = EDSS_worsening | relapse_T24 == 1 | lesioni_T24 == 1,
    event_30_C = EDSS_worsening | relapse_T30 == 1 | lesioni_T30 == 1,
    event_36_C = EDSS_worsening | relapse_T36 == 1 | lesioni_T36 == 1
  )

merged <- merged %>%
  mutate(
    event_time_C = case_when(
      event_6_C  ~ 6,
      event_12_C ~ 12,
      event_18_C ~ 18,
      event_24_C ~ 24,
      event_30_C ~ 30,
      event_36_C ~ 36,
      TRUE ~ NA_real_ ),
    event_status_C = ifelse(is.na(event_time_C), 0, 1),
    time_C = ifelse(is.na(event_time_C), 36, event_time_C))

surv_obj4 <- Surv(time = merged$time_C, event = merged$event_status_C)

fit4 <- survfit(surv_obj4 ~ GROUP,data = merged)

p <- ggsurvplot(
  fit4,
  data = merged,
  
  # Curve
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  
  # Assi
  title= "Kaplan-Meier survival curve",
  xlab = "Months",
  ylab = "Cumulative NEDA-3 survival",
  break.time.by = 6,
  xlim = c(0, 36),
  
  # Etichette curve
  legend.title = "Group",
  legend.labs = c("MAV", "TEC"),
  
  # Dimensione caratteri
  font.x = c(14, "plain"),
  font.y = c(14, "plain"),
  font.tickslab = c(12, "plain"),
  font.legend = c(12, "plain"),
  font.pval = c(12, "plain"),
  
  # Risk table
  risk.table.title = "Number at risk",
  risk.table.fontsize = 4,
  risk.table.y.text = FALSE,
  
  # Aspetto
  surv.median.line = "none")

print(p)

png(
  "C:\\Users\\gialu\\Desktop\\tesi_magistrale\\varie_tesi\\pdf_immagini_app\\foto_modificate\\NEDA3_GENERALE_scritte_allargate.png",
  width = 8,
  height = 7,
  units = "in",
  res = 300)

print(
  arrange_ggsurvplots(
    list(p),
    print = FALSE,
    ncol = 1,
    nrow = 1))

dev.off()


############altra survival subtentorial scritte ingrandite

merged_data <- merged_data %>%
  mutate(
    subtento = ifelse(subtentorial_lesions == "Yes", 1, 0),
    sub_group = case_when(
      subtento == 1 ~ "MAV subtentorial +",
      subtento == 0 ~ "MAV subtentorial -"
    )
  )

merged_data$therapy_group <- factor(
  merged_data$sub_group,
  levels = c(
    "MAV subtentorial +",
    "MAV subtentorial -"
  )
)

merged_data <- merged_data %>%
  rowwise() %>%
  mutate(
    EDSS_last = {
      x <- c(
        EDSS_T36, EDSS_T30, EDSS_T24,
        EDSS_T18, EDSS_T12, EDSS_T6, EDSS_T0
      )
      x[which(!is.na(x))[1]]
    }
  ) %>%
  ungroup()

merged_data <- merged_data %>%
  mutate(
    EDSS_worsening = case_when(
      EDSS_T0 == 0 ~ (EDSS_last - EDSS_T0) >= 1.5,
      EDSS_T0 >= 1 & EDSS_T0 <= 5 ~ (EDSS_last - EDSS_T0) >= 1,
      EDSS_T0 > 5 ~ (EDSS_last - EDSS_T0) >= 0.5,
      TRUE ~ FALSE
    )
  )


merged_data <- merged_data %>%
  mutate(
    event_6_C  = EDSS_worsening |
      relapse_T6 == "Yes" |
      lesioni_T6 == "Yes",
    
    event_12_C = EDSS_worsening |
      relapse_T12 == "Yes" |
      lesioni_T12 == "Yes",
    
    event_18_C = EDSS_worsening |
      relapse_T18 == "Yes" |
      lesioni_T18 == "Yes",
    
    event_24_C = EDSS_worsening |
      relapse_T24 == "Yes" |
      lesioni_T24 == "Yes",
    
    event_30_C = EDSS_worsening |
      relapse_T30 == "Yes" |
      lesioni_T30 == "Yes",
    
    event_36_C = EDSS_worsening |
      relapse_T36 == "Yes" |
      lesioni_T36 == "Yes")


merged_data <- merged_data %>%
  mutate(
    event_time_C = case_when(
      event_6_C  ~ 6,
      event_12_C ~ 12,
      event_18_C ~ 18,
      event_24_C ~ 24,
      event_30_C ~ 30,
      event_36_C ~ 36,
      TRUE ~ NA_real_),
    
    event_status_C = ifelse(
      is.na(event_time_C),
      0,
      1),
    
    time_C = ifelse(
      is.na(event_time_C),
      36,
      event_time_C
    )
  )


surv_obj <- Surv(time = merged_data$time_C, event = merged_data$event_status_C)

fit <- survfit(surv_obj ~ sub_group,data = merged_data)

p <- ggsurvplot(
  fit,
  data = merged_data,
  
  # Curve e risk table
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  
  # Titolo
  title = "Kaplan-Meier survival curve for MAV stratified by subtentorial lesions",
  
  # Legenda
  legend.title = "Groups",
  legend.labs = c(
    "MAV subtentorial lesions -",
    "MAV subtentorial lesions +"
  ),
  
  # Assi
  xlab = "Months",
  ylab = "Cumulative NEDA-3 survival",
  break.time.by = 6,
  xlim = c(0, 36),
  
  # Dimensione caratteri
  font.title = c(14, "bold"),
  font.x = c(14, "plain"),
  font.y = c(14, "plain"),
  font.tickslab = c(12, "plain"),
  font.legend = c(12, "plain"),
  font.pval = c(12, "plain"),
  
  # Risk table
  risk.table.title = "Number at risk",
  risk.table.fontsize = 4,
  risk.table.y.text = FALSE,
  
  # Aspetto
  surv.median.line = "none")

print(p)

png(
  "C:\\Users\\gialu\\Desktop\\tesi_magistrale\\varie_tesi\\pdf_immagini_app\\foto_modificate\\NEDA3_subtentorial.png",
  width = 8,
  height = 7,
  units = "in",
  res = 300)

print(
  arrange_ggsurvplots(
    list(p),
    print = FALSE,
    ncol = 1,
    nrow = 1
  )
)

dev.off()
