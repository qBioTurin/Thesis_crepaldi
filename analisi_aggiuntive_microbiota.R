source("Libraries.R")
ps = readRDS("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\SUPERVISED_DEC\\Bacteria_Supervised_decontam.rds")

# Estrai componenti
otu <- as.data.frame(otu_table(ps))
tax <- as.data.frame(tax_table(ps))
meta <- as.data.frame(sample_data(ps))

# Se i taxa sono nelle colonne li trasponiamo
if (!taxa_are_rows(ps)) {
  otu <- t(otu)
}

otu <- as.data.frame(otu)

# Mantieni lo stesso ordine dei campioni
otu <- otu[, rownames(meta)]

# Calcola la prevalenza per ogni timepoint
prevalence <- lapply(unique(meta$tp), function(tp){
  
  samples_tp <- rownames(meta)[meta$tp == tp]
  
  prev <- rowMeans(otu[, samples_tp] > 0)
  
  data.frame(
    Taxon = rownames(otu),
    tp = tp,
    prevalence = prev
  )
  
}) %>%
  bind_rows()

# Aggiungi la tassonomia
tax$Taxon <- rownames(tax)

prevalence <- prevalence %>%
  left_join(tax, by = "Taxon")

# Tabella finale

prevalence_table <- prevalence %>%
  group_by(Species, tp) %>%
  summarise(
    prevalence = max(prevalence),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = tp,
    values_from = prevalence
  )


#vedo quale specie ce solo in un timepoint e non negli altri

#T0
prevalence_table %>%
  filter(
    T0 > 0,
    T3 == 0,
    T12 == 0,
    T15 == 0,
    T24 == 0,
    T36 == 0
  )
#sp. ZQ21, desulfuricans

#T3
prevalence_table %>%
  filter(
    T0 == 0,
    T3 > 0,
    T12 == 0,
    T15 == 0,
    T24 == 0,
    T36 == 0
  )
#nessuna

#T12
prevalence_table %>%
  filter(
    T0 == 0,
    T3 == 0,
    T12 > 0,
    T15 == 0,
    T24 == 0,
    T36 == 0
  )
#nessuna

#T15
prevalence_table %>%
  filter(
    T0 == 0,
    T3 == 0,
    T12 == 0,
    T15 > 0,
    T24 == 0,
    T36 == 0
  )
#flexneri

#T24
prevalence_table %>%
  filter(
    T0 == 0,
    T3 == 0,
    T12 == 0,
    T15 == 0,
    T24 > 0,
    T36 == 0
  )
#kerstersii, sp. RHBSTW-00821 

#36
prevalence_table %>%
  filter(
    T0 == 0,
    T3 == 0,
    T12 == 0,
    T15 == 0,
    T24 == 0,
    T36 > 0
  )
#nessuna


##sia T3 che T15
prevalence_table %>%
  filter(
    T0 == 0,
    T3 > 0,
    T12 == 0,
    T15 > 0,
    T24 == 0,
    T36 == 0
  )
#nessuna

##sia T0 che T12
prevalence_table %>%
  filter(
    T0 > 0,
    T3 == 0,
    T12 > 0,
    T15 == 0,
    T24 == 0,
    T36 == 0
  )
#nessuna





#######tentativo claude


ps <- readRDS("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\SUPERVISED_DEC\\Bacteria_Supervised_decontam.rds")

otu <- as.data.frame(otu_table(ps))
tax <- as.data.frame(tax_table(ps))
meta <- as.data.frame(sample_data(ps))

if (!taxa_are_rows(ps)) otu <- t(otu)
otu <- as.data.frame(otu)[, rownames(meta)]

# n campioni per timepoint (utile per giudicare l'affidabilità)
n_per_tp <- table(meta$tp)
print(n_per_tp)

# soglia minima di prevalenza per considerare un taxon "presente" in un TP
# (es. almeno il 15-20% dei campioni di quel TP, invece di >0) (IO METTO 5%)
prev_threshold <- 0.10

prevalence <- lapply(unique(meta$tp), function(tp) {
  samples_tp <- rownames(meta)[meta$tp == tp]
  prev <- rowMeans(otu[, samples_tp] > 0)
  data.frame(Taxon = rownames(otu), tp = tp, prevalence = prev)
}) %>% bind_rows()

tax$Taxon <- rownames(tax)
prevalence <- prevalence %>% left_join(tax, by = "Taxon")

# tabella wide: righe = taxon, colonne = tp, valori = prevalenza
prevalence_wide <- prevalence %>%
  select(Taxon, Species, tp, prevalence) %>%
  pivot_wider(names_from = tp, values_from = prevalence)

tp_cols <- unique(meta$tp)

# genera automaticamente il pattern binario (con soglia) per ogni taxon
pattern_df <- prevalence_wide %>%
  mutate(across(all_of(tp_cols), ~ ifelse(.x > prev_threshold, 1, 0))) %>%
  mutate(pattern = do.call(paste, c(select(., all_of(tp_cols)), sep = "-")))

# conta quanti taxa condividono ogni pattern, ordinato per rarità/interesse
pattern_summary <- pattern_df %>%
  count(pattern, sort = TRUE)

# taxa "esclusivi" di un solo timepoint (pattern con un solo 1)
exclusive_patterns <- pattern_summary %>%
  filter(lengths(regmatches(pattern, gregexpr("1", pattern))) == 1)

pattern_df %>%
  filter(pattern %in% exclusive_patterns$pattern) %>%
  select(Taxon, Species, pattern, all_of(tp_cols))

##########pseudolongum   solo al T24
#######reuteri       solo al T24
########variicola     solo t24
####### michiganensis  solo t36










################################################################
#analisi mav09 e mav10 mav 34
#################################################################


ps <- readRDS("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\SUPERVISED_DEC\\Bacteria_Supervised_decontam.rds")

####MAV09
#aggrego a livello di Genus
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)

sample_data(ps_genus)$original  
ps_mav09 <- subset_samples(ps_genus, original == "MAV09")  

#la media di abbondanza tra i suoi campioni, oppure
otu_mav09 <- as.data.frame(otu_table(ps_mav09))
if (!taxa_are_rows(ps_mav09)) otu_mav09 <- t(otu_mav09)

mean_abundance <- rowMeans(otu_mav09)

#aggiungo la tassonomia (nome Genus)
tax_mav09 <- as.data.frame(tax_table(ps_mav09))

top10_MAV09 <- data.frame(
  Taxon = names(mean_abundance),
  Genus = tax_mav09$Genus[match(names(mean_abundance), rownames(tax_mav09))],
  Abundance = mean_abundance
) %>%
  arrange(desc(Abundance)) %>%
  slice_head(n = 10)

####MAV10
#aggrego a livello di Genus
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)


sample_data(ps_genus)$original  
ps_mav10 <- subset_samples(ps_genus, original == "MAV10") 

# la media di abbondanza tra i suoi campioni, oppure
otu_mav10 <- as.data.frame(otu_table(ps_mav10))
if (!taxa_are_rows(ps_mav10)) otu_mav10 <- t(otu_mav10)

mean_abundance <- rowMeans(otu_mav10)

#aggiungo la tassonomia (nome Genus)
tax_mav10 <- as.data.frame(tax_table(ps_mav10))

top10_MAV10 <- data.frame(
  Taxon = names(mean_abundance),
  Genus = tax_mav10$Genus[match(names(mean_abundance), rownames(tax_mav10))],
  Abundance = mean_abundance
) %>%
  arrange(desc(Abundance)) %>%
  slice_head(n = 10)

##MAV34
#aggrego a livello di Genus
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)


sample_data(ps_genus)$original  
ps_mav34 <- subset_samples(ps_genus, original == "MAV34") 

# la media di abbondanza tra i suoi campioni, oppure
otu_mav34 <- as.data.frame(otu_table(ps_mav34))
if (!taxa_are_rows(ps_mav34)) otu_mav34 <- t(otu_mav34)

mean_abundance <- rowMeans(otu_mav34)

#aggiungo la tassonomia (nome Genus)
tax_mav34 <- as.data.frame(tax_table(ps_mav34))

top10_MAV34 <- data.frame(
  Taxon = names(mean_abundance),
  Genus = tax_mav34$Genus[match(names(mean_abundance), rownames(tax_mav34))],
  Abundance = mean_abundance
) %>%
  arrange(desc(Abundance)) %>%
  slice_head(n = 10)


###################################################################
#plot della library size
####################################################################
######BATTERI
otu <- otu_table(ps)
#library size
if (taxa_are_rows(ps)) {
  library_size <- colSums(as(otu, "matrix"))
} else {
  library_size <- rowSums(as(otu, "matrix"))
}
# Metadata
meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
meta$Sample <- rownames(meta)

library_df <- data.frame(
  Sample = names(library_size),
  Library_size = as.numeric(library_size)
)
library_df <- dplyr::left_join(
  library_df,
  meta,
  by = "Sample"
)
#grafico
ggplot(library_df,
       aes(x = reorder(Sample, Library_size),
           y = Library_size,
           fill = tp)) +
  geom_col() +
  labs(
    x = "Sample",
    y = "Total reads",
    fill = "Timepoint"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

#####Archaea
ps1 <- readRDS("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\SUPERVISED_DEC\\Archaea_Supervised_decontam.rds")

otu <- otu_table(ps1)
#library size
if (taxa_are_rows(ps1)) {
  library_size <- colSums(as(otu, "matrix"))
} else {
  library_size <- rowSums(as(otu, "matrix"))
}
# Metadata
meta <- data.frame(sample_data(ps1), stringsAsFactors = FALSE)
meta$Sample <- rownames(meta)

library_df <- data.frame(
  Sample = names(library_size),
  Library_size = as.numeric(library_size)
)
library_df <- dplyr::left_join(
  library_df,
  meta,
  by = "Sample"
)
#grafico
ggplot(library_df,
       aes(x = reorder(Sample, Library_size),
           y = Library_size,
           fill = tp)) +
  geom_col() +
  labs(
    x = "Sample",
    y = "Total reads",
    fill = "Timepoint"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )


#####Eukaryota
ps2 <- readRDS("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\SUPERVISED_DEC\\Eukaryota_Supervised_decontam.rds")

otu <- otu_table(ps2)
#library size
if (taxa_are_rows(ps2)) {
  library_size <- colSums(as(otu, "matrix"))
} else {
  library_size <- rowSums(as(otu, "matrix"))
}
# Metadata
meta <- data.frame(sample_data(ps2), stringsAsFactors = FALSE)
meta$Sample <- rownames(meta)

library_df <- data.frame(
  Sample = names(library_size),
  Library_size = as.numeric(library_size)
)
library_df <- dplyr::left_join(
  library_df,
  meta,
  by = "Sample"
)
#grafico
ggplot(library_df,
       aes(x = reorder(Sample, Library_size),
           y = Library_size,
           fill = tp)) +
  geom_col() +
  labs(
    x = "Sample",
    y = "Total reads",
    fill = "Timepoint"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )






