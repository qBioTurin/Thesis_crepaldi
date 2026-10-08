#analisi delle due mutazioni piu alte in ognuno
setwd("C:/Users/gialu/Desktop/GWAS ANALISI/gwas_analysis")

# Elenco dei file di output PLINK
files_gwas <- list(
  "Bone_Marrow"      = "gwas_bone_marrow_con_correzione.assoc.logistic",
  "Subtentorial"     = "gwas_subtentorial_con_correzione.assoc.logistic",
  "Gadolinium"       = "gwas_gadolinium_con_correzione.assoc.logistic",
  "Lesion_Burden"    = "gwas_lesion_burden_con_correzione.assoc.logistic",
  "Shannon_7Specie"  = "gwas_shannon_linear_con_covariate.assoc.linear",
  "Shannon_Globale"  = "gwas_shannon_global_linear_cov.assoc.linear")

top_hits_list <- list()

for (trait_name in names(files_gwas)) {
  fpath <- files_gwas[[trait_name]]

  dt <- fread(fpath) #uso dt per comodità, carico output plink
  dt_clean <- dt[TEST == "ADD" & !is.na(P) & is.finite(P)]
  
  # dal più significativo al meno
  dt_clean <- dt_clean[order(P)]
  
  # top2 varianti
  n_top <- min(2, nrow(dt_clean))
  top_snps <- dt_clean[1:n_top]
  
  # -log10(P) e rank
  top_snps[, log10_P := -log10(P)]
  top_snps[, Analysis := trait_name]
  top_snps[, Rank := 1:n_top]
  
  top_hits_list[[trait_name]] <- top_snps
}

# unisco i result
top_hits_summary <- rbindlist(top_hits_list, fill = TRUE)

cols_order <- c("Analysis", "Rank", "CHR", "SNP", "BP", "A1", "NMISS", "P", "log10_P")
cols_to_keep <- intersect(cols_order, colnames(top_hits_summary))
res_table <- top_hits_summary[, ..cols_to_keep]
print(res_table)

fwrite(res_table, "Top2_Hits_Summary_Prof.tsv", sep = "\t", quote = FALSE, col.names = TRUE)


#analisi con le top hit
top2_df <- fread("Top2_Hits_Summary_Prof.tsv")

keep_mav <- fread("keep_mav.txt", header = FALSE)
mav_ids <- keep_mav$V2

#filtro perchè altrimenti è troppo pesante per leggere tutte le righe
vcf_file <- "filtrato_genomewide_high_moderate.vcf"
vcf_lines <- readLines(vcf_file)

header_line <- grep("^#CHROM", vcf_lines, value = TRUE)
sample_names <- unlist(strsplit(header_line, "\t"))[10:length(unlist(strsplit(header_line, "\t")))]

data_lines <- vcf_lines[!grepl("^#", vcf_lines)]

vcf_dt <- fread(text = paste(data_lines, collapse = "\n"), header = FALSE)
setnames(vcf_dt, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_names))

#cambio nomi chromosom
vcf_dt[, CHROM_CLEAN := gsub("[^0-9]", "", CHROM)]
vcf_dt[CHROM %in% c("X", "chrX"), CHROM_CLEAN := "23"]
vcf_dt[CHROM %in% c("Y", "chrY"), CHROM_CLEAN := "24"]

mav_samples_in_vcf <- intersect(sample_names, mav_ids)

#pazienti x mutazioni
results <- list()

for (i in 1:nrow(top2_df)) {
  target_chr <- as.character(top2_df$CHR[i])
  target_bp  <- as.numeric(top2_df$BP[i])
  analysis   <- top2_df$Analysis[i]
  rank       <- top2_df$Rank[i]
  p_val      <- top2_df$P[i]
  
  match_row <- vcf_dt[CHROM_CLEAN == target_chr & POS == target_bp]
  
  if (nrow(match_row) > 0) {
    genotypes_raw <- unlist(match_row[1, ..mav_samples_in_vcf])

    count_het <- sum(grepl("0/1|1/0|0\\|1|1\\|0", genotypes_raw))
    count_hom_alt <- sum(grepl("1/1|1\\|1", genotypes_raw))
    count_wt <- sum(grepl("0/0|0\\|0", genotypes_raw))
    total_mutated_pts <- count_het + count_hom_alt
    
    results[[i]] <- data.table(
      Analysis = analysis,
      Rank = rank,
      CHR = target_chr,
      BP = target_bp,
      REF = match_row$REF[1],
      ALT = match_row$ALT[1],
      P_value = p_val,
      Tot_Pazienti_MAV = length(mav_samples_in_vcf),
      Pazienti_Con_Mutazione = total_mutated_pts,
      Dettaglio_Genotipi = paste0("Etero(0/1): ", count_het, " | OmoMut(1/1): ", count_hom_alt, " | WT(0/0): ", count_wt)
    )
  }
}

final_report <- rbindlist(results)
print(final_report[, .(Analysis, Rank, CHR, BP, P_value, Pazienti_Con_Mutazione, Dettaglio_Genotipi)])

fwrite(final_report, "Risultato_Finale_Conteggio_Mutazioni_Prof.tsv", sep = "\t", quote = FALSE)



###########capire in che geni sono

top2_df <- fread("Top2_Hits_Summary_Prof.tsv")

vcf_file <- "filtrato_genomewide_high_moderate.vcf"
vcf_lines <- readLines(vcf_file)

#tengo solo le prime 8 colonne
vcf_dt <- fread(text = paste(data_lines, collapse = "\n"), select = 1:8, header = FALSE)
setnames(vcf_dt, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))

vcf_dt[, CHROM_CLEAN := gsub("[^0-9]", "", CHROM)]
vcf_dt[CHROM %in% c("X", "chrX"), CHROM_CLEAN := "23"]
vcf_dt[CHROM %in% c("Y", "chrY"), CHROM_CLEAN := "24"]

#estraggo gene ed effetto/impatto
annotated_results <- list()

for (i in 1:nrow(top2_df)) {
  target_chr <- as.character(top2_df$CHR[i])
  target_bp  <- as.numeric(top2_df$BP[i])
  
  match_row <- vcf_dt[CHROM_CLEAN == target_chr & POS == target_bp]
  
  gene_found <- "Non trovato nel VCF"
  effect_found <- "N/A"
  
  if (nrow(match_row) > 0) {
    info_str <- match_row$INFO[1]
  
    if (grepl("ANN=", info_str)) {
      ann_part <- sub(".*ANN=([^;]+).*", "\\1", info_str)
      fields <- unlist(strsplit(ann_part, "\\|"))
      
      # su SnpEff Campo 2 = Effetto, Campo 4 = Simbolo Gene
      if (length(fields) >= 4) {
        effect_found <- fields[2]
        gene_found   <- fields[4]
      }
    } else if (grepl("CSQ=", info_str)) {
      csq_part <- sub(".*CSQ=([^;]+).*", "\\1", info_str)
      fields <- unlist(strsplit(csq_part, "\\|"))
      if (length(fields) >= 4) {
        effect_found <- fields[2]
        gene_found   <- fields[4]
      }}}
  
  annotated_results[[i]] <- cbind(top2_df[i], Gene_VCF = gene_found, Effetto_VCF = effect_found)}

final_vcf_annotated <- rbindlist(annotated_results)

print(final_vcf_annotated[, .(Analysis, Rank, CHR, BP, P, Gene_VCF, Effetto_VCF)])
table(final_vcf_annotated[, .(Analysis, Rank, CHR, BP, P, Gene_VCF, Effetto_VCF)])

fwrite(final_vcf_annotated, "Top2_Hits_Geni_Dal_VCF.tsv", sep = "\t", quote = FALSE)



############àanalisi top hit per tesi #con tutte le variabili che abbiamo deciso e con le covariate sex age

#subtentorial
plink_out <- fread("gwas_subtentorial_con_correzione.assoc.logistic")
plink_top <- plink_out[TEST == "ADD" & BP %in% c(25284184, 23095953)]
print(plink_top[, .(CHR, SNP, BP, A1, TEST, NMISS, OR, P)])

vcf_file <- "filtrato_genomewide_high_moderate.vcf"
vcf_lines <- readLines(vcf_file)
data_lines <- vcf_lines[!grepl("^#", vcf_lines)]

header_line <- grep("^#CHROM", vcf_lines, value = TRUE)
sample_names <- unlist(strsplit(header_line, "\t"))[10:length(unlist(strsplit(header_line, "\t")))]

vcf_dt <- fread(text = paste(data_lines, collapse = "\n"), select = c(1:5, 10:length(unlist(strsplit(header_line, "\t")))), header = FALSE)
setnames(vcf_dt, c("CHROM", "POS", "ID", "REF", "ALT", sample_names))

vcf_dt[, CHROM_CLEAN := gsub("[^0-9]", "", CHROM)]
vcf_dt[CHROM %in% c("X", "chrX"), CHROM_CLEAN := "23"]

# PYGB (BP: 25284184)
pygb_row <- vcf_dt[POS == 25284184]
if(nrow(pygb_row) > 0) {
  gt_pygb <- unlist(pygb_row[1, ..sample_names])
  mut_samples_pygb <- names(gt_pygb[grepl("0/1|1/0|1/1|0\\|1|1\\|0|1\\|1", gt_pygb)])
  cat("PYGB (Chr 20, BP 25284184) - Presente nei campioni:", paste(mut_samples_pygb, collapse = ", "), "\n")
}

# GNAZ (BP: 23095953)
gnaz_row <- vcf_dt[POS == 23095953]
if(nrow(gnaz_row) > 0) {
  gt_gnaz <- unlist(gnaz_row[1, ..sample_names])
  mut_samples_gnaz <- names(gt_gnaz[grepl("0/1|1/0|1/1|0\\|1|1\\|0|1\\|1", gt_gnaz)])
  cat("GNAZ (Chr 22, BP 23095953) - Presente nei campioni:", paste(mut_samples_gnaz, collapse = ", "), "\n")
}



# fenotipo lesioni
cov <- fread("pheno_subtentorial.txt") # o il file fenotipo delle lesioni

# i pazienti PYGB
cat("=== presence of subtentorial x PYGB ===\n")
print(cov[IID %in% c("sample_MAV13T0_EKDN230029736-1A_HJGKNDSX7_L3_R1_hg38_pe",
                     "sample_MAV17T0_EKDN230029738-1A_HJGKNDSX7_L3_R1_hg38_pe",
                     "sample_MAV21T0_EKDN230029741-1A_HJGKNDSX7_L3_R1_hg38_pe",
                     "sample_MAV26T0_EKDN230029742-1A_HJGKNDSX7_L3_R1_hg38_pe",
                     "sample_MAV28T0_EKDN230029744-1A_HJGKNDSX7_L3_R1_hg38_pe",
                     "sample_MAV30T0_EKDN230029746-1A_HJGKNDSX7_L3_R1_hg38_pe",
                     "sample_MAV33T0_R1_hg38_pe",
                     "sample_MAV34T0_EKDN240045327-1A_22FNH5LT4_L3_R1_hg38_pe")])

# i pazienti GNAZ
cat("\n=== presence of subtentorial x GNAZ ===\n")
print(cov[IID %in% c("sample_MAV02T0_EKDN230029731-1A_HJJTWDSX7_L3_R1_hg38_pe",
                     "sample_MAV14T0_EKDN230029737-1A_HJGKNDSX7_L3_R1_hg38_pe",
                     "sample_MAV23T0_EKDN240045353-1A_22MYTFLT3_L8_R1_hg38_pe",
                     "sample_MAV26T0_EKDN230029742-1A_HJGKNDSX7_L3_R1_hg38_pe")])

# 1. Shannon lineare 8 das
pheno_sh <- fread("pheno_shannon_mav_linear.txt")
col_shannon <- colnames(pheno_sh)[ncol(pheno_sh)]
setnames(pheno_sh, c("FID", "IID", "Shannon"))

vcf_file <- "filtrato_genomewide_high_moderate.vcf"
vcf_lines <- readLines(vcf_file)
data_lines <- vcf_lines[!grepl("^#", vcf_lines)]

header_line <- grep("^#CHROM", vcf_lines, value = TRUE)
sample_names <- unlist(strsplit(header_line, "\t"))[10:length(unlist(strsplit(header_line, "\t")))]

vcf_dt <- fread(text = paste(data_lines, collapse = "\n"), select = c(1:5, 10:length(unlist(strsplit(header_line, "\t")))), header = FALSE)
setnames(vcf_dt, c("CHROM", "POS", "ID", "REF", "ALT", sample_names))

vcf_dt[, CHROM_CLEAN := gsub("[^0-9]", "", CHROM)]
vcf_dt[CHROM %in% c("X", "chrX"), CHROM_CLEAN := "23"]

top_shannon <- list(
  list(gene="ANKRD11", chr="16", bp=89280283),
  list(gene="XKRX", chr="23", bp=100922849))

for (item in top_shannon) {
  row_match <- vcf_dt[CHROM_CLEAN == item$chr & POS == item$bp]
  
  if (nrow(row_match) > 0) {
    gt_data <- data.table(
      IID = sample_names,
      Genotype = unlist(row_match[1, ..sample_names])
    )
    
    gt_data[, Status_Mutazione := ifelse(grepl("0/1|1/0|1/1|0\\|1|1\\|0|1\\|1", Genotype), "MUTATO", "WILD_TYPE")]
    
    merged_sh <- merge(gt_data, pheno_sh, by = "IID")
    
    cat("---> GENE:", item$gene, "(Chr", item$chr, "- Pos:", item$bp, ")\n") #mi dice il gene e sotto mi dice chi ha lesione e mutazione
    
    mut_samples <- merged_sh[Status_Mutazione == "MUTATO"]
    wt_samples <- merged_sh[Status_Mutazione == "WILD_TYPE"]
    
    cat("Pazienti MUTATI (n =", nrow(mut_samples), "):\n")
    print(mut_samples[, .(IID, Shannon)])
    #media shannon mutati e non per confronto
    cat("  - Media MUTATI:   ", round(mean(mut_samples$Shannon, na.rm=TRUE), 4), "\n")
    cat("  - Media WILD_TYPE:", round(mean(wt_samples$Shannon, na.rm=TRUE), 4), "\n")
      }}
#altra prova per avere anche deviazione standard
cat(sprintf("PAZIENTI MUTATI   (n=%2d)  |  Media ± SD: %6.4f ± %6.4f  |  Mediana: %6.4f\n", 
            nrow(mut_samples), mean_mut, sd_mut, med_mut))
cat(sprintf("PAZIENTI WILD-TYPE (n=%2d)  |  Media ± SD: %6.4f ± %6.4f  |  Mediana: %6.4f\n", 
            nrow(wt_samples), mean_wt, sd_wt, med_wt))



##############aggiunta

vcf_file <- "filtrato_genomewide_high_moderate.vcf"
vcf_lines <- readLines(vcf_file)
data_lines <- vcf_lines[!grepl("^#", vcf_lines)]

header_line <- grep("^#CHROM", vcf_lines, value = TRUE)
sample_names <- unlist(strsplit(header_line, "\t"))[10:length(unlist(strsplit(header_line, "\t")))]

vcf_dt <- fread(text = paste(data_lines, collapse = "\n"), select = c(1:5, 10:length(unlist(strsplit(header_line, "\t")))), header = FALSE)
setnames(vcf_dt, c("CHROM", "POS", "ID", "REF", "ALT", sample_names))

vcf_dt[, CHROM_CLEAN := gsub("[^0-9]", "", CHROM)]
vcf_dt[CHROM %in% c("X", "chrX"), CHROM_CLEAN := "23"]

# prendo i top hit delle analisi delle altre lesioni
top_hits_list <- list(
  list(analysis="LESION BURDEN (Continuo)", file="pheno_lesion_burden.txt", type="continuous", 
       hits=list(list(gene="ANKRD11", chr="16", bp=89280283),
                 list(gene="ATP2B3",  chr="23", bp=153559817))),
  
  list(analysis="GADOLINIUM CONTRAST (Dicotomico)", file="pheno_gadolinium.txt", type="binary", 
       hits=list(list(gene="IL1A",   chr="2",  bp=112783737),
                 list(gene="NCAPH2", chr="22", bp=50519258))),
  
  list(analysis="BONE MARROW LESIONS (Dicotomico)", file="pheno_bone_marrow.txt", type="binary", 
       hits=list(list(gene="FBXO41", chr="2",  bp=73266564),
                 list(gene="ABCB9", chr="12", bp=122930024)))
)

  pheno <- fread(trait$file)
  setnames(pheno, c(colnames(pheno)[1], colnames(pheno)[2], "Pheno"))
  pheno[, IID := as.character(IID)]
  
  for (item in trait$hits) {
    row_match <- vcf_dt[CHROM_CLEAN == item$chr & POS == item$bp]
    
    if (nrow(row_match) > 0) {
      gt_data <- data.table(
        IID = as.character(sample_names),
        Genotype = unlist(row_match[1, ..sample_names])
      )
      
      gt_data[, Status_Mutazione := ifelse(grepl("0/1|1/0|1/1|0\\|1|1\\|0|1\\|1", Genotype), "MUTATO", "WILD_TYPE")]
      merged <- merge(gt_data, pheno, by = "IID")
      
      cat("\n---> GENE:", item$gene, "(Chr", item$chr, "- Pos:", item$bp, ")\n")
      mut_samples <- merged[Status_Mutazione == "MUTATO"]
      wt_samples  <- merged[Status_Mutazione == "WILD_TYPE"]
      
      cat("Pazienti MUTATI (n =", nrow(mut_samples), "):\n")
      print(mut_samples[, .(IID, Pheno)])
      
      if (trait$type == "continuous") {  # Statistica per fenotipi continui (Lesion Burden)
        mean_mut <- mean(mut_samples$Pheno, na.rm=TRUE)
        sd_mut   <- sd(mut_samples$Pheno, na.rm=TRUE)
        mean_wt  <- mean(wt_samples$Pheno, na.rm=TRUE)
        sd_wt    <- sd(wt_samples$Pheno, na.rm=TRUE)
        
        cat("\n Lesion Burden Continuous:\n")
        cat(sprintf("  - MUTATI   (n=%2d): Media ± SD = %6.4f ± %6.4f\n", nrow(mut_samples), mean_mut, sd_mut))
        cat(sprintf("  - WILD_TYPE(n=%2d): Media ± SD = %6.4f ± %6.4f\n", nrow(wt_samples), mean_wt, sd_wt))
        
      } else {
        # Tabella 1=No, 2=Sì
        merged[, Status_Clinico := ifelse(Pheno == 2, "SI (Caso)", "NO (Controllo)")]
        cat("\nTabella Frequenze Fenotipo (1=No, 2=Sì):\n")
        print(table(Status_Clinico = merged$Status_Clinico, Genotipo = merged$Status_Mutazione))
      }
    } 
  }
  


#####shannon globale fenotipi per tabella quindi su tutte le specie
pheno_sh_global <- fread("pheno_shannon_global_linear.txt")
setnames(pheno_sh_global, c("FID", "IID", "Shannon_Global"))

vcf_file <- "filtrato_genomewide_high_moderate.vcf"
vcf_lines <- readLines(vcf_file)
data_lines <- vcf_lines[!grepl("^#", vcf_lines)]

header_line <- grep("^#CHROM", vcf_lines, value = TRUE)
sample_names <- unlist(strsplit(header_line, "\t"))[10:length(unlist(strsplit(header_line, "\t")))]

vcf_dt <- fread(text = paste(data_lines, collapse = "\n"), select = c(1:5, 10:length(unlist(strsplit(header_line, "\t")))), header = FALSE)
setnames(vcf_dt, c("CHROM", "POS", "ID", "REF", "ALT", sample_names))

vcf_dt[, CHROM_CLEAN := gsub("[^0-9]", "", CHROM)]
vcf_dt[CHROM %in% c("X", "chrX"), CHROM_CLEAN := "23"]

# Top Hits per Shannon Global:
# Rank 1: WWC3   (Chr 23/X, BP 10079762)
# Rank 2: FAM24A (Chr 10,   BP 122912854)

top_shannon_global <- list(
  list(gene="WWC3",   chr="23", bp=10079762),
  list(gene="FAM24A", chr="10", bp=122912854))

mean_gen <- mean(pheno_sh_global$Shannon_Global, na.rm=TRUE)
sd_gen   <- sd(pheno_sh_global$Shannon_Global, na.rm=TRUE)
cat(sprintf("pazienti (N=%d)  |  media ± SD: %6.4f ± %6.4f\n", 
            nrow(pheno_sh_global), mean_gen, sd_gen))


for (item in top_shannon_global) {
  row_match <- vcf_dt[CHROM_CLEAN == item$chr & POS == item$bp]
  
  if (nrow(row_match) > 0) {
    gt_data <- data.table(
      IID = sample_names,
      Genotype = unlist(row_match[1, ..sample_names])
    )
    
    gt_data[, Status_Mutazione := ifelse(grepl("0/1|1/0|1/1|0\\|1|1\\|0|1\\|1", Genotype), "CARRIER", "WILD_TYPE")]
    
    merged_sh <- merge(gt_data, pheno_sh_global, by = "IID")
    
    cat("---> GENE:", item$gene, "(Chr", item$chr, "- Pos:", item$bp, ")\n")
    
    carrier_samples <- merged_sh[Status_Mutazione == "CARRIER"]
    wt_samples      <- merged_sh[Status_Mutazione == "WILD_TYPE"]
    
    cat("\nPazienti CARRIER (n =", nrow(carrier_samples), "):\n")
    print(carrier_samples[, .(IID, Shannon_Global)])
    
#calcolo sd e media
    mean_car <- mean(carrier_samples$Shannon_Global, na.rm=TRUE)
    sd_car   <- sd(carrier_samples$Shannon_Global, na.rm=TRUE)
    
    mean_wt  <- mean(wt_samples$Shannon_Global, na.rm=TRUE)
    sd_wt    <- sd(wt_samples$Shannon_Global, na.rm=TRUE)

    cat(sprintf("  - PAZIENTI CARRIER   (n=%2d): Media ± SD = %6.4f ± %6.4f\n", 
                nrow(carrier_samples), mean_car, sd_car))
    cat(sprintf("  - PAZIENTI WILD-TYPE (n=%2d): Media ± SD = %6.4f ± %6.4f\n", 
                nrow(wt_samples), mean_wt, sd_wt))
    cat(sprintf("  - COORTE GENERALE    (N=%2d): Media ± SD = %6.4f ± %6.4f\n", 
                nrow(pheno_sh_global), mean_gen, sd_gen))
   } 
}

