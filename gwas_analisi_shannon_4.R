merged_data = read.csv("C:\\Users\\gialu\\Desktop\\GWAS ANALISI\\gwas_analysis\\merged_data.csv")
setwd("C:/Users/gialu/Desktop/GWAS ANALISI/gwas_analysis")

dt_merged <- as.data.table(merged_data)

#ID_SHORT per il matching
if(!"ID_SHORT" %in% colnames(dt_merged)) {
  if("ID" %in% colnames(dt_merged)) {
    dt_merged[, ID_SHORT := ID]
  }
}

#keep_mav contiene solo i miei pazienti mavs enza i sani
# lista degli ID reali dei MAV dal VCF
mav_samples <- fread("keep_mav.txt", header = FALSE)
setnames(mav_samples, c("FID", "IID"))
mav_samples[, ID_SHORT := sub(".*_(MAV[0-9]+)T0.*", "\\1", IID)]

#mediana età
mav_merged <- dt_merged[ID_SHORT %in% mav_samples$ID_SHORT]
median_age <- median(mav_merged$age_at_t0, na.rm = TRUE)
cat("Mediana Età sui MAV:", median_age, "anni\n\n")

#convertire in binario le variabili per plink
export_plink_tsv <- function(dt, var_name, output_filename, type = "binary") {
  
  sub_dt <- dt[, .(ID_SHORT, RAW_VAL = get(var_name))]
  sub_dt[, PHENOTYPE := -9]
  
  if (type == "binary") {
    # 1 = No/Assente, 2 = Yes/Presente
    sub_dt[RAW_VAL %in% c(0, "0", "No", "FALSE", FALSE), PHENOTYPE := 1]
    sub_dt[RAW_VAL %in% c(1, "1", "Yes", "TRUE", TRUE),  PHENOTYPE := 2]
  } else if (type == "sex") {
    # 1 = Male, 2 = Female
    sub_dt[RAW_VAL %in% c("M", "Male", 1),   PHENOTYPE := 1]
    sub_dt[RAW_VAL %in% c("F", "Female", 2), PHENOTYPE := 2]
  } else if (type == "age") {
    # 1 = Giovani (<= mediana), 2 = Vecchi (> mediana)
    sub_dt[!is.na(RAW_VAL) & RAW_VAL <= median_age, PHENOTYPE := 1]
    sub_dt[!is.na(RAW_VAL) & RAW_VAL > median_age,  PHENOTYPE := 2]
  }
  
# Merge con gli ID del VCF
  res <- merge(mav_samples[, .(FID, IID, ID_SHORT)], sub_dt, by = "ID_SHORT", all.x = TRUE)
  res[is.na(PHENOTYPE), PHENOTYPE := -9]
  
  plink_out <- res[, .(FID, IID, PHENOTYPE)]
  
#salvo tsv
  fwrite(plink_out, output_filename, sep = "\t", quote = FALSE, col.names = TRUE)

#Generazione dei 6 file tsv con valori numerici 1 e 2 per ogni variabile e covariata

export_plink_tsv(dt_merged, "subtentorial_lesions", "pheno_subtentorial.tsv", type = "binary")
export_plink_tsv(dt_merged, "bone_marrow_lesions",  "pheno_bone_marrow.tsv",  type = "binary")
export_plink_tsv(dt_merged, "lesion_burden",       "pheno_lesion_burden.tsv", type = "binary")
export_plink_tsv(dt_merged, "gadolinium_contrast", "pheno_gadolinium.tsv",    type = "binary")

export_plink_tsv(dt_merged, "sex",                 "pheno_sex.tsv",            type = "sex")
export_plink_tsv(dt_merged, "age_at_t0",           "pheno_age_cat.tsv",        type = "age")

#############################################################
###############creazione covariata
#cambio nome per semplificare codice guida github
dt_merged <- as.data.table(merged_data)

mav_samples <- fread("keep_mav.txt", header = FALSE)
setnames(mav_samples, c("FID", "IID"))

# estraggo id breve per match con mio mergedd data
mav_samples[, ID_SHORT := sub(".*_(MAV[0-9]+)T0.*", "\\1", IID)]

if(!"ID_SHORT" %in% colnames(dt_merged)) {
  if("ID" %in% colnames(dt_merged)) {
    dt_merged[, ID_SHORT := ID]
  }
}
# sesso nel formato PLINK (1 = M, 2 = F)
dt_merged[, SEX_PLINK := -9]
dt_merged[sex %in% c("M", "Male", 1), SEX_PLINK := 1]
dt_merged[sex %in% c("F", "Female", 2), SEX_PLINK := 2]

# mediana età
mav_merged <- dt_merged[ID_SHORT %in% mav_samples$ID_SHORT]

median_age <- median(mav_merged$age_at_t0, na.rm = TRUE)

# 1 = Giovani, 2 = Vecchi
dt_merged[, AGE_CAT := ifelse(age_at_t0 > median_age, 2, 1)]

# merge
cov_data <- merge(mav_samples[, .(FID, IID, ID_SHORT)], 
                  dt_merged[, .(ID_SHORT, SEX_PLINK, AGE_CAT)], 
                  by = "ID_SHORT", 
                  all.x = TRUE)

cov_data[is.na(SEX_PLINK), SEX_PLINK := -9]
cov_data[is.na(AGE_CAT), AGE_CAT := -9]

# formattazione x plink
plink_covariates <- cov_data[, .(FID, IID, SEX = SEX_PLINK, AGE_CAT)]

fwrite(plink_covariates, "covariates.txt", sep = "\t", quote = FALSE, col.names = TRUE)
#controllo conteggi
cat("\nConteggio Sesso: Male (1) =", sum(plink_covariates$SEX == 1), "| Female (2) =", sum(plink_covariates$SEX == 2))
cat("\nConteggio Età: Giovani (1) =", sum(plink_covariates$AGE_CAT == 1), "| Vecchi (2) =", sum(plink_covariates$AGE_CAT == 2))




########################creo i txt per fare plink
dt_merged <- as.data.table(merged_data)

mav_samples <- fread("keep_mav.txt", header = FALSE)
setnames(mav_samples, c("FID", "IID"))
mav_samples[, ID_SHORT := sub(".*_(MAV[0-9]+)T0.*", "\\1", IID)]

#formato plink(1 = No/Assente, 2 = Yes/Presente, -9 = NA)
create_plink_pheno <- function(dt, var_name, output_filename) {
  sub_dt <- dt[, .(ID_SHORT, RAW_VAL = get(var_name))]
  sub_dt[, PHENO_PLINK := -9]
  sub_dt[RAW_VAL %in% c(1, "1", "Yes", "TRUE", TRUE), PHENO_PLINK := 2]
  sub_dt[RAW_VAL %in% c(0, "0", "No", "FALSE", FALSE), PHENO_PLINK := 1]
  res <- merge(mav_samples[, .(FID, IID, ID_SHORT)], sub_dt, by = "ID_SHORT", all.x = TRUE)# Merge
  res[is.na(PHENO_PLINK), PHENO_PLINK := -9]
  
  plink_out <- res[, .(FID, IID, PHENOTYPE = PHENO_PLINK)]
  
  fwrite(plink_out, output_filename, sep = "\t", quote = FALSE, col.names = TRUE)
}

create_plink_pheno(dt_merged, "subtentorial_lesions", "pheno_subtentorial.txt")
create_plink_pheno(dt_merged, "bone_marrow_lesions",  "pheno_bone_marrow.txt")
create_plink_pheno(dt_merged, "lesion_burden",       "pheno_lesion_burden.txt")
create_plink_pheno(dt_merged, "gadolinium_contrast", "pheno_gadolinium.txt")



#############################################
#####manhattan
#########à

# Funzione per pulire i dati e generare il Manhattan Plot
plot_logistic_gwas <- function(file_path, plot_title, output_png) {
  gwas_raw <- fread(file_path)
  gwas_clean <- gwas_raw[TEST == "ADD" & !is.na(P) & !is.na(BP) & !is.na(CHR)]
  gwas_clean[, CHR := as.character(CHR)] #pulizia cromosomi e nomi
  gwas_clean[CHR %in% c("X", "23"), CHR := "23"]
  gwas_clean[CHR %in% c("Y", "24"), CHR := "24"]
  gwas_clean[, CHR := gsub("[^0-9]", "", CHR)]
  gwas_clean[, CHR := as.numeric(CHR)]
  gwas_clean[, BP  := as.numeric(BP)]
  gwas_clean[, P   := as.numeric(P)]
  gwas_clean <- gwas_clean[!is.na(CHR) & CHR >= 1 & CHR <= 24 & is.finite(P) & P > 0]
  gwas_clean[SNP == "." | is.na(SNP), SNP := paste0("chr", CHR, ":", BP)]
  gwas_clean <- gwas_clean[order(CHR, BP)]
  
  gwas_plot_data <- as.data.frame(gwas_clean)
  
#mi da warning possibile soluzione
  chroms_presenti <- sort(unique(gwas_plot_data$CHR))
  labs_dinamiche <- as.character(chroms_presenti)
  labs_dinamiche[labs_dinamiche == "23"] <- "X"
  labs_dinamiche[labs_dinamiche == "24"] <- "Y"

png(output_png, width = 1600, height = 700, res = 130)
  
  manhattan(
    gwas_plot_data, 
    chr = "CHR", 
    bp = "BP", 
    p = "P", 
    snp = "SNP",
    main = plot_title,
    suggestiveline = -log10(1e-5),   # Linea blu Suggestive (p = 1e-5)
    genomewideline = -log10(3.7e-7), # Linea rossa Bonferroni (p = 3.7e-7)
    col = c("navy", "skyblue"),
    chrlabs = labs_dinamiche,
    cex.axis = 0.8)
  
  dev.off()
}

#############
#4 file con covariata

# 1. Gadolinium Contrast
plot_logistic_gwas(
  file_path   = "gwas_gadolinium_con_correzione.assoc.logistic",
  plot_title  = "GWAS Manhattan Plot - Gadolinium Contrast (Adjusted for Age & Sex)",
  output_png  = "Manhattan_Plot_Gadolinium_Corretto.png")

# 2. Subtentorial Lesions
plot_logistic_gwas(
  file_path   = "gwas_subtentorial_con_correzione.assoc.logistic",
  plot_title  = "GWAS Manhattan Plot - Subtentorial Lesions (Adjusted for Age & Sex)",
  output_png  = "Manhattan_Plot_Subtentorial_Corretto.png")

# 3. Bone Marrow Lesions
plot_logistic_gwas(
  file_path   = "gwas_bone_marrow_con_correzione.assoc.logistic",
  plot_title  = "GWAS Manhattan Plot - Bone Marrow Lesions (Adjusted for Age & Sex)",
  output_png  = "Manhattan_Plot_BoneMarrow_Corretto.png")

# 4. Lesion Burden
plot_logistic_gwas(
  file_path   = "gwas_lesion_burden_con_correzione.assoc.logistic",
  plot_title  = "GWAS Manhattan Plot - Lesion Burden (Adjusted for Age & Sex)",
  output_png  = "Manhattan_Plot_LesionBurden_Corretto.png")



############################àcreazione ultima shannon
#file shannon originle primo
shannon_df <- fread("pheno_shannon.txt")

if(ncol(shannon_df) == 3) {
  setnames(shannon_df, c("FID", "IID", "SHANNON"))
} else {
  setnames(shannon_df, c("IID", "SHANNON"))
  shannon_df[, FID := IID]
}

# solo mav per tenere solo mav da quellka che aveva anche i sani
keep_mav <- fread("keep_mav.txt", header = FALSE)
setnames(keep_mav, c("FID", "IID"))
mav_shannon <- merge(keep_mav, shannon_df, by = c("FID", "IID"), all.x = TRUE)

fwrite(mav_shannon[, .(FID, IID, SHANNON)], "pheno_shannon_mav_linear.txt", sep = "\t", quote = FALSE, col.names = TRUE)

#mediana di nuovo ma solo mav
mediana_shannon <- median(mav_shannon$SHANNON, na.rm = TRUE)
cat("Mediana dell'Indice di Shannon MAV:", mediana_shannon, "\n")

#Shannon Binario (1 = Low, 2 = High)
mav_shannon_bin <- copy(mav_shannon)
mav_shannon_bin[, PHENOTYPE := ifelse(SHANNON > mediana_shannon, 2, 1)]
mav_shannon_bin[is.na(SHANNON), PHENOTYPE := -9]

fwrite(mav_shannon_bin[, .(FID, IID, PHENOTYPE)], "pheno_shannon_mav_binary.txt", sep = "\t", quote = FALSE, col.names = TRUE)
#stampa
cat("Pazienti Low Shannon (1):", sum(mav_shannon_bin$PHENOTYPE == 1), "\n")
cat("Pazienti High Shannon (2):", sum(mav_shannon_bin$PHENOTYPE == 2), "\n")
cat("Dati mancanti (-9):", sum(mav_shannon_bin$PHENOTYPE == -9), "\n")


#######manhattan delle 2 shannon

plot_gwas_cov <- function(file_path, plot_title, output_png) {

  gwas_raw <- fread(file_path)

  gwas_clean <- gwas_raw[TEST == "ADD" & !is.na(P) & !is.na(BP) & !is.na(CHR)]
  # formattazione cromosomi
  gwas_clean[, CHR := as.character(CHR)]
  gwas_clean[CHR %in% c("X", "23"), CHR := "23"]
  gwas_clean[CHR %in% c("Y", "24"), CHR := "24"]
  gwas_clean[, CHR := gsub("[^0-9]", "", CHR)]
  
  gwas_clean[, CHR := as.numeric(CHR)]
  gwas_clean[, BP  := as.numeric(BP)]
  gwas_clean[, P   := as.numeric(P)]
  
  gwas_clean <- gwas_clean[!is.na(CHR) & CHR >= 1 & CHR <= 24 & is.finite(P) & P > 0]
  gwas_clean[SNP == "." | is.na(SNP), SNP := paste0("chr", CHR, ":", BP)]
  gwas_clean <- gwas_clean[order(CHR, BP)]
  
  gwas_plot_data <- as.data.frame(gwas_clean)
  # etichette dei cromosomi
  chroms_presenti <- sort(unique(gwas_plot_data$CHR))
  labs_dinamiche <- as.character(chroms_presenti)
  labs_dinamiche[labs_dinamiche == "23"] <- "X"
  labs_dinamiche[labs_dinamiche == "24"] <- "Y"
  
  png(output_png, width = 1600, height = 700, res = 130)
  
  manhattan(
    gwas_plot_data, 
    chr = "CHR", 
    bp = "BP", 
    p = "P", 
    snp = "SNP",
    main = plot_title,
    suggestiveline = -log10(1e-5),   # Linea blu Suggestive (p = 1e-5)
    genomewideline = -log10(3.7e-7), # Linea rossa Bonferroni (p = 3.7e-7)
    col = c("navy", "skyblue"),
    chrlabs = labs_dinamiche,
    cex.axis = 0.8)
  
  dev.off()
}

# Plot Shannon Continuo
plot_gwas_cov(
  file_path  = "gwas_shannon_linear_con_covariate.assoc.linear",
  plot_title = "mGWAS Manhattan Plot - Shannon Index Continuous (Soli MAV, Adj for Age & Sex)",
  output_png = "Manhattan_Plot_Shannon_Linear_Cov.png")

# Plot Shannon Dicotomico
plot_gwas_cov(
  file_path  = "gwas_shannon_binary_con_covariate.assoc.logistic",
  plot_title = "mGWAS Manhattan Plot - High vs Low Shannon Binary (Soli MAV, Adj for Age & Sex)",
  output_png = "Manhattan_Plot_Shannon_Binary_Cov.png")





#############################################################à
####################
#shannon totale e non solo specie della das -> global
# phyloseq del microbioma
rds_data <- readRDS("C:\\Users\\gialu\\Desktop\\GWAS ANALISI\\gwas_analysis\\Bacteria_Supervised_decontam.rds")

# calcolo dello Shannon  su tutti i campioni
otu_tab_t <- as.data.frame(t(as.data.frame(otu_table(rds_data))))
shannon_vec <- diversity(otu_tab_t, index = "shannon")

# data.table microbioma
df_microbiome <- data.table(
  RAW_ID = names(shannon_vec),
  SHANNON_VAL = as.numeric(shannon_vec))

#escludo i snai
df_microbiome <- df_microbiome[grepl("MAV", RAW_ID, ignore.case = TRUE)]

# id 
df_microbiome[, ID_SHORT := gsub("[_ -]", "", toupper(RAW_ID))]

# prendo da keep gli id vcf
keep_mav <- fread("keep_mav.txt", header = FALSE)
setnames(keep_mav, c("FID", "IID"))
mav_matched_list <- list()

#non lo trova prova così
for (i in 1:nrow(keep_mav)) {
  vcf_iid <- keep_mav$IID[i]
  vcf_fid <- keep_mav$FID[i]
  
  # Cerca tra i soli MAV del microbioma quale ID breve combacia con la stringa VCF
  match_idx <- which(sapply(df_microbiome$ID_SHORT, function(x) grepl(x, vcf_iid, fixed = TRUE)))
  
  if (length(match_idx) > 0) {
    val_shannon <- df_microbiome$SHANNON_VAL[match_idx[1]]
  } else {
    val_shannon <- NA_real_
  }
  
  mav_matched_list[[i]] <- data.table(FID = vcf_fid, IID = vcf_iid, SHANNON_VAL = val_shannon)
}

mav_shannon_all <- rbindlist(mav_matched_list)

#file shannon global tx
plink_linear <- mav_shannon_all[, .(FID, IID, SHANNON = SHANNON_VAL)]
fwrite(plink_linear, "pheno_shannon_global_linear.txt", sep = "\t", quote = FALSE, col.names = TRUE)

#file binario con mediana #questo poi l'ho scartato on veniv una bella analisi
mediana_val <- median(mav_shannon_all$SHANNON_VAL, na.rm = TRUE)

mav_shannon_all[, PHENOTYPE := -9L]
mav_shannon_all[SHANNON_VAL <= mediana_val, PHENOTYPE := 1L]
mav_shannon_all[SHANNON_VAL >  mediana_val, PHENOTYPE := 2L]

plink_binary <- mav_shannon_all[, .(FID, IID, PHENOTYPE)]
fwrite(plink_binary, "pheno_shannon_global_binary.txt", sep = "\t", quote = FALSE, col.names = TRUE)


#############i 2 manhattan globale
plot_gwas_cov <- function(file_path, plot_title, output_png) {
  
  if(!file.exists(file_path)) {
    warning(paste("File non trovato:", file_path))
    return(NULL)
  }
  
  gwas_raw <- fread(file_path)
  gwas_clean <- gwas_raw[TEST == "ADD" & !is.na(P) & !is.na(BP) & !is.na(CHR)]
  gwas_clean[, CHR := as.character(CHR)]
  gwas_clean[CHR %in% c("X", "23"), CHR := "23"]
  gwas_clean[CHR %in% c("Y", "24"), CHR := "24"]
  gwas_clean[, CHR := gsub("[^0-9]", "", CHR)]
  
  gwas_clean[, CHR := as.numeric(CHR)]
  gwas_clean[, BP  := as.numeric(BP)]
  gwas_clean[, P   := as.numeric(P)]
  
  gwas_clean <- gwas_clean[!is.na(CHR) & CHR >= 1 & CHR <= 24 & is.finite(P) & P > 0]
  gwas_clean[SNP == "." | is.na(SNP), SNP := paste0("chr", CHR, ":", BP)]
  gwas_clean <- gwas_clean[order(CHR, BP)]
  
  gwas_plot_data <- as.data.frame(gwas_clean)
  chroms_presenti <- sort(unique(gwas_plot_data$CHR))
  labs_dinamiche <- as.character(chroms_presenti)
  labs_dinamiche[labs_dinamiche == "23"] <- "X"
  labs_dinamiche[labs_dinamiche == "24"] <- "Y"
  
  png(output_png, width = 1600, height = 700, res = 130)
  
  manhattan(
    gwas_plot_data, 
    chr = "CHR", 
    bp = "BP", 
    p = "P", 
    snp = "SNP",
    main = plot_title,
    suggestiveline = -log10(1e-5),   # Linea blu Suggestive (p = 1e-5)
    genomewideline = -log10(3.7e-7), # Linea rossa Bonferroni (p = 3.7e-7)
    col = c("navy", "skyblue"),
    chrlabs = labs_dinamiche,
    cex.axis = 0.8)
  
  dev.off()
}

# Plot Shannon Globale Continuo
plot_gwas_cov(
  file_path  = "gwas_shannon_global_linear_cov.assoc.linear",
  plot_title = "mGWAS Manhattan Plot - Global Shannon Continuous (Soli MAV, Adj for Age & Sex)",
  output_png = "Manhattan_Plot_Shannon_Global_Linear_Covariate.png")

# Plot Shannon Globale Dicotomico
plot_gwas_cov(
  file_path  = "gwas_shannon_global_binary_cov.assoc.logistic",
  plot_title = "mGWAS Manhattan Plot - Global High vs Low Shannon Binary (Soli MAV, Adj for Age & Sex)",
  output_png = "Manhattan_Plot_Shannon_Global_Binary_Covariate.png")

