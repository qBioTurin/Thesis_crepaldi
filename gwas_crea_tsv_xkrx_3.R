#creazione tsv per tool sandro #prova metriche high moderate impact
setwd("C:\\Users\\gialu\\Desktop\\GWAS ANALISI\\gwas_analysis")

# file dei fenotipi originale
pheno <- fread("pheno.txt", header = FALSE)
setnames(pheno, c("FID", "IID", "PHENOTYPE"))

# Seleziona solo IID e Fenotipo
tsv_sandro <- pheno[, .(IID, Pheno_Sandro)]

# TSV 
fwrite(tsv_sandro, "samples_pheno_sandro.tsv", sep = "\t", col.names = FALSE)

print(head(tsv_sandro))

vcf_file <- "filtrato_genomewide_high_moderate.vcf"


#####shannon 8 das
###########################lettura shannon
gwas_res <- fread("gwas_shannon_results.assoc.linear")

# elimina i valori mancanti
gwas_add <- gwas_res[TEST == "ADD" & !is.na(P)]

# Ordina per P dal più significativo
gwas_add <- gwas_add[order(P)]

# tabella filtrata in CSV
fwrite(gwas_add, "GWAS_Shannon_Top_Hits.csv")

print(head(gwas_add[, .(CHR, SNP, BP, A1, BETA, STAT, P)], 20))



#################

gwas_top <- fread("GWAS_Shannon_Top_Hits.csv")

vcf_data <- fread("filtrato_genomewide_high_moderate.vcf", skip = "#CHROM", header = FALSE, sep = "\t")

# colonne dei campioni
sample_cols <- paste0("V", 10:ncol(vcf_data))
setnames(vcf_data, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_cols))

# Estrai Gene ed Effect
vcf_data[, Gene := sub(".*ANN=[^|]*\\|[^|]*\\|[^|]*\\|([^|]*).*", "\\1", INFO)]
vcf_data[, Effect := sub(".*ANN=[^|]*\\|([^|]*).*", "\\1", INFO)]

vcf_data[, POS := as.integer(POS)]
vcf_data[, CHROM := as.character(CHROM)]

gwas_top[, BP := as.integer(BP)]
gwas_top[, CHR := as.character(CHR)]

# Merge per posizione e cromosoma 
top_annotated <- merge(
  gwas_top[P < 0.005, .(CHR, BP, A1, BETA, STAT, P)],
  vcf_data[, .(CHROM, POS, Gene, Effect)],
  by.x = c("CHR", "BP"),
  by.y = c("CHROM", "POS")
)

top_annotated <- top_annotated[order(P)]

#tab annotata
fwrite(top_annotated, "GWAS_Shannon_Annotata.csv")

print(head(top_annotated[, .(CHR, BP, Gene, Effect, BETA, P)], 15))


#####---- analisi top hit ultima analisi shannon 8 das
gwas_top <- fread("GWAS_Shannon_Top_Hits.csv")

#primi 10 segnali più forti
top_hits <- head(gwas_top, 10)

vcf_lines <- readLines("filtrato_genomewide_high_moderate.vcf")

#stampa in console cosi più facile da copiare
cat("RISULTATI TOP HITS GWAS MICROBIOTA  \n")

for(i in 1:nrow(top_hits)) {
  pos_target <- as.character(top_hits$BP[i])
  p_val <- top_hits$P[i]
  beta_val <- top_hits$BETA[i]
# Cerca la riga corrispondente nel VCF tramite la posizione
  matching_line <- grep(paste0("\t", pos_target, "\t"), vcf_lines, value = TRUE)
  
  if(length(matching_line) > 0) {
    gene <- sub(".*ANN=[^|]*\\|[^|]*\\|[^|]*\\|([^|]*).*", "\\1", matching_line[1]) # Estrai il gene
    effect <- sub(".*ANN=[^|]*\\|([^|]*).*", "\\1", matching_line[1])
    
    cat("  Posizione Genomica : Chr ", top_hits$CHR[i], " - Pos: ", pos_target, "\n", sep = "")
    cat("  GENE               : ", gene, "\n", sep = "")
    cat("  Effetto variante   : ", effect, "\n", sep = "")
    cat("  Statistica         : P-value = ", p_val, " | BETA = ", beta_val, "\n", sep = "")
    cat("---------------------------------------------------------\n")
  }
}

#######crea manhatan shannon sani vs malati
# 1. Imposta la cartella di lavoro
setwd("C:/Users/gialu/Desktop/GWAS ANALISI/gwas_analysis")
gwas_big <- fread("gwas_shannon_results.assoc.linear")

# (ADD)
gwas_clean <- gwas_big[TEST == "ADD" & !is.na(P) & !is.na(BP) & !is.na(CHR)]

# conversione tipologia colonne e X/Y
gwas_clean[, CHR := as.character(CHR)]
gwas_clean[CHR %in% c("X", "23"), CHR := "23"]
gwas_clean[CHR %in% c("Y", "24"), CHR := "24"]
gwas_clean[, CHR := gsub("[^0-9]", "", CHR)]

# formato numerico
gwas_clean[, CHR := as.numeric(CHR)]
gwas_clean[, BP  := as.numeric(BP)]
gwas_clean[, P   := as.numeric(P)]

# P-value numerici positivi
gwas_clean <- gwas_clean[!is.na(CHR) & CHR >= 1 & CHR <= 24 & is.finite(P) & P > 0]

# Se gli SNP hanno notazione . genera la posizione univoca
gwas_clean[SNP == "." | is.na(SNP), SNP := paste0("chr", CHR, ":", BP)]

# ordino per cromosoma e posizione
gwas_clean <- gwas_clean[order(CHR, BP)]
gwas_plot_data <- as.data.frame(gwas_clean)
#manhattan
png("Manhattan_Plot_Shannon_Definitivo.png", width = 1200, height = 600, res = 120)

manhattan(
  gwas_plot_data, 
  chr = "CHR", 
  bp = "BP", 
  p = "P", 
  snp = "SNP",
  main = "mGWAS Manhattan Plot - Shannon Diversity Index",
  suggestiveline = -log10(1e-4), # Linea blu per segnali suggestivi (p = 0.0001)
  genomewideline = -log10(3.7e-7), # Linea rossa Bonferroni
  col = c("navy", "skyblue"),
  chrlabs = c(1:22, "X", "Y")
)

dev.off()






#####################due mut significative
vcf_data <- fread("filtrato_genomewide_high_moderate.vcf", skip = "#CHROM", header = FALSE, sep = "\t")

sample_cols <- paste0("V", 10:ncol(vcf_data))
setnames(vcf_data, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_cols))

con <- file("filtrato_genomewide_high_moderate.vcf", "r")
while(TRUE) {
  line <- readLines(con, n = 1)
  if(startsWith(line, "#CHROM")) {
    real_sample_names <- unlist(strsplit(line, "\t"))[10:ncol(vcf_data)]
    break
  }
}
close(con)
setnames(vcf_data, sample_cols, real_sample_names)

pheno_orig <- fread("pheno.txt", header = TRUE)
setnames(pheno_orig, c("FID", "IID", "PHENOTYPE"))
pheno_orig[, Gruppo := ifelse(PHENOTYPE == 2, "MAV", "HV (Controllo)")]

#contare i fenotipi
get_genotype_counts <- function(target_pos) {
  row_var <- vcf_data[POS == target_pos]
  if(nrow(row_var) == 0) return(paste("Posizione", target_pos, "non trovata nel VCF"))
  
  # genotipo di ogni campione
  gts <- sapply(row_var[, .SD, .SDcols = real_sample_names], function(x) substr(x, 1, 3))
  
  df_gt <- data.table(IID = real_sample_names, Genotipo = gts)
  merged <- merge(df_gt, pheno_orig, by = "IID")
  
  merged[, Posizione := target_pos]
  merged[, Stato := ifelse(Genotipo %in% c("0/1", "0|1", "1/1", "1|1", "1"), "Portatore Mutazione", "WT (Riferimento)")]
  
  res <- merged[, .(
    Totale = .N,
    MAV = sum(Gruppo == "MAV"),
    HV_Controlli = sum(Gruppo == "HV (Controllo)")
  ), by = .(Posizione, Genotipo, Stato)]
  
  return(res)
}

# stampa console distribuzione
cat("=== DISTRIBUZIONE PAZIENTI - POSIZIONE 100922855 (XKRX) ===\n")
print(get_genotype_counts(100922855))

cat("\n=== DISTRIBUZIONE PAZIENTI - POSIZIONE 100922858 (XKRX) ===\n")
print(get_genotype_counts(100922858))


#####---#####---
#manhattan sani malati
setwd("C:/Users/gialu/Desktop/GWAS ANALISI/gwas_analysis")
gwas_big <- fread("gwas_mav_hv_results.assoc")

# dati mancanti
gwas_clean <- gwas_big[!is.na(P) & !is.na(BP) & !is.na(CHR)]

# pulizia
gwas_clean[, CHR := as.character(CHR)]
gwas_clean[CHR %in% c("X", "23"), CHR := "23"]
gwas_clean[CHR %in% c("Y", "24"), CHR := "24"]
gwas_clean[, CHR := gsub("[^0-9]", "", CHR)]

# Converto
gwas_clean[, CHR := as.numeric(CHR)]
gwas_clean[, BP  := as.numeric(BP)]
gwas_clean[, P   := as.numeric(P)]

#filtro
gwas_clean <- gwas_clean[!is.na(CHR) & CHR >= 1 & CHR <= 24 & is.finite(P) & P > 0]

#snp .
gwas_clean[SNP == "." | is.na(SNP), SNP := paste0("chr", CHR, ":", BP)]

gwas_clean <- gwas_clean[order(CHR, BP)]

gwas_plot_data <- as.data.frame(gwas_clean)

#manhattan
png("Manhattan_Plot_MAV_vs_HV_Definitivo.png", width = 1200, height = 600, res = 120)

manhattan(
  gwas_plot_data, 
  chr = "CHR", 
  bp = "BP", 
  p = "P", 
  snp = "SNP",
  main = "GWAS Manhattan Plot - MAV vs Controlli (High/Moderate Impact)",
  suggestiveline = -log10(1e-4),   # Linea blu: Suggestive standard (p = 0.00001)
  genomewideline = -log10(3.7e-7), # Linea rossa: Bonferroni (p = 3.7e-7)
  col = c("navy", "skyblue"),
  chrlabs = c(1:22, "X", "Y")
)

dev.off()



#######################################
#######################################
#########################################
###########################################
vcf_data <- fread("filtrato_genomewide_high_moderate.vcf", skip = "#CHROM", header = FALSE, sep = "\t")

# Recupera le intestazioni e i nomi campioni reali
con <- file("filtrato_genomewide_high_moderate.vcf", "r")
while(TRUE) {
  line <- readLines(con, n = 1)
  if(startsWith(line, "#CHROM")) {
    real_sample_names <- unlist(strsplit(line, "\t"))[10:ncol(vcf_data)]
    break
  }
}
close(con)

# riga della posizione top hit (XKRX)
row_var <- vcf_data[V2 == 100922858]

#genotipi di ciascun paziente
gts <- sapply(row_var[, 10:ncol(vcf_data), with = FALSE], function(x) substr(x, 1, 3))

df_gt <- data.table(IID = real_sample_names, Genotipo = gts)

# pheno per abbinare la diagnosi
pheno <- fread("pheno.txt", header = TRUE)
setnames(pheno, c("FID", "IID", "PHENOTYPE"))
pheno[, Gruppo := ifelse(PHENOTYPE == 2, "MAV", "HV (Controllo)")]

res <- merge(df_gt, pheno, by = "IID")
res[, Stato := ifelse(Genotipo %in% c("0/1", "0|1", "1/1", "1|1", "1"), "MUTATO", "WT")]

#stampa 5 pazienti portatori
cat("=== I 5 PAZIENTI PORTATORI (POS 100922855 - XKRX) ===\n")
print(res[Stato == "MUTATO", .(IID, Gruppo, Genotipo)])
