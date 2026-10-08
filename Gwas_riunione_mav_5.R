#prova con subtentorial
setwd("C:\\Users\\gialu\\Desktop\\GWAS ANALISI\\gwas_analysis")

# Verifica dei nomi campioni nel VCF
con <- file("filtrato_genomewide_high_moderate.vcf", "r")
while(TRUE) {
  line <- readLines(con, n = 1)
  if(startsWith(line, "#CHROM")) {
    vcf_ids <- unlist(strsplit(line, "\t"))[10:15] # Mostra i primi 5 campioni
    break
  }
}
close(con)

#--------------------------------------------------------
old_pheno <- fread("pheno.txt")
setnames(old_pheno, c("FID", "IID", "PHENOTYPE_OLD"))
mav_samples <- old_pheno[PHENOTYPE_OLD == 2]

#chi ha le lesioni
lesions_map <- data.table(
  ID_SHORT = c("MAV02", "MAV05", "MAV07", "MAV08", "MAV09", "MAV10", "MAV11", "MAV12", 
               "MAV13", "MAV14", "MAV15", "MAV16", "MAV17", "MAV18", "MAV19", "MAV20", 
               "MAV21", "MAV22", "MAV23", "MAV24", "MAV25", "MAV26", "MAV27", "MAV28", 
               "MAV29", "MAV30", "MAV31", "MAV32", "MAV33", "MAV34"),
  
  # 2 = Presente (Yes), 1 = Assente (No), -9 = Mancante
  PHENOTYPE = c(2, 1, 1, 2, 2, 2, 1, 2, 
                1, 2, -9, -9, 1, 2, 1, 1, 
                1, 2, 1, 2, 1, 2, 2, 1, 
                1, 1, 2, 1, 1, -9))

mav_samples[, ID_SHORT := sub(".*_(MAV[0-9]+)T0.*", "\\1", IID)]

# Merge
new_pheno <- merge(mav_samples[, .(FID, IID, ID_SHORT)], lesions_map, by = "ID_SHORT", all.x = TRUE)

#se non ce corrispondenza -9
new_pheno[is.na(PHENOTYPE), PHENOTYPE := -9]

plink_pheno <- new_pheno[, .(FID, IID, PHENOTYPE)]

# salvo xt per plink
fwrite(plink_pheno, "pheno_subtentorial.txt", sep = "\t", quote = FALSE, col.names = TRUE)

#campioni mav da includere -> analsi solo mav
fwrite(plink_pheno[, .(FID, IID)], "keep_mav.txt", sep = "\t", quote = FALSE, col.names = FALSE)

print(plink_pheno)
cat("\nCampioni MAV con lesioni (2):", sum(plink_pheno$PHENOTYPE == 2))
cat("\nCampioni MAV senza lesioni (1):", sum(plink_pheno$PHENOTYPE == 1))
cat("\nCampioni MAV con dato mancante (-9):", sum(plink_pheno$PHENOTYPE == -9), "\n")




##########################################manhattan prova 
gwas_sub <- fread("gwas_subtentorial_unfiltered_results_SUBTENTORIAL.assoc")

gwas_sub <- gwas_sub[!is.na(P) & !is.na(BP) & !is.na(CHR)]

gwas_sub[, CHR := as.character(CHR)]
gwas_sub[CHR %in% c("X", "23"), CHR := "23"]
gwas_sub[CHR %in% c("Y", "24"), CHR := "24"]
gwas_sub[, CHR := gsub("[^0-9]", "", CHR)]

gwas_sub[, CHR := as.numeric(CHR)]
gwas_sub[, BP  := as.numeric(BP)]
gwas_sub[, P   := as.numeric(P)]

gwas_sub_clean <- gwas_sub[!is.na(CHR) & CHR >= 1 & CHR <= 24 & is.finite(P) & P > 0]
gwas_sub_clean[SNP == "." | is.na(SNP), SNP := paste0("chr", CHR, ":", BP)]
gwas_sub_clean <- gwas_sub_clean[order(CHR, BP)]

gwas_plot_data <- as.data.frame(gwas_sub_clean)

# manhattan
png("Manhattan_Plot_Subtentorial_Lesions.png", width = 1200, height = 600, res = 120)

manhattan(
  gwas_plot_data, 
  chr = "CHR", 
  bp = "BP", 
  p = "P", 
  snp = "SNP",
  main = "GWAS Manhattan Plot - Lesioni Subtentoriali (Soli Pazienti MAV)",
  suggestiveline = -log10(1e-5),   # Linea blu: Suggestive standard (P = 1e-5)
  genomewideline = -log10(3.7e-7), # Linea rossa: Bonferroni (P = 3.7e-7)
  col = c("navy", "skyblue"),
  chrlabs = c(1:22, "X", "Y"))

dev.off()



###############manhattan shannon binario #poi l'abbiamo escluso
gwas_bin <- fread("gwas_shannon_binary_results.assoc")


gwas_bin <- gwas_bin[!is.na(P) & !is.na(BP) & !is.na(CHR)]

gwas_bin[, CHR := as.character(CHR)]
gwas_bin[CHR %in% c("X", "23"), CHR := "23"]
gwas_bin[CHR %in% c("Y", "24"), CHR := "24"]
gwas_bin[, CHR := gsub("[^0-9]", "", CHR)]

gwas_bin[, CHR := as.numeric(CHR)]
gwas_bin[, BP  := as.numeric(BP)]
gwas_bin[, P   := as.numeric(P)]

gwas_bin_clean <- gwas_bin[!is.na(CHR) & CHR >= 1 & CHR <= 24 & is.finite(P) & P > 0]
gwas_bin_clean[SNP == "." | is.na(SNP), SNP := paste0("chr", CHR, ":", BP)]
gwas_bin_clean <- gwas_bin_clean[order(CHR, BP)]

gwas_plot_data <- as.data.frame(gwas_bin_clean)


png("Manhattan_Plot_Shannon_Binary_MAV.png", width = 1600, height = 700, res = 130)

chroms_presenti <- sort(unique(gwas_plot_data$CHR))

labs_dinamiche <- as.character(chroms_presenti)
labs_dinamiche[labs_dinamiche == "23"] <- "X"
labs_dinamiche[labs_dinamiche == "24"] <- "Y"

manhattan(
  gwas_plot_data, 
  chr = "CHR", 
  bp = "BP", 
  p = "P", 
  snp = "SNP",
  main = "GWAS Manhattan Plot - High vs Low Shannon (Soli Pazienti MAV)",
  suggestiveline = -log10(1e-5),   # Linea blu: Suggestive standard (P = 1e-5)
  genomewideline = -log10(3.7e-7), # Linea rossa: Bonferroni (P = 3.7e-7)
  col = c("navy", "skyblue"),
  chrlabs = labs_dinamiche,        # Etichette perfettamente allineate ai dati
  cex.axis = 0.8)

dev.off()





#######################prova con il sesso
old_pheno <- fread("pheno.txt")
setnames(old_pheno, c("FID", "IID", "PHENOTYPE_OLD"))
mav_samples <- old_pheno[PHENOTYPE_OLD == 2]

# (1 = M, 2 = F)
sex_map <- data.table(
  ID_SHORT = c("MAV02", "MAV05", "MAV07", "MAV08", "MAV09", "MAV10", "MAV11", "MAV12", 
               "MAV13", "MAV14", "MAV15", "MAV16", "MAV17", "MAV18", "MAV19", "MAV20", 
               "MAV21", "MAV22", "MAV23", "MAV24", "MAV25", "MAV26", "MAV27", "MAV28", 
               "MAV29", "MAV30", "MAV31", "MAV32", "MAV33", "MAV34"),
  
  PHENOTYPE = c(1, 1, 2, 1, 2, 1, 2, 2, 
                2, 2, 2, 2, 1, 2, 1, 2, 
                2, 1, 2, 2, 2, 2, 2, 2, 
                2, 1, 2, 1, 2, 2))

# mach id
mav_samples[, ID_SHORT := sub(".*_(MAV[0-9]+)T0.*", "\\1", IID)]
new_pheno <- merge(mav_samples[, .(FID, IID, ID_SHORT)], sex_map, by = "ID_SHORT", all.x = TRUE)
new_pheno[is.na(PHENOTYPE), PHENOTYPE := -9]

plink_sex <- new_pheno[, .(FID, IID, PHENOTYPE)]

#txt
fwrite(plink_sex, "pheno_sex_mav.txt", sep = "\t", quote = FALSE, col.names = TRUE)
fwrite(plink_sex[, .(FID, IID)], "keep_mav_sex.txt", sep = "\t", quote = FALSE, col.names = FALSE)


###########manhattan sesso
gwas_sex <- fread("gwas_sex_mav_results.assoc")
gwas_sex <- gwas_sex[!is.na(P) & !is.na(BP) & !is.na(CHR)]

gwas_sex[, CHR := as.character(CHR)]
gwas_sex[CHR %in% c("X", "23"), CHR := "23"]
gwas_sex[CHR %in% c("Y", "24"), CHR := "24"]
gwas_sex[, CHR := gsub("[^0-9]", "", CHR)]

gwas_sex[, CHR := as.numeric(CHR)]
gwas_sex[, BP  := as.numeric(BP)]
gwas_sex[, P   := as.numeric(P)]

gwas_sex_clean <- gwas_sex[!is.na(CHR) & CHR >= 1 & CHR <= 24 & is.finite(P) & P > 0]
gwas_sex_clean[SNP == "." | is.na(SNP), SNP := paste0("chr", CHR, ":", BP)]
gwas_sex_clean <- gwas_sex_clean[order(CHR, BP)]

gwas_plot_data <- as.data.frame(gwas_sex_clean)

png("Manhattan_Plot_Sex_MAV.png", width = 1600, height = 700, res = 130)

chroms_presenti <- sort(unique(gwas_plot_data$CHR))
labs_dinamiche <- as.character(chroms_presenti)
labs_dinamiche[labs_dinamiche == "23"] <- "X"
labs_dinamiche[labs_dinamiche == "24"] <- "Y"

manhattan(
  gwas_plot_data, 
  chr = "CHR", 
  bp = "BP", 
  p = "P", 
  snp = "SNP",
  main = "GWAS Manhattan Plot - sex (only MAV)",
  suggestiveline = -log10(1e-5),   # Linea blu: Suggestive standard (P = 1e-5)
  genomewideline = -log10(3.7e-7), # Linea rossa: Bonferroni (P = 3.7e-7)
  col = c("navy", "skyblue"),
  chrlabs = labs_dinamiche,
  cex.axis = 0.8)

dev.off()


