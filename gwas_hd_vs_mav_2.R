
setwd("C:\\Users\\gialu\\Desktop\\GWAS ANALISI\\gwas_analysis")

#sani vs malati
#ANALISI HIGH MODERATE IMPACT filtato su snpsfit
##################################################################################àà

# 1. Carica i fenotipi
pheno <- fread("pheno.txt", header = FALSE)
setnames(pheno, c("FID", "IID", "PHENOTYPE"))
pheno[, Group := ifelse(PHENOTYPE == 2, "MAV", "Controllo")]

# il VCF filtrato da SnpSift
vcf_file <- "filtrato_genomewide_high_moderate.vcf"

#Estrazione dell'intestazione dei pazienti perchè sono diversi tra il mio vcf e il mio dataset (#CHROM)
con <- file(vcf_file, "r")
header_line <- ""
while(TRUE) {
  line <- readLines(con, n = 1)
  if(length(line) == 0) break
  if(startsWith(line, "#CHROM")) {
    header_line <- line
    break
  }
}
close(con)

#nomi dei 41 pazienti sia sani che malati
sample_names <- unlist(strsplit(header_line, "\t"))[10:length(unlist(strsplit(header_line, "\t")))]

# Carica solo le righe delle varianti
vcf_data <- fread(vcf_file, skip = "#CHROM", header = FALSE, sep = "\t")

# Analisi per ogni variante trovata
  results <- list()
  
  for (i in 1:nrow(vcf_data)) {
    row <- vcf_data[i]
    info_str <- row$INFO
    ann_part <- sub(".*ANN=", "", info_str)# Estrazione Gene e Impatto da SnpEff (ANN)
    gene_info <- unlist(strsplit(ann_part, "\\|"))
    impact <- gene_info[2]
    gene <- gene_info[4]
    feature <- gene_info[11] # che cambio è avvenuto tipo un amminoacido
    gts <- sapply(row[, .SD, .SDcols = sample_names], function(x) substr(x, 1, 3))
    carriers <- sample_names[gts %in% c("0/1", "0|1", "1/0", "1|0", "1/1", "1|1")]# Identifica i portatori (genotipi 0/1, 0|1, 1/1, 1|1)
    mav_count <- sum(carriers %in% pheno[Group == "MAV", IID])# conte fenotipi
    ctrl_count <- sum(carriers %in% pheno[Group == "Controllo", IID])
    
    results[[i]] <- data.table(
      Chr = row$CHROM,
      Pos = row$POS,
      Ref = row$REF,
      Alt = row$ALT,
      Gene = gene,
      Impatto = impact,
      Modifica = feature,
      Totale_Portatori = length(carriers),
      Pazienti_MAV = mav_count,
      Pazienti_Controlli = ctrl_count)
  }
  
  summary_table <- rbindlist(results) #tabella riassuntiva dei pz e delle varianti di ciascuno
  print(summary_table)
  
# i risultati in CSV per la tesi
  fwrite(summary_table, "Varianti_High_Moderate_Sottogruppi.csv")


##############################################################################
##prova 2
#################

# Carica i fenotipi
pheno <- fread("pheno.txt", header = FALSE)
setnames(pheno, c("FID", "IID", "PHENOTYPE"))
pheno[, Group := ifelse(PHENOTYPE == 2, "MAV", "Controllo")]

mav_samples  <- pheno[Group == "MAV", IID]
ctrl_samples <- pheno[Group == "Controllo", IID]

# VCF
vcf_file <- "filtrato_genomewide_high_moderate.vcf"

# l'intestazione dei pazienti
con <- file(vcf_file, "r")
header_line <- ""
while(TRUE) {
  line <- readLines(con, n = 1)
  if(length(line) == 0) break
  if(startsWith(line, "#CHROM")) {
    header_line <- line
    break
  }
}
close(con)

sample_names <- unlist(strsplit(header_line, "\t"))[10:length(unlist(strsplit(header_line, "\t")))]

# solo i campioni presenti anche nel file dei fenotipi
mav_in_vcf  <- intersect(mav_samples, sample_names)
ctrl_in_vcf <- intersect(ctrl_samples, sample_names)

# Caricamento del VCF escludendo le righe di intestazione (altrimenti troppo pesante x analisi)
vcf_data <- fread(vcf_file, skip = "#CHROM", header = FALSE, sep = "\t")
setnames(vcf_data, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_names))

# Estrazione dei genotipi
gt_matrix <- vcf_data[, .SD, .SDcols = sample_names]
gt_matrix <- gt_matrix[, lapply(.SD, function(x) substr(x, 1, 3))]

# Matrice logica: TRUE se il paziente è portatore della mutazione (0/1 o 1/1)
is_carrier <- gt_matrix[, lapply(.SD, function(x) x %in% c("0/1", "0|1", "1/0", "1|0", "1/1", "1|1"))]

# Conteggio per gruppo
vcf_data[, Pazienti_MAV := rowSums(is_carrier[, .SD, .SDcols = mav_in_vcf])]
vcf_data[, Pazienti_Controlli := rowSums(is_carrier[, .SD, .SDcols = ctrl_in_vcf])]
vcf_data[, Totale_Portatori := Pazienti_MAV + Pazienti_Controlli]

# Estrazione veloce di Gene e Impatto da INFO
vcf_data[, Impatto := sub(".*ANN=[^|]*\\|([^|]*).*", "\\1", INFO)]
vcf_data[, Gene := sub(".*ANN=[^|]*\\|[^|]*\\|[^|]*\\|([^|]*).*", "\\1", INFO)]
vcf_data[, Modifica := sub(".*ANN=[^|]*\\|[^|]*\\|[^|]*\\|[^|]*\\|[^|]*\\|[^|]*\\|[^|]*\\|[^|]*\\|[^|]*\\|([^|]*).*", "\\1", INFO)] #questo non capisco ma senza non va

# Filtra varianti portate da at least 1 paziente
res <- vcf_data[Totale_Portatori > 0, .(CHROM, POS, REF, ALT, Gene, Impatto, Modifica, Totale_Portatori, Pazienti_MAV, Pazienti_Controlli)]

# Ordina per varianti più frequenti MAV
res <- res[order(-Pazienti_MAV, -Totale_Portatori)]

fwrite(res, "Varianti_High_Moderate_Sottogruppi.csv")

print(head(res, 10))



############################################################################àà
###############àprova filtro per mutazioni principalmente in mav

res <- fread("Varianti_High_Moderate_Sottogruppi.csv")

# ESCLUSIVE MAV: Varianti presenti in almeno 3 pazienti MAV e ZERO controlli
esclusive_mav <- res[Pazienti_MAV >= 3 & Pazienti_Controlli == 0]
esclusive_mav <- esclusive_mav[order(-Pazienti_MAV)]

# ARRICCHITE IN MAV: Rapporto MAV/Controlli molto sbilanciato
arricchite_mav <- res[Pazienti_MAV >= 5 & (Pazienti_MAV / (Pazienti_Controlli + 1)) >= 2]
arricchite_mav <- arricchite_mav[order(-Pazienti_MAV)]

# SOLO IMPATTO "HIGH" (Stop-gain, frameshift, splice site)
impatto_high <- res[Impatto == "HIGH" & Pazienti_MAV >= 2]
impatto_high <- impatto_high[order(-Pazienti_MAV)]

# Salvataggio dei sottoinsiemi utili
fwrite(esclusive_mav, "Varianti_Esclusive_MAV.csv")
fwrite(arricchite_mav, "Varianti_Arricchite_MAV.csv")
fwrite(impatto_high, "Varianti_High_Impact_MAV.csv")

#per vederlo stampato in console che faccio prima
cat("1. Varianti Esclusive nei MAV (>=3 casi, 0 controlli):", nrow(esclusive_mav), "\n")
cat("2. Varianti Arricchite nei MAV (>=5 casi):", nrow(arricchite_mav), "\n")
cat("3. Varianti ad impatto HIGH nei MAV (>=2 casi):", nrow(impatto_high), "\n")

#giga arricchiti
arricchite_mav_giga <- res[Pazienti_MAV >= 6 & (Pazienti_MAV / (Pazienti_Controlli + 1)) >= 2]
arricchite_mav_giga <- arricchite_mav[order(-Pazienti_MAV)]
arricchite_mav_giga$Gene

###################prova per geni sclerosi #quelli di ale del paper introduzione tesi
res <- fread("Varianti_High_Moderate_Sottogruppi.csv")

# Vettore dei geni target pulito
target_genes <- c(
  "ABI3BP", "ADCY3", "AFF1", "ALPK2", "ANAPC1", "ANKRD55", "ARHGAP27", "ASF1A", "ATXN1", 
  "B4GALT5", "BACH2", "BATF", "BATF3", "BCAS1", "BCL10", "BCL6", "BCL9L", "C10orf55", 
  "C1orf106", "C1orf52", "C7orf72", "CAMK2G", "CARD11", "CBLB", "CCDC88B", "CD200R1", 
  "CD226", "CD27-AS1", "CD28", "CD37", "CD40", "CD48", "CD5", "CD58", "CD6", "CD86", 
  "CDHR3", "CHST12", "CLCN6", "CLEC16A", "CLEC2D", "CLECL1", "CMC1", "CNR2", "COPB1", 
  "CSF2RB", "CTLA4", "CTSH", "CXCR4", "CXCR5", "CYP24A1", "DDX6", "DLEU1", "DRAP1", 
  "ELMO1", "ENTHD1", "EOMES", "EPHA4", "EPPK1", "EPS15L1", "ERG", "ESPN", "ETS1", 
  "ETV7", "EVI5", "EXOC6", "EXTL2", "FAM213B", "FAM69A", "FAM76B", "FBXO48", "FCRL2", 
  "FCRL3", "FOXP1", "FRMD6", "FUCA2", "GABARAPL3", "GATA3", "GFI1", "GNG2", "GRAP2", 
  "GRB2", "HBS1L", "HDGFRP2", "HHEX", "HOXA1", "IFI30", "IFITM3", "IFNGR2", "IGSF23", 
  "IKZF1", "IKZF3", "IL12A-AS1", "IL20RA", "IL22RA2", "IL2RA", "IL7R", "ILDR1", "IMPG2", 
  "INTS8", "IQCB1", "IRF5", "IRF8", "JADE2", "JAK1", "JARID2", "JAZF1", "KCNH8", "KEAP1", 
  "KCP", "KPNB1", "LBH", "LCK", "LEF1", "LIMK2", "LINC00271", "LINC00603", "LINC00824", 
  "LINC01082", "LINC01108", "LINC01146", "LOC100129203", "LOC100129603", "LOC100130476", 
  "LOC100505625", "LOC100506457", "LOC101060498", "LOC101929282", "LOC102467146", 
  "LOC102723649", "LOC102723854", "LOC102724874", "LOC152225", "LOC285626", "LOC285740", 
  "LPP", "LPP-AS2", "LRRC34", "LRRIQ4", "LTBR", "L3MBTL3", "LY86", "MAF", "MALT1", 
  "MAP3K14", "MAPK1", "MARCH1", "MERTK", "MIR100HG", "MIR1204", "MIR1208", "MIR4435-1HG", 
  "MIR4778", "MIR548AN", "MPHOSPH9", "MPV17L2", "MTHFR", "MYB", "MYBPC3", "NCF4", 
  "NCOA5", "NDFIP1", "NGRN", "NLRC5", "NPEPPS", "NR1D1", "NSL1", "ODF3B", "OLIG3", 
  "OS9", "PAPD7", "PDE3B", "PDE4A", "PHGDH", "PHLDB1", "PITPNM2", "PKIA", "PLEC", 
  "PLEK", "PLEKHG6", "PLXNC1", "PPA2", "PRDX5", "PRM1", "PRNCR1", "PRR5L", "PTGER4", 
  "PTPRK", "PUS10", "PVR", "PVT1", "RAD51B", "RASGRF1", "RBM17", "RBM43", "RCOR1", 
  "RGS1", "RGS14", "RGS21", "RHOH", "RMI2", "RRAS2", "RREB1", "RUNX3", "SAE1", 
  "SATB1-AS1", "SH2B3", "SLAMF1", "SLAMF7", "SLC2A4RG", "SLC30A7", "SLC6A16", "SLC9A8", 
  "SLC9B1", "SMARCA4", "SOCS1", "SOX8", "SP140", "SPRED2", "SPRY4", "SRRM1", "SSTR5-AS1", 
  "STAT3", "STAT4", "SYDE2", "SYPL1", "TAGAP", "TBX6", "TCF7", "TET2", "THEMIS", 
  "THSD7B", "TIMMDC1", "TNFAIP8", "TNFRSF1A", "TNFSF14", "TNIP3", "TNP2", "TOX2", 
  "TRAF3", "TRIM14", "TRIM71", "TSGA10IP", "TXK", "TYK2", "UBASH3B", "UBXN6", "UPK2", 
  "VANGL2", "VCAM1", "VMP1", "WWOX", "ZC3HAV1", "ZEB1-AS1", "ZBTB38", "ZBTB7B", 
  "ZFP36L1", "ZFP36L2", "ZMIZ1", "ZNF365", "ZNF438", "ZNF767P", "ZNHIT3", "ZPLD1")

# Filtra la tabella
Filtrato_200mut <- res[Gene %in% target_genes]

# Ordina per frequenza 
Filtrato_200mut <- Filtrato_200mut[order(-Pazienti_MAV, -Totale_Portatori)]

fwrite(Filtrato_200mut, "Varianti_Geni_Target_SM.csv")

#nulla di buono


fwrite(esclusive_mav, "esclusive_mav.csv")
fwrite(arricchite_mav, "arricchite_mav.csv")
fwrite(Filtrato_200mut, "Filtrato_200mut.csv")




