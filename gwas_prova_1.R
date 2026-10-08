#manattan plot gwas prova iniziale

#file generato da PLINK
gwas <- read.table("risultati_gwas.assoc", header = TRUE)

head(gwas)

#Pulizia dati
gwas$CHR <- as.numeric(as.character(gwas$CHR))
gwas$BP  <- as.numeric(as.character(gwas$BP))
gwas$P   <- as.numeric(as.character(gwas$P))
gwas_clean <- na.omit(gwas)

# Manhattan Plot
png("Manhattan_Plot.png", width = 1000, height = 500, res = 120)
manhattan(gwas_clean, 
          chr = "CHR", 
          bp = "BP", 
          p = "P", 
          snp = "SNP",
          main = "Manhattan Plot - GWAS MAV vs Controlli",
          suggestiveline = -log10(0.05), # treshold significatività (0.05) #prova come paper esplorativi
          genomewideline = FALSE, 
          col = c("blue4", "orange3"),
          ylim = c(0, 3))
dev.off()

#comando più rapido, ma viene meno chiaro visivamente
manhattan(gwas)

##############
##########prova lettura gwas completo senza filtri prima
gwas_results <- read.table("gwas_completo_prova.assoc.adjusted", header = TRUE)

manhattan(gwas_results, chr="CHR", bp="BP", p="P", snp="SNP")
print(manhattan(gwas_results, chr="CHR", bp="BP", p="P", snp="SNP"))

#per salvare
png("Manhattan_Plot_Definitivo.png", width = 1000, height = 500, res = 120)

manhattan(gwas_results)

dev.off()
#cosi viene male perchè senza filtro ci sono troppe varianti e na


#################### prova omettendo na e con un po' di pulizia prima del manhattan
#questo non è nel git ignore quindi metto percorso
setwd("C:/Users/gialu/Desktop/GWAS ANALISI/gwas_analysis")
gwas_big <- fread("gwas_completo_prova.assoc", select = c("CHR", "BP", "SNP", "P"))

gwas_big <- na.omit(gwas_big)

gwas_big[, CHR := as.numeric(as.character(gsub("[^0-9]", "", CHR)))]
gwas_big[, BP  := as.numeric(BP)]
gwas_big[, P   := as.numeric(P)]

# Mantiene solo i p validi e cromosomi riconosciuti
gwas_clean <- gwas_big[!is.na(CHR) & CHR >= 1 & CHR <= 22 & is.finite(P) & P > 0]

gwas_sig <- gwas_clean[P < 0.05]
gwas_nonsig <- gwas_clean[P >= 0.05]

set.seed(123)
# Se lo sfondo è molto grande prende l'1%, altrimenti usa tutte le righe disponibili (modifica per renderlo visibile in immagine)
if (nrow(gwas_nonsig) > 10000) {
  gwas_nonsig_sampled <- gwas_nonsig[sample(.N, .N * 0.01)]
} else {
  gwas_nonsig_sampled <- gwas_nonsig
}

gwas_plot_data <- rbind(gwas_sig, gwas_nonsig_sampled)
gwas_plot_data <- as.data.frame(gwas_plot_data)

#manhattan
png("Manhattan_Plot_GenomeWide.png", width = 1200, height = 600, res = 120)
manhattan(gwas_plot_data, 
          chr = "CHR", 
          bp = "BP", 
          p = "P", 
          snp = "SNP",
          main = "Manhattan Plot - GWAS Genome-Wide (MAV vs Controlli)",
          suggestiveline = -log10(1e-5),
          genomewideline = -log10(5e-8),
          col = c("navy", "skyblue"))
dev.off()




##########varianti più significative dall'analisi senza filtro

gwas_significativo <- gwas_clean[gwas_clean$P <= 1e-5, ]

# Filtra le varianti con P-value significativo dal dataset di PLINK
top_snps <- gwas_clean[P < 1e-5, .(CHR, BP, SNP, P)]

# Ordina dalla più significativa (P-value più basso)
top_snps <- top_snps[order(P)]

print(top_snps)
top_snps <- gwas_big[P < 1e-5, .(CHR, BP, SNP, P)]

