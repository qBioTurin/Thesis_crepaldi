#manattan plot analisi gwas
install.packages("qqman")
library(qqman)
# Imposta la tua cartella di lavoro
setwd("C:/Users/gialu/Desktop/GWAS ANALISI/gwas_analysis")

# Leggi il file generato da PLINK
gwas_results <- read.table("C:/Users/gialu/Desktop/GWAS ANALISI/gwas_analysis/risultati_gwas.assoc", header = TRUE)

# Rimuovi eventuali righe con p-value mancanti (NA)
gwas_results <- na.omit(gwas_results)

# Salva il grafico come immagine PNG
png("Manhattan_Plot.png", width = 1200, height = 600)

# Genera il Manhattan plot
manhattan(gwas_results, 
          chr = "CHR", 
          bp = "BP", 
          p = "P", 
          snp = "SNP",
          main = "Manhattan Plot - GWAS MAV vs Controlli",
          suggestiveline = -log10(1e-5), # Linea blu per significatività suggestiva
          genomewideline = -log10(5e-8), # Linea rossa per significatività genomica
          col = c("blue4", "orange3"))

dev.off()