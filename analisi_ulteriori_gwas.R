#analisi GWAS


#################################################################
library(data.table)
library(vegan)

# 1. Carica il file RDS
rds_data <- readRDS("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_app_healty_plus_MAV\\shared_results\\SUPERVISED_DEC\\Bacteria_Supervised_decontam.rds")
library(phyloseq)
library(data.table)
library(vegan)



otu_tab <- as.data.frame(otu_table(rds_data))
tax_tab <- as.data.frame(tax_table(rds_data))
otu_tab_t <- as.data.frame(t(otu_tab))

# 2. Lista completa e precisa dei TaxID identificati
# Includiamo i 6 precedenti + Collinsella (3132705) + le specie di Parabacteroides presenti
target_ids <- c(
  "328813",  # Alistipes onderdonkii
  "84112",   # Eggerthella lenta
  "471189",  # Gordonibacter pamelaea
  "165179",  # Segatella copri
  "116085",  # Coprococcus catus
  "1683",    # Bifidobacterium angulatum
  "3132705" # Collinsella sp. i05-0019-G5
  
)

# Verifica la presenza dei TaxID nella OTU table
valid_ids <- intersect(target_ids, colnames(otu_tab_t))
cat("TaxID finali inseriti nel calcolo:", length(valid_ids), "\n")

# 3. Filtra la OTU table per questi TaxID e calcola Shannon
otu_filtered <- otu_tab_t[, valid_ids, drop = FALSE]
shannon_val <- diversity(otu_filtered, index = "shannon")

df_shannon <- data.table(
  Short_ID = names(shannon_val),
  SHANNON = as.numeric(shannon_val)
)

# 4. Abbina con pheno.txt per ottenere gli ID lunghi del VCF
pheno_orig <- fread("D:\\GWAS_ANALISI\\GWAS\\pheno.txt", header = TRUE)
setnames(pheno_orig, c("FID", "IID", "PHENOTYPE"))

df_shannon[, IID := sapply(Short_ID, function(x) {
  m <- pheno_orig[grepl(x, IID), IID]
  if (length(m) > 0) return(m[1]) else return(NA)
})]

pheno_plink <- merge(pheno_orig[, .(FID, IID)], df_shannon, by = "IID")
pheno_plink <- pheno_plink[, .(FID, IID, SHANNON)]

# 5. Salva il file per PLINK
fwrite(pheno_plink, "pheno_shannon.txt", sep = "\t", col.names = FALSE)

cat("\n=========================================================\n")
cat("File 'pheno_shannon.txt' GENERATO CON SUCCESSO!\n")
cat("=========================================================\n")
print(head(pheno_plink))




#########prova gwas con shannon senza filtro per le specie
library(phyloseq)
library(data.table)
library(vegan)


# 1. Carica l'oggetto RDS phyloseq
rds_data <- readRDS("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_app_healty_plus_MAV\\shared_results\\SUPERVISED_DEC\\Bacteria_Supervised_decontam.rds")

# Estrazione della tabella OTU (matrice campioni x specie)
otu_tab <- as.data.frame(otu_table(rds_data))
otu_tab_t <- as.data.frame(t(otu_tab))

cat("Numero totale di specie/taxa nell'analisi globale:", ncol(otu_tab_t), "\n")

# 2. CALCOLO SHANNON GLOBALE (su TUTTE le specie)
shannon_global <- diversity(otu_tab_t, index = "shannon")

df_shannon_global <- data.table(
  Short_ID = names(shannon_global),
  SHANNON  = as.numeric(shannon_global)
)

# 3. Mappatura degli ID brevi con gli ID complessi del VCF (tramite keep_mav.txt)
keep_mav <- fread("keep_mav.txt", header = FALSE)
setnames(keep_mav, c("FID", "IID"))
keep_mav[, ID_SHORT := sub(".*_(MAV[0-9]+)T0.*", "\\1", IID)]

# Merge degli ID complessi MAV con i valori dello Shannon Globale
mav_shannon_global <- merge(keep_mav[, .(FID, IID, ID_SHORT)], 
                            df_shannon_global, 
                            by.x = "ID_SHORT", 
                            by.y = "Short_ID", 
                            all.x = TRUE)

# ==============================================================================
# OUTPUT 1: FILE PER REGRESSIONE LINEARE (Shannon Continuo Globale)
# ==============================================================================
plink_linear <- mav_shannon_global[, .(FID, IID, SHANNON)]

fwrite(plink_linear, "pheno_shannon_global_linear.txt", sep = "\t", quote = FALSE, col.names = TRUE)
cat("\n1. File 'pheno_shannon_global_linear.txt' salvato correttamente!\n")

# ==============================================================================
# OUTPUT 2: FILE PER REGRESSIONE LOGISTICA (Shannon Binario Globale)
# ==============================================================================
# Calcolo della mediana dello Shannon Globale sui soli pazienti MAV
mediana_global <- median(mav_shannon_global$SHANNON, na.rm = TRUE)
cat("Mediana dello Shannon Globale nei MAV:", mediana_global, "\n")

# Discretizzazione: 1 = Low Shannon (<= mediana), 2 = High Shannon (> mediana)
mav_shannon_global[, PHENOTYPE := ifelse(SHANNON > mediana_global, 2, 1)]
mav_shannon_global[is.na(SHANNON), PHENOTYPE := -9]

plink_binary <- mav_shannon_global[, .(FID, IID, PHENOTYPE)]

fwrite(plink_binary, "pheno_shannon_global_binary.txt", sep = "\t", quote = FALSE, col.names = TRUE)
cat("2. File 'pheno_shannon_global_binary.txt' salvato correttamente!\n")

cat("\n=== CONTEGGI SHANNON GLOBALE BINARIO ===\n")
cat("Pazienti Low Shannon (1):", sum(plink_binary$PHENOTYPE == 1), "\n")
cat("Pazienti High Shannon (2):", sum(plink_binary$PHENOTYPE == 2), "\n")
cat("Dati mancanti (-9):", sum(plink_binary$PHENOTYPE == -9), "\n")