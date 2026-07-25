source("Libraries.R")

rainbow_pastel <- c("#FFADAD",  "#FFD6A5",  "#FDFFB6", "#CAFFBF","#A0C4FF","#CDB4DB" ) 


#########
#batteri
##########
Bacteria_alpha_metadata = read.csv("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\ALPHA\\Bacteria_alpha_metadata.csv")

Bacteria_alpha_metadata$tp = factor(Bacteria_alpha_metadata$tp, levels=c("T0", "T3", "T12", "T15", "T24", "T36"))

#grafico Observed
ggplot(Bacteria_alpha_metadata, aes(x = tp, y = Observed)) +
  geom_boxplot(fill = c("#FFADAD",  "#FFD6A5",  "#FDFFB6", "#CAFFBF","#A0C4FF","#CDB4DB" ), color = "black") +
  labs(
    x = "Timepoint",
    y = "Alpha diversity Observed",
    title = "Bacteria Observed index by timepoint"
  ) +
  theme_classic()

#grafico Shannon
ggplot(Bacteria_alpha_metadata, aes(x = tp, y = Shannon)) +
  geom_boxplot(fill = c("#FFADAD",  "#FFD6A5",  "#FDFFB6", "#CAFFBF","#A0C4FF","#CDB4DB" ), color = "black") +
  labs(
    x = "Timepoint",
    y = "Alpha diversity Shannon",
    title = "Bacteria Shannon index by timepoint"
  ) +
  theme_classic()

#grafico Simpson
ggplot(Bacteria_alpha_metadata, aes(x = tp, y = Simpson)) +
  geom_boxplot(fill = c("#FFADAD",  "#FFD6A5",  "#FDFFB6", "#CAFFBF","#A0C4FF","#CDB4DB" ), color = "black") +
  labs(
    x = "Timepoint",
    y = "Alpha diversity Simpson",
    title = "Bacteria Simpson index by timepoint"
  ) +
  theme_classic()


##########
#archaea
###########
Archaea_alpha_metadata = read.csv("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\ALPHA\\Archaea_alpha_metadata.csv")

Archaea_alpha_metadata$tp = factor(Archaea_alpha_metadata$tp, levels=c("T0", "T3", "T12", "T15", "T24", "T36"))

#grafico Observed
ggplot(Archaea_alpha_metadata, aes(x = tp, y = Observed)) +
  geom_boxplot(fill = c("#FFADAD",  "#FFD6A5",  "#FDFFB6", "#CAFFBF","#A0C4FF","#CDB4DB" ), color = "black") +
  labs(
    x = "Timepoint",
    y = "Alpha diversity Observed",
    title = "Archaea Observed index by timepoint"
  ) +
  theme_classic()

#grafico Shannon
ggplot(Archaea_alpha_metadata, aes(x = tp, y = Shannon)) +
  geom_boxplot(fill = c("#FFADAD",  "#FFD6A5",  "#FDFFB6", "#CAFFBF","#A0C4FF","#CDB4DB" ), color = "black") +
  labs(
    x = "Timepoint",
    y = "Alpha diversity Shannon",
    title = "Archaea Shannon index by timepoint"
  ) +
  theme_classic()

#grafico Simpson
ggplot(Archaea_alpha_metadata, aes(x = tp, y = Simpson)) +
  geom_boxplot(fill = c("#FFADAD",  "#FFD6A5",  "#FDFFB6", "#CAFFBF","#A0C4FF","#CDB4DB" ), color = "black") +
  labs(
    x = "Timepoint",
    y = "Alpha diversity Simpson",
    title = "Archaea Simpson index by timepoint"
  ) +
  theme_classic()

###########
#eukaryota
##########
Eukaryota_alpha_metadata = read.csv("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\ALPHA\\Eukaryota_alpha_metadata.csv")

Eukaryota_alpha_metadata$tp = factor(Eukaryota_alpha_metadata$tp, levels=c("T0", "T3", "T12", "T15", "T24", "T36"))

#grafico Observed
ggplot(Eukaryota_alpha_metadata, aes(x = tp, y = Observed)) +
  geom_boxplot(fill = c("#FFADAD",  "#FFD6A5",  "#FDFFB6", "#CAFFBF","#A0C4FF","#CDB4DB" ), color = "black") +
  labs(
    x = "Timepoint",
    y = "Alpha diversity Observed",
    title = "Eukaryota Observed index by timepoint"
  ) +
  theme_classic()

#grafico Shannon
ggplot(Eukaryota_alpha_metadata, aes(x = tp, y = Shannon)) +
  geom_boxplot(fill = c("#FFADAD",  "#FFD6A5",  "#FDFFB6", "#CAFFBF","#A0C4FF","#CDB4DB" ), color = "black") +
  labs(
    x = "Timepoint",
    y = "Alpha diversity Shannon",
    title = "Eukaryota Shannon index by timepoint"
  ) +
  theme_classic()

#grafico Simpson
ggplot(Eukaryota_alpha_metadata, aes(x = tp, y = Simpson)) +
  geom_boxplot(fill = c("#FFADAD",  "#FFD6A5",  "#FDFFB6", "#CAFFBF","#A0C4FF","#CDB4DB" ), color = "black") +
  labs(
    x = "Timepoint",
    y = "Alpha diversity Simpson",
    title = "Eukaryota Simpson index by timepoint"
  ) +
  theme_classic()
