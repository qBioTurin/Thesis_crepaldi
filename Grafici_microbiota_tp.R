source("Libraries.R")

##########
#batteri
###########
Bacteria_alpha_metadata <- read.csv("C:/Users/gialu/Desktop/tesi_magistrale/Analisi_pts_completi_app/Risultati_Maaslin/Risultati_MAV_app/shared_results/ALPHA/Bacteria_alpha_metadata.csv")

Bacteria_alpha_metadata$tp <- factor(Bacteria_alpha_metadata$tp, levels = c("T0","T3","T12","T15","T24","T36"))

metrics <- c("Observed","Shannon","Simpson")


# Tutte le coppie di time point
timepoints <- as.character(levels(Bacteria_alpha_metadata$tp))
timepoints <- timepoints[timepoints %in% unique(as.character(Bacteria_alpha_metadata$tp))]

pairs <- combn(timepoints,2,simplify = FALSE)

# Wilcoxon paired
all_results_batteri <- map_dfr(metrics,function(metric){
  map_dfr(pairs,function(p){
    tp1 <- as.character(p[[1]])
    tp2 <- as.character(p[[2]])
    dat <- Bacteria_alpha_metadata %>%
      filter(tp %in% c(tp1,tp2)) %>%
      select(original,tp,value = all_of(metric)) %>%
      pivot_wider(names_from = tp,
                  values_from = value) %>%
      filter(!is.na(.data[[tp1]]),
             !is.na(.data[[tp2]]))
    if(nrow(dat) < 2){
      return(
        tibble(
          Metric = metric,
          group1 = tp1,
          group2 = tp2,
          N = nrow(dat),
          Statistic = NA,
          P_value = NA))}
    
test <- wilcox.test(dat[[tp1]], dat[[tp2]], paired = TRUE, exact = FALSE )
    
tibble(Metric = metric,
      group1 = tp1,
      group2 = tp2,
      N = nrow(dat),
      Statistic = unname(test$statistic),
      P_value = test$p.value)})
  })

# Funzione per fare il grafico
plot_alpha <- function(metric){
  stat <- all_results_batteri %>%
    filter(Metric == metric,
           !is.na(P_value),
           P_value < 0.05)
ymax <- max(Bacteria_alpha_metadata[[metric]], na.rm = TRUE)
  if(nrow(stat)>0){
    stat$y.position <- seq(
      ymax*1.05,
      ymax*1.35,
      length.out = nrow(stat) )
    
stat$p.signif <- case_when(
      stat$P_value <= 0.0001 ~ "****",
      stat$P_value <= 0.001 ~ "***",
      stat$P_value <= 0.01 ~ "**",
      stat$P_value < 0.05 ~ "*" )
  
stat$label <- paste0(stat$p.signif, " (p = ", formatC(stat$P_value, digits = 3, format = "f"), ")" )}
  
  p <- ggplot(Bacteria_alpha_metadata, aes(x = tp, y = .data[[metric]])) +
    geom_boxplot(fill = c("#FFADAD",  "#FFD6A5",  "#FDFFB6", "#CAFFBF","#A0C4FF","#CDB4DB" ), color = "black") +
    theme_classic(base_size = 14) +
    labs(
      x = "Timepoint",
      y = ifelse(metric == "Simpson", "Gini-Simpson", metric))
  
  if(nrow(stat)>0){ p <- p + stat_pvalue_manual( stat, label = "label", tip.length = 0.01)}
  return(p)}

bacteria_observed=plot_alpha("Observed")
bacteria_shannon=plot_alpha("Shannon")
bacteria_simpson=plot_alpha("Simpson")

batteri_composita_alpha <- plot_grid(bacteria_observed, bacteria_shannon, bacteria_simpson, ncol = 3)
ggsave("batteri_composita_alpha.jpeg", batteri_composita_alpha, width = 12, height = 8, dpi = 300)




##########
#archaea
###########
Archaea_alpha_metadata = read.csv("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\ALPHA\\Archaea_alpha_metadata.csv")

Archaea_alpha_metadata$tp = factor(Archaea_alpha_metadata$tp, levels=c("T0", "T3", "T12", "T15", "T24", "T36"))

metrics <- c("Observed","Shannon","Simpson")


# Tutte le coppie di time point
timepoints <- as.character(levels(Archaea_alpha_metadata$tp))
timepoints <- timepoints[timepoints %in% unique(as.character(Archaea_alpha_metadata$tp))]

pairs <- combn(timepoints,2,simplify = FALSE)

# Wilcoxon paired
all_results_archaea <- map_dfr(metrics,function(metric){
  map_dfr(pairs,function(p){
    tp1 <- as.character(p[[1]])
    tp2 <- as.character(p[[2]])
    dat <- Archaea_alpha_metadata %>%
      filter(tp %in% c(tp1,tp2)) %>%
      select(original,tp,value = all_of(metric)) %>%
      pivot_wider(names_from = tp,
                  values_from = value) %>%
      filter(!is.na(.data[[tp1]]),
             !is.na(.data[[tp2]]))
    if(nrow(dat) < 2){ return(tibble(
          Metric = metric,
          group1 = tp1,
          group2 = tp2,
          N = nrow(dat),
          Statistic = NA,
          P_value = NA))}
    
test <- wilcox.test(dat[[tp1]], dat[[tp2]], paired = TRUE, exact = FALSE )
    
tibble(Metric = metric,
      group1 = tp1,
      group2 = tp2,
      N = nrow(dat),
      Statistic = unname(test$statistic),
      P_value = test$p.value)})
})

# Funzione per fare il grafico
plot_alpha <- function(metric){
  stat <- all_results_archaea %>%
    filter(Metric == metric,
           !is.na(P_value),
           P_value < 0.05)
ymax <- max(Archaea_alpha_metadata[[metric]], na.rm = TRUE)
  if(nrow(stat)>0){
    stat$y.position <- seq(
      ymax*1.05,
      ymax*1.35,
      length.out = nrow(stat) )
    
stat$p.signif <- case_when(
      stat$P_value <= 0.0001 ~ "****",
      stat$P_value <= 0.001 ~ "***",
      stat$P_value <= 0.01 ~ "**",
      stat$P_value < 0.05 ~ "*" )
    
stat$label <- paste0(stat$p.signif,  " (p = ", formatC(stat$P_value, digits = 3, format = "f"), ")" )}
  
p <- ggplot(Archaea_alpha_metadata, aes(x = tp, y = .data[[metric]])) +
  geom_boxplot(fill = c("#FFADAD",  "#FFD6A5",  "#FDFFB6", "#CAFFBF","#A0C4FF","#CDB4DB" ), color = "black") +
    theme_classic(base_size = 14) +
    labs(
      x = "Timepoint",
      y = ifelse(metric == "Simpson", "Gini-Simpson", metric)
    )
  
if(nrow(stat)>0){ p <- p + stat_pvalue_manual( stat, label = "label", tip.length = 0.01)}
  return(p)}

archaea_observed = plot_alpha("Observed")
archaea_shannon = plot_alpha("Shannon")
archaea_simpson = plot_alpha("Simpson")

archaea_composita_alpha <- plot_grid(archaea_observed, archaea_shannon, archaea_simpson, ncol = 3)
ggsave("archaea_composita_alpha.jpeg", archaea_composita_alpha, width = 12, height = 8, dpi = 300)


###########
#eukaryota
##########
Eukaryota_alpha_metadata = read.csv("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\ALPHA\\Eukaryota_alpha_metadata.csv")

Eukaryota_alpha_metadata$tp = factor(Eukaryota_alpha_metadata$tp, levels=c("T0", "T3", "T12", "T15", "T24", "T36"))

metrics <- c("Observed","Shannon","Simpson")


# Tutte le coppie di time point
timepoints <- as.character(levels(Eukaryota_alpha_metadata$tp))
timepoints <- timepoints[timepoints %in% unique(as.character(Eukaryota_alpha_metadata$tp))]

pairs <- combn(timepoints,2,simplify = FALSE)

# Wilcoxon paired
all_results_eukaryota <- map_dfr(metrics,function(metric){
  map_dfr(pairs,function(p){
    tp1 <- as.character(p[[1]])
    tp2 <- as.character(p[[2]])
    dat <- Eukaryota_alpha_metadata %>%
      filter(tp %in% c(tp1,tp2)) %>%
      select(original,tp,value = all_of(metric)) %>%
      pivot_wider(names_from = tp,
                  values_from = value) %>%
      filter(!is.na(.data[[tp1]]),
             !is.na(.data[[tp2]]))
    if(nrow(dat) < 2){
      return(
        tibble(
          Metric = metric,
          group1 = tp1,
          group2 = tp2,
          N = nrow(dat),
          Statistic = NA,
          P_value = NA))}
    
test <- wilcox.test(dat[[tp1]], dat[[tp2]], paired = TRUE, exact = FALSE )
    
tibble(Metric = metric,
      group1 = tp1,
      group2 = tp2,
      N = nrow(dat),
      Statistic = unname(test$statistic),
      P_value = test$p.value)})
})

# Funzione per fare il grafico
plot_alpha <- function(metric){
  stat <- all_results_eukaryota %>%
    filter(Metric == metric,
           !is.na(P_value),
           P_value < 0.05)
  ymax <- max(Eukaryota_alpha_metadata[[metric]], na.rm = TRUE)
  if(nrow(stat)>0){
    stat$y.position <- seq(
      ymax*1.05,
      ymax*1.35,
      length.out = nrow(stat) )
    
    stat$p.signif <- case_when(
      stat$P_value <= 0.0001 ~ "****",
      stat$P_value <= 0.001 ~ "***",
      stat$P_value <= 0.01 ~ "**",
      stat$P_value < 0.05 ~ "*" )
    
stat$label <- paste0(stat$p.signif,  " (p = ", formatC(stat$P_value, digits = 3, format = "f"), ")" )}
  
p <- ggplot(Eukaryota_alpha_metadata, aes(x = tp, y = .data[[metric]])) +
    geom_boxplot(fill = c("#FFADAD",  "#FFD6A5",  "#FDFFB6", "#CAFFBF","#A0C4FF","#CDB4DB" ), color = "black") +
    theme_classic(base_size = 14) +
    labs(
      x = "Timepoint",
      y = ifelse(metric == "Simpson", "Gini-Simpson", metric)
    )
  
if(nrow(stat)>0){ p <- p + stat_pvalue_manual( stat, label = "label", tip.length = 0.01)}
  return(p)}

eukaryoti_observed = plot_alpha("Observed")
eukaryoti_shannon = plot_alpha("Shannon")
eukaryoti_simpson = plot_alpha("Simpson")

eukaryota_composita_alpha <- plot_grid(eukaryoti_observed, eukaryoti_shannon, eukaryoti_simpson, ncol = 3)
ggsave("eukaryota_composita_alpha.jpeg", eukaryota_composita_alpha, width = 12, height = 8, dpi = 300)






######################
######calcolo IQR SANI VS MALATI DA INSERIRE IN TESI X OBSERVED
library(dplyr)
Bacteria_alpha_metadata <- read.csv("C:/Users/gialu/Desktop/tesi_magistrale/Analisi_pts_completi_app/Risultati_Maaslin/Risultati_app_healty_plus_MAV/shared_results/ALPHA/Bacteria_alpha_metadata.csv")
# MAV = ID che iniziano con "MAV"
mav_t0 <- Bacteria_alpha_metadata %>%
  filter(grepl("^MAV", id))

# HD = tutti gli altri
hd_t0 <- Bacteria_alpha_metadata %>%
  filter(!grepl("^MAV", id))

# Statistiche MAV
mav_stats <- mav_t0 %>%
  summarise(
    median_observed = median(Observed, na.rm = TRUE),
    IQR_observed = IQR(Observed, na.rm = TRUE),
    Q1 = quantile(Observed, 0.25, na.rm = TRUE),
    Q3 = quantile(Observed, 0.75, na.rm = TRUE),
    n = n())

# Statistiche HD
hd_stats <- hd_t0 %>%
  summarise(
    median_observed = median(Observed, na.rm = TRUE),
    IQR_observed = IQR(Observed, na.rm = TRUE),
    Q1 = quantile(Observed, 0.25, na.rm = TRUE),
    Q3 = quantile(Observed, 0.75, na.rm = TRUE),
    n = n())

mav_stats
hd_stats
