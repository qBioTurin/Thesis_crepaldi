source("Libraries.R")

#############################
#calcolo pvalue alpha - PER TP!!
#############################
####batteri
Bacteria_alpha_metadata = read.csv("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\ALPHA\\Bacteria_alpha_metadata.csv")

# Metriche da analizzare
metrics <- c("Observed", "Shannon", "Simpson")

# Timepoint presenti nel dataset
timepoints <- sort(unique(Bacteria_alpha_metadata$tp))

# Tutte le coppie di timepoint
pairs <- combn(timepoints, 2, simplify = FALSE)

# Dataframe dei risultati
all_results_batteri <- map_dfr(metrics, function(metric) {
  
  map_dfr(pairs, function(p) {
    
    tp1 <- p[1]
    tp2 <- p[2]
    
    # Mantiene solo i pazienti con entrambe le misurazioni
    dat <- Bacteria_alpha_metadata %>%
      filter(tp %in% c(tp1, tp2)) %>%
      select(original, tp, value = all_of(metric)) %>%
      pivot_wider(names_from = tp, values_from = value) %>%
      filter(!is.na(.data[[tp1]]) & !is.na(.data[[tp2]]))
    
    # Se ci sono meno di 2 coppie non è possibile fare il test
    if (nrow(dat) < 2) {
      return(data.frame(
        Metric = metric,
        Comparison = paste(tp1, "vs", tp2),
        N = nrow(dat),
        Statistic = NA,
        P_value = NA
      ))
    }
    
    # Wilcoxon signed-rank test paired
    test <- wilcox.test(
      dat[[tp1]],
      dat[[tp2]],
      paired = TRUE,
      exact = FALSE
    )
    
    data.frame(
      Metric = metric,
      Comparison = paste(tp1, "vs", tp2),
      N = nrow(dat),
      Statistic = unname(test$statistic),
      P_value = test$p.value
    )
    
  })
  
})



############archaea
Archaea_alpha_metadata = read.csv("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\ALPHA\\Archaea_alpha_metadata.csv")

# Metriche da analizzare
metrics <- c("Observed", "Shannon", "Simpson")

# Timepoint presenti nel dataset
timepoints <- sort(unique(Archaea_alpha_metadata$tp))

# Tutte le coppie di timepoint
pairs <- combn(timepoints, 2, simplify = FALSE)

# Dataframe dei risultati
all_results_archaea <- map_dfr(metrics, function(metric) {
  
  map_dfr(pairs, function(p) {
    
    tp1 <- p[1]
    tp2 <- p[2]
    
    # Mantiene solo i pazienti con entrambe le misurazioni
    dat <- Archaea_alpha_metadata %>%
      filter(tp %in% c(tp1, tp2)) %>%
      select(original, tp, value = all_of(metric)) %>%
      pivot_wider(names_from = tp, values_from = value) %>%
      filter(!is.na(.data[[tp1]]) & !is.na(.data[[tp2]]))
    
    # Se ci sono meno di 2 coppie non è possibile fare il test
    if (nrow(dat) < 2) {
      return(data.frame(
        Metric = metric,
        Comparison = paste(tp1, "vs", tp2),
        N = nrow(dat),
        Statistic = NA,
        P_value = NA
      ))
    }
    
    # Wilcoxon signed-rank test paired
    test <- wilcox.test(
      dat[[tp1]],
      dat[[tp2]],
      paired = TRUE,
      exact = FALSE
    )
    
    data.frame(
      Metric = metric,
      Comparison = paste(tp1, "vs", tp2),
      N = nrow(dat),
      Statistic = unname(test$statistic),
      P_value = test$p.value
    )
    
  })
  
})

######eukaryota
Eukaryota_alpha_metadata = read.csv("C:\\Users\\gialu\\Desktop\\tesi_magistrale\\Analisi_pts_completi_app\\Risultati_Maaslin\\Risultati_MAV_app\\shared_results\\ALPHA\\Eukaryota_alpha_metadata.csv")

# Metriche da analizzare
metrics <- c("Observed", "Shannon", "Simpson")

# Timepoint presenti nel dataset
timepoints <- sort(unique(Eukaryota_alpha_metadata$tp))

# Tutte le coppie di timepoint
pairs <- combn(timepoints, 2, simplify = FALSE)

# Dataframe dei risultati
all_results_archaea <- map_dfr(metrics, function(metric) {
  
  map_dfr(pairs, function(p) {
    
    tp1 <- p[1]
    tp2 <- p[2]
    
    # Mantiene solo i pazienti con entrambe le misurazioni
    dat <- Eukaryota_alpha_metadata %>%
      filter(tp %in% c(tp1, tp2)) %>%
      select(original, tp, value = all_of(metric)) %>%
      pivot_wider(names_from = tp, values_from = value) %>%
      filter(!is.na(.data[[tp1]]) & !is.na(.data[[tp2]]))
    
    # Se ci sono meno di 2 coppie non è possibile fare il test
    if (nrow(dat) < 2) {
      return(data.frame(
        Metric = metric,
        Comparison = paste(tp1, "vs", tp2),
        N = nrow(dat),
        Statistic = NA,
        P_value = NA
      ))
    }
    
    # Wilcoxon signed-rank test paired
    test <- wilcox.test(
      dat[[tp1]],
      dat[[tp2]],
      paired = TRUE,
      exact = FALSE
    )
    
    data.frame(
      Metric = metric,
      Comparison = paste(tp1, "vs", tp2),
      N = nrow(dat),
      Statistic = unname(test$statistic),
      P_value = test$p.value
    )
    
  })
  
})



