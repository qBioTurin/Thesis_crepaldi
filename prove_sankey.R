merged_data <- read.csv("Input/merged_data.csv")

#####prova con aurelia per rendere più fluido il sankey
 # Creazione dataset con fasce EDSS
 sankey_data <- merged_data %>%
   transmute(
     ID,
     T0 = case_when(
       EDSS_T0 < 3 ~ "<3",
       EDSS_T0 < 6 ~ "3-5.5",
       EDSS_T0 > 6 ~ "≥6"
     ),
     T6 = case_when(
       EDSS_T6 < 3 ~ "<3",
       EDSS_T6 < 6 ~ "3-5.5",
       EDSS_T6 > 6 ~ "≥6"
     ),
     T12 = case_when(
       EDSS_T12 < 3 ~ "<3",
       EDSS_T12 < 6 ~ "3-5.5",
       EDSS_T12 > 6 ~ "≥6"
     ),
     T18 = case_when(
       EDSS_T18 < 3 ~ "<3",
       EDSS_T18 < 6 ~ "3-5.5",
       EDSS_T18 > 6 ~ "≥6"
     ),
     T24 = case_when(
       EDSS_T24 < 3 ~ "<3",
       EDSS_T24 < 6 ~ "3-5.5",
       EDSS_T24 > 6 ~ "≥6"
     ),
     T30 = case_when(
       EDSS_T30 < 3 ~ "<3",
       EDSS_T30 < 6 ~ "3-5.5",
       EDSS_T30 > 6 ~ "≥6"
     ),
     T36 = case_when(
       EDSS_T36 < 3 ~ "<3",
       EDSS_T36 < 6 ~ "3-5.5",
       EDSS_T36 > 6 ~ "≥6"
     )
   )
 
 
 links <- bind_rows(
   sankey_data %>% filter(!is.na(T0),  !is.na(T6))  %>% transmute(source = paste0("T0_",  T0),  target = paste0("T6_",  T6)),
   sankey_data %>% filter(!is.na(T6),  !is.na(T12)) %>% transmute(source = paste0("T6_",  T6),  target = paste0("T12_", T12)),
   sankey_data %>% filter(!is.na(T12), !is.na(T18)) %>% transmute(source = paste0("T12_", T12), target = paste0("T18_", T18)),
   sankey_data %>% filter(!is.na(T18), !is.na(T24)) %>% transmute(source = paste0("T18_", T18), target = paste0("T24_", T24)),
   sankey_data %>% filter(!is.na(T24), !is.na(T30)) %>% transmute(source = paste0("T24_", T24), target = paste0("T30_", T30)),
   sankey_data %>% filter(!is.na(T30), !is.na(T36)) %>% transmute(source = paste0("T30_", T30), target = paste0("T36_", T36))
 ) %>%
   count(source, target, name = "value")
 
 # Creazione nodi
 nodes <- data.frame(name = unique(c(links$source, links$target))) %>%
   separate(name, into = c("time", "edss"), sep = "_", remove = FALSE) %>%
   mutate(
     group = case_when(
       grepl("≥6", edss)    ~ "≥6",
       grepl("3-5.5", edss) ~ "3-5.5",
       grepl("<3", edss)    ~ "<3"
     ),
     time = factor(time, levels = c("T0","T6","T12","T18","T24","T30","T36")),
     group_order = factor(group, levels = c("≥6","3-5.5","<3"))
   ) %>%
   arrange(time, group_order) %>%
   mutate(name = paste(time, edss, sep = "_")) %>%
   select(name, group)
 
 # Controllo rapido
 print(table(nodes$group))
 
 # Rimappatura indici
 links$source <- match(links$source, nodes$name) - 1
 links$target <- match(links$target, nodes$name) - 1
 
 # Assegna gruppo ai link in base al nodo
 links <- links %>%
   mutate(group = nodes$group[target + 1])
 
 # Sankey plot
 sankeyNetwork(
   Links = links,
   Nodes = nodes,
   Source = "source",
   Target = "target",
   Value = "value",
   NodeID = "name",
   NodeGroup = "group",
   LinkGroup = "group",
   fontSize = 12,
   nodeWidth = 25,
   colourScale = JS("d3.scaleOrdinal().domain(['≥6','3-5.5','<3']).range(['#E74C3C','#F1C40F','#2ECC71'])"),
   nodePadding = 20,
   iterations = 0,
   sinksRight = FALSE
 )
 
 
 sankey_data %>%
   filter(if_any(T0:T36, is.na)) %>%
   nrow()
 
 
 # Creazione dataset con fasce EDSS
 sankey_data <- merged_data %>%
   transmute(
     ID,
     T0 = case_when(
       EDSS_T0 < 3 ~ "<3",
       EDSS_T0 < 6 ~ "3-5.5",
       EDSS_T0 >= 6 ~ "≥6"
     ),
     T6 = case_when(
       EDSS_T6 < 3 ~ "<3",
       EDSS_T6 < 6 ~ "3-5.5",
       EDSS_T6 >= 6 ~ "≥6"
     ),
     T12 = case_when(
       EDSS_T12 < 3 ~ "<3",
       EDSS_T12 < 6 ~ "3-5.5",
       EDSS_T12 >= 6 ~ "≥6"
     ),
     T18 = case_when(
       EDSS_T18 < 3 ~ "<3",
       EDSS_T18 < 6 ~ "3-5.5",
       EDSS_T18 >= 6 ~ "≥6"
     ),
     T24 = case_when(
       EDSS_T24 < 3 ~ "<3",
       EDSS_T24 < 6 ~ "3-5.5",
       EDSS_T24 >= 6 ~ "≥6"
     ),
     T30 = case_when(
       EDSS_T30 < 3 ~ "<3",
       EDSS_T30 < 6 ~ "3-5.5",
       EDSS_T30 >= 6 ~ "≥6"
     ),
     T36 = case_when(
       EDSS_T36 < 3 ~ "<3",
       EDSS_T36 < 6 ~ "3-5.5",
       EDSS_T36 >= 6 ~ "≥6"
     )
   ) %>%
   # porta avanti l'ultima fascia osservata dove c'è NA
   pivot_longer(-ID, names_to = "time", values_to = "band") %>%
   mutate(time = factor(time, levels = c("T0","T6","T12","T18","T24","T30","T36"))) %>%
   arrange(ID, time) %>%
   group_by(ID) %>%
   fill(band, .direction = "down") %>%
   ungroup() %>%
   pivot_wider(names_from = time, values_from = band) %>%
   filter(!is.na(T0))   # esclude chi non ha EDSS basale (niente da portare avanti)
 
 links <- bind_rows(
   sankey_data %>% filter(!is.na(T0),  !is.na(T6))  %>% transmute(source = paste0("T0_",  T0),  target = paste0("T6_",  T6)),
   sankey_data %>% filter(!is.na(T6),  !is.na(T12)) %>% transmute(source = paste0("T6_",  T6),  target = paste0("T12_", T12)),
   sankey_data %>% filter(!is.na(T12), !is.na(T18)) %>% transmute(source = paste0("T12_", T12), target = paste0("T18_", T18)),
   sankey_data %>% filter(!is.na(T18), !is.na(T24)) %>% transmute(source = paste0("T18_", T18), target = paste0("T24_", T24)),
   sankey_data %>% filter(!is.na(T24), !is.na(T30)) %>% transmute(source = paste0("T24_", T24), target = paste0("T30_", T30)),
   sankey_data %>% filter(!is.na(T30), !is.na(T36)) %>% transmute(source = paste0("T30_", T30), target = paste0("T36_", T36))
 ) %>%
   count(source, target, name = "value")
 
 # Creazione nodi
 nodes <- data.frame(name = unique(c(links$source, links$target))) %>%
   separate(name, into = c("time", "edss"), sep = "_", remove = FALSE) %>%
   mutate(
     group = case_when(
       grepl("≥6", edss)    ~ "≥6",
       grepl("3-5.5", edss) ~ "3-5.5",
       grepl("<3", edss)    ~ "<3"
     ),
     time = factor(time, levels = c("T0","T6","T12","T18","T24","T30","T36")),
     group_order = factor(group, levels = c("≥6","3-5.5","<3"))
   ) %>%
   arrange(time, group_order) %>%
   mutate(name = paste(time, edss, sep = "_")) %>%
   select(name, group)
 
 # Controllo rapido
 print(table(nodes$group))
 
 # Rimappatura indici JS
 links$source <- match(links$source, nodes$name) - 1
 links$target <- match(links$target, nodes$name) - 1
 
 # Assegna gruppo ai link in base al nodo di destinazione
 links <- links %>%
   mutate(group = nodes$group[target + 1])
 
 # Sankey plot
 sankeyNetwork(
   Links = links,
   Nodes = nodes,
   Source = "source",
   Target = "target",
   Value = "value",
   NodeID = "name",
   NodeGroup = "group",
   LinkGroup = "group",
   fontSize = 12,
   nodeWidth = 25,
   colourScale = JS("d3.scaleOrdinal().domain(['≥6','3-5.5','<3']).range(['#E74C3C','#F1C40F','#2ECC71'])"),
   nodePadding = 20,
   iterations = 0,
   sinksRight = FALSE
 )
 
 
 
 
 
 #DEFINITIVO CON I PAZIENTI CHE SI FERMANO ALL'NA
 ###################################################
 # Fasce EDSS troncate al primo NA
 sankey_long <- merged_data %>%
   transmute(
     ID,
     T0 = case_when(
       EDSS_T0 < 3 ~ "<3",
       EDSS_T0 < 6 ~ "3-5.5",
       EDSS_T0 >= 6 ~ "≥6"
     ),
     T6 = case_when(
       EDSS_T6 < 3 ~ "<3",
       EDSS_T6 < 6 ~ "3-5.5",
       EDSS_T6 >= 6 ~ "≥6"
     ),
     T12 = case_when(
       EDSS_T12 < 3 ~ "<3",
       EDSS_T12 < 6 ~ "3-5.5",
       EDSS_T12 >= 6 ~ "≥6"
     ),
     T18 = case_when(
       EDSS_T18 < 3 ~ "<3",
       EDSS_T18 < 6 ~ "3-5.5",
       EDSS_T18 >= 6 ~ "≥6"
     ),
     T24 = case_when(
       EDSS_T24 < 3 ~ "<3",
       EDSS_T24 < 6 ~ "3-5.5",
       EDSS_T24 >= 6 ~ "≥6"
     ),
     T30 = case_when(
       EDSS_T30 < 3 ~ "<3",
       EDSS_T30 < 6 ~ "3-5.5",
       EDSS_T30 >= 6 ~ "≥6"
     ),
     T36 = case_when(
       EDSS_T36 < 3 ~ "<3",
       EDSS_T36 < 6 ~ "3-5.5",
       EDSS_T36 >= 6 ~ "≥6"
     )
   ) %>%
   pivot_longer(-ID, names_to = "time", values_to = "band") %>%
   mutate(time = factor(time, levels = c("T0","T6","T12","T18","T24","T30","T36"))) %>%
   arrange(ID, time) %>%
   group_by(ID) %>%
   filter(cumsum(is.na(band)) == 0) %>%   # tiene le osservazioni fino al primo NA (escluso)
   ungroup()
 
 # Link tra timepoint consecutivi
 links <- sankey_long %>%
   group_by(ID) %>%
   mutate(
     source = paste0(time, "_", band),
     target = lead(paste0(time, "_", band))
   ) %>%
   ungroup() %>%
   filter(!is.na(target)) %>%        # l'ultima osservazione di ogni paziente non genera link
   count(source, target, name = "value")
 
 
 # Creazione nodi
 nodes <- data.frame(name = unique(c(links$source, links$target))) %>%
   separate(name, into = c("time", "edss"), sep = "_", remove = FALSE) %>%
   mutate(
     group = case_when(
       grepl("≥6", edss)    ~ "≥6",
       grepl("3-5.5", edss) ~ "3-5.5",
       grepl("<3", edss)    ~ "<3"
     ),
     time = factor(time, levels = c("T0","T6","T12","T18","T24","T30","T36")),
     group_order = factor(group, levels = c("≥6","3-5.5","<3"))
   ) %>%
   arrange(time, group_order) %>%
   mutate(name = paste(time, edss, sep = "_")) %>%
   select(name, group)
 
 # Controllo rapido
 print(table(nodes$group))
 
 # Rimappatura indici JS
 links$source <- match(links$source, nodes$name) - 1
 links$target <- match(links$target, nodes$name) - 1
 
 # Assegna gruppo ai link in base al nodo di destinazione
 links <- links %>%
   mutate(group = nodes$group[target + 1])
 
 
 # Sankey plot
 sankeyNetwork(
   Links = links,
   Nodes = nodes,
   Source = "source",
   Target = "target",
   Value = "value",
   NodeID = "name",
   NodeGroup = "group",
   LinkGroup = "group",
   fontSize = 12,
   nodeWidth = 25,
   colourScale = JS("d3.scaleOrdinal().domain(['≥6','3-5.5','<3']).range(['#E74C3C','#F1C40F','#2ECC71'])"),
   nodePadding = 20,
   iterations = 0,
   sinksRight = FALSE
 )
 
 
 
 
 # Numero di pazienti totali nel dataset
 n_totali <- merged_data %>%
   summarise(n = n_distinct(ID)) %>%
   pull(n)
 
 # Numero di pazienti considerati dal Sankey
 n_sankey <- sankey_long %>%
   summarise(n = n_distinct(ID)) %>%
   pull(n)
 
 cat("Pazienti totali in merged_data:", n_totali, "\n")
 cat("Pazienti considerati dal Sankey:", n_sankey, "\n")
 cat("Pazienti esclusi completamente dal Sankey:", n_totali - n_sankey, "\n")
 