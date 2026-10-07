source("Libraries.R")
#filtro il dataset da mettere sull'app che deve corrispondere all'id dei sequenziamenti
modified_metadata <- read.csv("Input/merged_data.csv")

batch_data <- readxl::read_excel("Input/batch_sequenziamento.xlsx",skip = 1) %>%
  select(ID, sequencing_batch) %>%
  rename(id = ID)
batch_data$id[batch_data$id == "MAV5T15"] <- "MAV05T15"
#rinomino MAV23 con i report giusti
batch_data$id[batch_data$id == "MAV23T12"] <- "MAV23T0"
batch_data$id[batch_data$id == "MAV23T24"] <- "MAV23T12"
batch_data$id[batch_data$id == "MAV23T36"] <- "MAV23T24"

modified_metadata <- modified_metadata %>%
  filter(!ID %in% c("MAV08", "MAV24"))

modified_metadata <- modified_metadata %>%
  pivot_longer(
    cols = matches("^(EDSS|lesioni|relapse|cortisone)_"),
    names_to = c(".value", "tp"),
    names_sep = "_"
  ) %>%
  filter(tp %in% c("T0", "T6", "T12", "T18", "T24", "T36")) %>%
  bind_rows(
    filter(., tp == "T6") %>% mutate(tp = "T3"),
    filter(., tp == "T18") %>% mutate(tp = "T15")
  ) %>%
  filter(!tp %in% c("T6", "T18")) %>%   # rimuove T6 e T18
  mutate(tp = factor(tp, levels = c("T0", "T3", "T12", "T15", "T24", "T36"))) %>%
  arrange(ID, tp) %>%
  group_by(ID) %>%
  mutate(
    EDSS_trend = case_when(
      tp == "T3" ~ case_when(
        EDSS[tp == "T12"] == EDSS[tp == "T0"] ~ "stable",
        EDSS[tp == "T12"] > EDSS[tp == "T0"]  ~ "worsening",
        EDSS[tp == "T12"] < EDSS[tp == "T0"]  ~ "improvement",
        TRUE ~ NA_character_
      ),
      tp == "T15" ~ case_when(
        EDSS[tp == "T24"] == EDSS[tp == "T12"] ~ "stable",
        EDSS[tp == "T24"] > EDSS[tp == "T12"]  ~ "worsening",
        EDSS[tp == "T24"] < EDSS[tp == "T12"]  ~ "improvement",
        TRUE ~ NA_character_
      ),
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  mutate(sample = paste0(ID, tp)) %>%
  rename(original = ID) %>%
  rename(id = sample) %>%
  select(id, everything()) %>%
  select(where(~ !all(is.na(.)))) %>%
  select(-sequencing_batch) %>%
  left_join(batch_data, by = "id") %>%
  filter(!is.na(sequencing_batch))



#mettere in ordine nei grafici app aggiungendo lo 0 perchè non supporta la decisione dell'ordine
modified_metadata$id <- replace(
  modified_metadata$id,
  modified_metadata$id %in% c("MAV05T3","MAV07T3","MAV09T3","MAV11T3","MAV13T3","MAV14T3","MAV17T3","MAV18T3","MAV23T3","MAV26T3","MAV27T3","MAV28T3","MAV29T3","MAV30T3","MAV31T3","MAV32T3","MAV34T3"),
  c("MAV05T03","MAV07T03","MAV09T03","MAV11T03","MAV13T03","MAV14T03","MAV17T03","MAV18T03","MAV23T03","MAV26T03","MAV27T03","MAV28T03","MAV29T03","MAV30T03","MAV31T03","MAV32T03","MAV34T03")[match(modified_metadata$id[modified_metadata$id %in% c("MAV05T3","MAV07T3","MAV09T3","MAV11T3","MAV13T3","MAV14T3","MAV17T3","MAV18T3","MAV23T3","MAV26T3","MAV27T3","MAV28T3","MAV29T3","MAV30T3","MAV31T3","MAV32T3","MAV34T3")], c("MAV05T3","MAV07T3","MAV09T3","MAV11T3","MAV13T3","MAV14T3","MAV17T3","MAV18T3","MAV23T3","MAV26T3","MAV27T3","MAV28T3","MAV29T3","MAV30T3","MAV31T3","MAV32T3","MAV34T3"))]
)

modified_metadata <- modified_metadata %>%
  mutate(
    medas_group = if_else(MEDAS >= 46, "High_adherence", "Low_adherence"),
    bmi_group = if_else(bmi >= 24, "High_bmi", "Low_bmi"))

###scrivi file csv come input per app
write.csv(modified_metadata, "C:/Users/gialu/Desktop/tesi_magistrale/Git.hub/Thesis_crepaldi/Input/modified_metadata.csv", row.names = FALSE)


############
#creazione dataset con t3 e t15 uniti per fare la beta diversity di prova
######
modified_metadata <- read.csv("Input/merged_data.csv")

batch_data <- readxl::read_excel("Input/batch_sequenziamento.xlsx", skip = 1) %>%
  select(ID, sequencing_batch) %>%
  rename(id = ID)

batch_data$id[batch_data$id == "MAV5T15"] <- "MAV05T15"
#rinomino MAV23 con i report giusti
batch_data$id[batch_data$id == "MAV23T12"] <- "MAV23T0"
batch_data$id[batch_data$id == "MAV23T24"] <- "MAV23T12"
batch_data$id[batch_data$id == "MAV23T36"] <- "MAV23T24"
modified_metadata <- modified_metadata %>%
  filter(!ID %in% c("MAV08", "MAV24"))

modified_metadata <- modified_metadata %>%
  pivot_longer(
    cols = matches("^(EDSS|lesioni|relapse|cortisone)_"),
    names_to = c(".value", "tp"),
    names_sep = "_"
  ) %>%
  filter(tp %in% c("T0", "T6", "T12", "T18", "T24", "T36")) %>%
  bind_rows(
    filter(., tp == "T6") %>% mutate(tp = "T3"),
    filter(., tp == "T18") %>% mutate(tp = "T15")
  ) %>%
  filter(!tp %in% c("T6", "T18")) %>%   # rimuove T6 e T18
  mutate(tp = factor(tp, levels = c("T0", "T3", "T12", "T15", "T24", "T36"))) %>%
  arrange(ID, tp) %>%
  group_by(ID) %>%
  mutate(
    EDSS_trend = case_when(
      tp == "T3" ~ case_when(
        EDSS[tp == "T12"] == EDSS[tp == "T0"] ~ "stable",
        EDSS[tp == "T12"] > EDSS[tp == "T0"]  ~ "worsening",
        EDSS[tp == "T12"] < EDSS[tp == "T0"]  ~ "improvement",
        TRUE ~ NA_character_
      ),
      tp == "T15" ~ case_when(
        EDSS[tp == "T24"] == EDSS[tp == "T12"] ~ "stable",
        EDSS[tp == "T24"] > EDSS[tp == "T12"]  ~ "worsening",
        EDSS[tp == "T24"] < EDSS[tp == "T12"]  ~ "improvement",
        TRUE ~ NA_character_
      ),
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  mutate(sample = paste0(ID, tp)) %>%          # id "congela" qui il tp originale (T3/T15)
  mutate(
    tp = case_when(
      tp %in% c("T3", "T15") ~ "Post_treatment",
      TRUE ~ as.character(tp)
    ),
    tp = factor(tp, levels = c("T0", "Post_treatment", "T12", "T24", "T36"))
  ) %>%
  rename(original = ID) %>%
  rename(id = sample) %>%
  select(id, everything()) %>%
  select(where(~ !all(is.na(.)))) %>%
  select(-sequencing_batch) %>%
  left_join(batch_data, by = "id") %>%
  filter(!is.na(sequencing_batch))

#mettere in ordine nei grafici app
modified_metadata$id <- replace(
  modified_metadata$id,
  modified_metadata$id %in% c("MAV05T3","MAV07T3","MAV09T3","MAV11T3","MAV13T3","MAV14T3","MAV17T3","MAV18T3","MAV23T3","MAV26T3","MAV27T3","MAV28T3","MAV29T3","MAV30T3","MAV31T3","MAV32T3","MAV34T3"),
  c("MAV05T03","MAV07T03","MAV09T03","MAV11T03","MAV13T03","MAV14T03","MAV17T03","MAV18T03","MAV23T03","MAV26T03","MAV27T03","MAV28T03","MAV29T03","MAV30T03","MAV31T03","MAV32T03","MAV34T03")[match(modified_metadata$id[modified_metadata$id %in% c("MAV05T3","MAV07T3","MAV09T3","MAV11T3","MAV13T3","MAV14T3","MAV17T3","MAV18T3","MAV23T3","MAV26T3","MAV27T3","MAV28T3","MAV29T3","MAV30T3","MAV31T3","MAV32T3","MAV34T3")], c("MAV05T3","MAV07T3","MAV09T3","MAV11T3","MAV13T3","MAV14T3","MAV17T3","MAV18T3","MAV23T3","MAV26T3","MAV27T3","MAV28T3","MAV29T3","MAV30T3","MAV31T3","MAV32T3","MAV34T3"))]
)

modified_metadata <- modified_metadata %>%
  mutate(
    medas_group = if_else(MEDAS >= 46, "High_adherence", "Low_adherence"),
    bmi_group = if_else(bmi >= 24, "High_bmi", "Low_bmi")
  )
###scrivi file csv per beta diversity
write.csv(modified_metadata, "C:/Users/gialu/Desktop/tesi_magistrale/Git.hub/Thesis_crepaldi/Input/data_x_beta.csv", row.names = FALSE)




##-----------------------------------------------##
######dataset con anche i sani per MAV vs HD############

merged_data <- read.csv("Input/merged_data.csv")

batch_data <- readxl::read_excel("Input/batch_sequenziamento.xlsx", skip = 1) 
batch_data$ID[batch_data$ID == "MAV5T15"] <- "MAV05T15"
#rinomino MAV23 con i report giusti
batch_data$ID[batch_data$ID == "MAV23T12"] <- "MAV23T0"
batch_data$ID[batch_data$ID == "MAV23T24"] <- "MAV23T12"
batch_data$ID[batch_data$ID == "MAV23T36"] <- "MAV23T24"


batch_data <- batch_data %>%
select(ID, sequencing_batch) %>%
  rename(id = ID) %>%
  filter(grepl("^MAV.*T0$", id)) %>%
  mutate(
    id = sub("T0$", "", id))

cleaned_metadata <- read.csv(
  "Input/metadata_cleaned_completo.csv",
  sep=";"
)

cleaned_metadata_SANI <- cleaned_metadata %>%
  filter(grepl("^(HV|AMCTRL|MBCTRL|NBCTRL|RRCTRL)", ID))

names(cleaned_metadata_SANI)[names(cleaned_metadata_SANI)=="data_nascita"] <- "Data_di_nascita"
names(cleaned_metadata_SANI)[names(cleaned_metadata_SANI)=="age_at_T0"] <- "age_at_t0"

cols_keep <- c(
  "ID",
  "Data_di_nascita",
  "naive",
  "sex",
  "bmi",
  "age_at_t0",
  "sequencing_batch")

merged_data_reduced <- merged_data %>%
  select(any_of(cols_keep))

cleaned_metadata_SANI_reduced <- cleaned_metadata_SANI %>%
  select(any_of(cols_keep))

merged_data_complete <- bind_rows(
  merged_data_reduced,
  cleaned_metadata_SANI_reduced
) %>%
  rename(id = ID)

## aggiorna SOLO i MAV
merged_data_complete <- merged_data_complete %>%
  left_join(
    batch_data %>%
      rename(sequencing_batch_T0 = sequencing_batch),
    by = "id"
  ) %>%
  mutate(
    sequencing_batch = if_else(
      startsWith(id, "MAV"),
      sequencing_batch_T0,
      sequencing_batch
    )
  ) %>%
  select(-sequencing_batch_T0)

## elimina i MAV senza T0
merged_data_complete <- merged_data_complete %>%
  filter(!id %in% c(
    "MAV08",
    "MAV12",
    "MAV16",
    "MAV22",
    "MAV24",
    "MAV25",
    "MAV31"))

## aggiungi T0 ai pazienti
merged_data_complete <- merged_data_complete %>%
  mutate(
    id = if_else(
      startsWith(id, "MAV"),
      paste0(id, "T0"),
      id
    )
  )

## gruppo
merged_data_complete <- merged_data_complete %>%
  mutate(
    group = if_else(
      startsWith(id, "MAV"),
      "MS patients",
      "Healthy donors"
    )
  )
merged_data_complete <- merged_data_complete %>%
  relocate(group, .after = id)

###scrivi file csv come input per app + SANI
write.csv(merged_data_complete, "C:/Users/gialu/Desktop/tesi_magistrale/Git.hub/Thesis_crepaldi/Input/merged_data_complete_T0.csv", row.names = FALSE)






###########################################
#tenendo anche T24 e T36 per prova prof rolla paper da pubblicare
###############################
###### dataset con anche i sani 
modified_metadata <- read.csv("Input/merged_data.csv")

batch_data <- readxl::read_excel("Input/batch_sequenziamento.xlsx",skip = 1) %>%
  select(ID, sequencing_batch) %>%
  rename(id = ID)
batch_data$id[batch_data$id == "MAV5T15"] <- "MAV05T15"
#rinomino MAV23 con i report giusti
batch_data$id[batch_data$id == "MAV23T12"] <- "MAV23T0"
batch_data$id[batch_data$id == "MAV23T24"] <- "MAV23T12"
batch_data$id[batch_data$id == "MAV23T36"] <- "MAV23T24"

modified_metadata <- modified_metadata %>%
  filter(!ID %in% c("MAV08", "MAV24"))

modified_metadata <- modified_metadata %>%
  pivot_longer(
    cols = matches("^(EDSS|lesioni|relapse|cortisone)_"),
    names_to = c(".value", "tp"),
    names_sep = "_"
  ) %>%
  filter(tp %in% c("T0", "T6", "T12", "T18", "T24", "T36")) %>%
  bind_rows(
    filter(., tp == "T6") %>% mutate(tp = "T3"),
    filter(., tp == "T18") %>% mutate(tp = "T15")
  ) %>%
  filter(!tp %in% c("T6", "T18")) %>%   # rimuove T6 e T18
  mutate(tp = factor(tp, levels = c("T0", "T3", "T12", "T15", "T24", "T36"))) %>%
  arrange(ID, tp) %>%
  group_by(ID) %>%
  mutate(
    EDSS_trend = case_when(
      tp == "T3" ~ case_when(
        EDSS[tp == "T12"] == EDSS[tp == "T0"] ~ "stable",
        EDSS[tp == "T12"] > EDSS[tp == "T0"]  ~ "worsening",
        EDSS[tp == "T12"] < EDSS[tp == "T0"]  ~ "improvement",
        TRUE ~ NA_character_
      ),
      tp == "T15" ~ case_when(
        EDSS[tp == "T24"] == EDSS[tp == "T12"] ~ "stable",
        EDSS[tp == "T24"] > EDSS[tp == "T12"]  ~ "worsening",
        EDSS[tp == "T24"] < EDSS[tp == "T12"]  ~ "improvement",
        TRUE ~ NA_character_
      ),
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  mutate(sample = paste0(ID, tp)) %>%
  rename(original = ID) %>%
  rename(id = sample) %>%
  select(id, everything()) %>%
  select(where(~ !all(is.na(.)))) %>%
  select(-sequencing_batch) %>%
  left_join(batch_data, by = "id") %>%
  filter(!is.na(sequencing_batch))



#mettere in ordine nei grafici app
modified_metadata$id <- replace(
  modified_metadata$id,
  modified_metadata$id %in% c("MAV05T3","MAV07T3","MAV09T3","MAV11T3","MAV13T3","MAV14T3","MAV17T3","MAV18T3","MAV23T3","MAV26T3","MAV27T3","MAV28T3","MAV29T3","MAV30T3","MAV31T3","MAV32T3","MAV34T3"),
  c("MAV05T03","MAV07T03","MAV09T03","MAV11T03","MAV13T03","MAV14T03","MAV17T03","MAV18T03","MAV23T03","MAV26T03","MAV27T03","MAV28T03","MAV29T03","MAV30T03","MAV31T03","MAV32T03","MAV34T03")[match(modified_metadata$id[modified_metadata$id %in% c("MAV05T3","MAV07T3","MAV09T3","MAV11T3","MAV13T3","MAV14T3","MAV17T3","MAV18T3","MAV23T3","MAV26T3","MAV27T3","MAV28T3","MAV29T3","MAV30T3","MAV31T3","MAV32T3","MAV34T3")], c("MAV05T3","MAV07T3","MAV09T3","MAV11T3","MAV13T3","MAV14T3","MAV17T3","MAV18T3","MAV23T3","MAV26T3","MAV27T3","MAV28T3","MAV29T3","MAV30T3","MAV31T3","MAV32T3","MAV34T3"))]
)

cleaned_metadata <- read.csv(
  "Input/metadata_cleaned_completo.csv",
  sep = ";"
)

healthy_metadata <- cleaned_metadata %>%
  filter(grepl("^(HV|AMCTRL|MBCTRL|NBCTRL|RRCTRL)", ID)) %>%
  rename(
    original = ID,
    age_at_t0 = age_at_T0,
    Data_di_nascita = data_nascita
  ) %>%
  mutate(
    id = original,
    tp = factor("T0", levels = levels(modified_metadata$tp)),
    sequencing_batch = sequencing_batch,
    medas_group = NA_character_,
    bmi_group = NA_character_,
    EDSS_trend = NA_character_
  )

## Mantieni solo le colonne presenti in modified_metadata
healthy_metadata <- healthy_metadata %>%
  select(any_of(names(modified_metadata)))

## Aggiungi i sani al dataset finale

# Converte ogni colonna di healthy_metadata nello stesso tipo di modified_metadata
for (col in intersect(names(modified_metadata), names(healthy_metadata))) {
  
  if (is.character(modified_metadata[[col]])) {
    healthy_metadata[[col]] <- as.character(healthy_metadata[[col]])
    
  } else if (is.factor(modified_metadata[[col]])) {
    healthy_metadata[[col]] <- factor(
      healthy_metadata[[col]],
      levels = levels(modified_metadata[[col]])
    )
    
  } else if (is.numeric(modified_metadata[[col]])) {
    healthy_metadata[[col]] <- as.numeric(healthy_metadata[[col]])
    
  } else if (is.integer(modified_metadata[[col]])) {
    healthy_metadata[[col]] <- as.integer(healthy_metadata[[col]])
    
  } else if (is.logical(modified_metadata[[col]])) {
    healthy_metadata[[col]] <- as.logical(healthy_metadata[[col]])
  }
}

modified_metadata <- bind_rows(modified_metadata, healthy_metadata)
modified_metadata <- modified_metadata %>%
  mutate(group = if_else( startsWith(id, "MAV"),"Affected MAV", "Healthy"))

write.csv(modified_metadata, "C:/Users/gialu/Desktop/tesi_magistrale/Git.hub/Thesis_crepaldi/Input/merged_data_complete_T0.csv", row.names = FALSE)


