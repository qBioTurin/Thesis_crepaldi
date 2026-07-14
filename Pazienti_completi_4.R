source("Libraries.R")

modified_metadata <- read.csv("Input/merged_data.csv")

batch_data <- readxl::read_excel("Input/batch_sequenziamento.xlsx",skip = 1) %>%
  select(ID, sequencing_batch) %>%
  rename(id = ID)
batch_data$id[batch_data$id == "MAV5T15"] <- "MAV05T15"

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

modified_metadata <- modified_metadata %>%
  mutate(
    medas_group = if_else(MEDAS >= 46, "High_adherence", "Low_adherence"),
    bmi_group = if_else(bmi >= 24, "High_bmi", "Low_bmi")
  )


write.csv(modified_metadata, "C:/Users/gialu/Desktop/tesi_magistrale/Git.hub/Thesis_crepaldi/Input/modified_metadata.csv", row.names = FALSE)

