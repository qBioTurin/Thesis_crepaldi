source("Libraries.R")

merged_data <- read.csv("Input/merged_data.csv")

batch_data <- readxl::read_excel("Input/batch_sequenziamento.xlsx",skip = 1) %>%
  select(ID, sequencing_batch) %>%
  rename(id = ID)

filtered_patients <- merged_data %>%
  filter(ID %in% c("MAV05", "MAV07", "MAV09", "MAV11", "MAV13", "MAV14", "MAV34")) %>%
  distinct(ID, .keep_all = TRUE)

filtered_patients <- filtered_patients %>%
  pivot_longer(
    cols = matches("^(EDSS|lesioni|relapse|cortisone)_"),
    names_to = c(".value", "tp"),
    names_sep = "_"
  ) %>%
  filter(tp %in% c("T0", "T12")) %>%
  bind_rows(filter(., tp == "T0") %>% mutate(tp = "T3")) %>%
  mutate(tp = factor(tp, levels = c("T0", "T3", "T12"))) %>%
  arrange(ID, tp) %>%
  group_by(ID) %>%
  mutate(
    EDSS_trend = case_when(
      tp != "T3" ~ NA_character_,
      EDSS[tp == "T12"] == EDSS[tp == "T0"] ~ "stable",
      EDSS[tp == "T12"] > EDSS[tp == "T0"]  ~ "worsening",
      EDSS[tp == "T12"] < EDSS[tp == "T0"]  ~ "improvement",
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
  left_join(batch_data, by = "id")



write.csv(filtered_patients, "C:/Users/gialu/Desktop/tesi_magistrale/Git.hub/Thesis_crepaldi/Input/filtered_patients.csv", row.names = FALSE)