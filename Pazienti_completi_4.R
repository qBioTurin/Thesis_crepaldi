source("Libraries.R")

modified_metadata <- read.csv("Input/merged_data.csv")

batch_data <- readxl::read_excel("Input/batch_sequenziamento.xlsx",skip = 1) %>%
  select(ID, sequencing_batch) %>%
  rename(id = ID)

modified_metadata <- modified_metadata %>%
  pivot_longer(
    cols = matches("^(EDSS|lesioni|relapse|cortisone)_"),
    names_to = c(".value", "tp"),
    names_sep = "_"
  ) %>%
  filter(tp %in% c("T0", "T12", "T24", "T36")) %>%
  bind_rows(
    filter(., tp == "T0") %>% mutate(tp = "T3"),
    filter(., tp == "T12") %>% mutate(tp = "T15")
  ) %>%
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
  left_join(batch_data, by = "id")


write.csv(modified_metadata, "C:/Users/gialu/Desktop/tesi_magistrale/Git.hub/Thesis_crepaldi/Input/modified_metadata.csv", row.names = FALSE)


