source("Libraries.R")
merged_data <- read.csv("Input/merged_data.csv")
data <- merged_data

data <- data %>%
  rowwise() %>%
  mutate(
    EDSS_finale = last(na.omit(c(
      EDSS_T0,
      EDSS_T6,
      EDSS_T12,
      EDSS_T18,
      EDSS_T24,
      EDSS_T30,
      EDSS_T36
    )))
  ) %>%
  ungroup()

data <- data %>%
  mutate(
    delta_EDSS = EDSS_finale - EDSS_T0
  )

data <- data %>%
  mutate(
    andamento_EDSS = case_when(
      
      delta_EDSS >= 0.5 ~ "Peggioramento",
      
      delta_EDSS <= -0.5 ~ "Miglioramento",
      
      abs(delta_EDSS) < 0.5 & EDSS_T0 < 4 ~
        "Stabile basso",
      
      abs(delta_EDSS) < 0.5 & EDSS_T0 >= 4 ~
        "Stabile alto"
    )
  )

str(data)
