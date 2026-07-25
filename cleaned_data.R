source("Libraries.R")
cleaned_metadata <- read.csv("Input/metadata_cleaned_completo.csv", sep=";")
merged_data <- read.csv("Input/merged_data.csv")
merged <- read.csv2("C:\\Users\\gialu\\Desktop\\merged.csv")


cleaned_metadata= cleaned_metadata[startsWith(cleaned_metadata$ID, "MAV"), ]

cleaned_metadata$ID <- sub("T0$", "", cleaned_metadata$ID)

confronto <- CECI %>%
  select(ID, cambio_terapia) %>%
  rename(ceci = cambio_terapia) %>%
  left_join(
    merged_data %>%
      select(ID, cambio_terapia) %>%
      rename(merged_data = cambio_terapia),by = "ID"
  )


