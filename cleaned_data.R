source("Libraries.R")
cleaned_metadata <- read.csv("Input/metadata_cleaned_completo.csv", sep=";")
merged_data <- read.csv("Input/merged_data.csv")
CECI <- read.csv2("C:\\Users\\gialu\\Desktop\\merged.csv")


cleaned_metadata= cleaned_metadata[startsWith(cleaned_metadata$ID, "MAV"), ]

cleaned_metadata$ID <- sub("T0$", "", cleaned_metadata$ID)

confronto <- cleaned_metadata %>%
  select(ID, gadolinium_contrast) %>%
  rename(cleaned = gadolinium_contrast) %>%
  left_join(
    merged_data %>%
      select(ID, gadolinium_contrast) %>%
      rename(merged = gadolinium_contrast),by = "ID"
  )

