source("Libraries.R")
cleaned_metadata <- read.csv("Input/metadata_cleaned_completo.csv", sep=";")
merged_data <- read.csv("Input/merged_data.csv")
merged <- read.csv2("C:\\Users\\gialu\\Desktop\\merged.csv")

table(merged$GROUP, merged$gadolinium_contrast)

#####
merged= merged[startsWith(merged$ID, "MAV"), ]

cleaned_metadata$ID <- sub("T0$", "", cleaned_metadata$ID)
dati_radiologici = read_excel("Input/merged_dati_radiologici.xlsx")


confronto <- merged %>%
  select(ID, lesioni_T12) %>%
  rename(ceci = lesioni_T12) %>%
  left_join(
    merged_data %>%
      select(ID, lesioni_T12) %>%
      rename(merged_data = lesioni_T12),by = "ID"
  )

str(merged_data$gc_treatment)

########


########################################
#creazione excel per ale
#########################################
merged_subset <- merged %>%
  select(
    ID, CENTRO, data_prelievo_T0, bmi, age_at_diagnosis, MEDAS,
    smoking_habit, physical_activity, gc_treatment, lesion_burden,
    bone_marrow_lesions, subtentorial_lesions, gadolinium_contrast
  )


library("writexl")
write_xlsx(merged_subset, "C:\\Users\\gialu\\Desktop\\merged_subset.xlsx")




