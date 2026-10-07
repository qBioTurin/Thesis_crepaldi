source("Libraries.R")
merged_data <- read.csv("Input/merged_data.csv")

merged_data <- merged_data %>%
  filter(!ID %in% c("MAV08", "MAV12", "MAV16", "MAV22",
                    "MAV24", "MAV25", "MAV31"))
summary(merged_data$age_at_t0)

merged_data <- merged_data %>%
  mutate(
    age_group = if_else(age_at_t0 >= 38, "old", "young"),
    bmi_group = if_else(bmi >= 24, "High", "Low")
  )


table(merged_data$sex)                    
table(merged_data$bmi_group)
table(merged_data$age_group)
table(merged_data$lesion_burden)
table(merged_data$gadolinium_contrast)
table(merged_data$bone_marrow_lesions)
table(merged_data$subtentorial_lesions)
