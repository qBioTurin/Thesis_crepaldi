source("Libraries.R")
cleaned_metadata <- read.csv("Input/metadata_cleaned_completo.csv", sep=";")
merged_data <- read.csv("Input/merged_data.csv")

summary(cleaned_metadata$bmi)

##gc_treatment
#########################################################################
merged <- read.csv2("C:/Users/gialu/Desktop/merged.csv", stringsAsFactors = FALSE)
merged[merged$ID == "MAV32", c("bmi", "smoking_activity", "physical_activity")] <- list(33.20, "yes", "sedentary")
merged[merged$ID == "TEC25", c("bmi", "smoking_activity", "physical_activity")] <- list(21.88, "no", "sedentary")
merged[merged$ID == "MAV31", c("bmi", "smoking_activity", "physical_activity")] <- list(27.43, "no", "sedentary")
merged[merged$ID == "TEC26", c("bmi", "smoking_activity", "physical_activity")] <- list(22.46, "yes", "active")
merged[merged$ID == "TEC27", c("bmi", "smoking_activity", "physical_activity")] <- list(25.17, "no", "partially_active")
merged[merged$ID == "TEC32", c( "physical_activity")] <- list(NA)

merged$bmi

#dataset mio x fare surv con MAV
merged_data <- read.csv("Input/merged_data.csv")

## COX GC TREATMENT SOLO MAV 

merged_data <- merged_data %>%
  mutate(
    gadolinium_contrast = case_when(
      gadolinium_contrast %in% c("Yes", "1") ~ 1,
      gadolinium_contrast %in% c("No", "0")  ~ 0,
      TRUE ~ NA_real_),
    
    lesion_burden = case_when(
      lesion_burden %in% c("Yes", "1") ~ 1,
      lesion_burden %in% c("No", "0")  ~ 0,
      TRUE ~ NA_real_),
    
    bone_marrow_lesions = case_when(
      bone_marrow_lesions %in% c("Yes", "1") ~ 1,
      bone_marrow_lesions %in% c("No", "0")  ~ 0,
      TRUE ~ NA_real_),
    
    subtentorial_lesions = case_when(
      subtentorial_lesions %in% c("Yes", "1") ~ 1,
      subtentorial_lesions %in% c("No", "0")  ~ 0,
      TRUE ~ NA_real_),
    
    gadolinium_contrast = factor(gadolinium_contrast, levels = c(0, 1)),
    lesion_burden = factor(lesion_burden, levels = c(0, 1)),
    bone_marrow_lesions = factor(bone_marrow_lesions, levels = c(0, 1)),
    subtentorial_lesions = factor(subtentorial_lesions, levels = c(0, 1)),
    gc_treatment = factor(gc_treatment))

merged_data <- merged_data %>%
  rowwise() %>%
  mutate(
    time_followup = case_when(
      !is.na(EDSS_T36) ~ 36,
      !is.na(EDSS_T30) ~ 30,
      !is.na(EDSS_T24) ~ 24,
      !is.na(EDSS_T18) ~ 18,
      !is.na(EDSS_T12) ~ 12,
      !is.na(EDSS_T6)  ~ 6,
      TRUE ~ 0)
  ) %>%
  ungroup()

# NEDA-3
merged_data <- merged_data %>%
  rowwise() %>%
  mutate(
    EDSS_finale = last(na.omit(c(
      EDSS_T0, EDSS_T6, EDSS_T12, EDSS_T18,
      EDSS_T24, EDSS_T30, EDSS_T36)))
  ) %>%
  ungroup() %>%
  mutate(delta_EDSS = as.numeric(as.character(EDSS_finale)) -
      as.numeric(as.character(EDSS_T0)),  EDSS_progression = case_when(
      EDSS_T0 == 0                ~ delta_EDSS >= 1.5,
      EDSS_T0 >= 1 & EDSS_T0 <= 5 ~ delta_EDSS >= 1,
      EDSS_T0 > 5                 ~ delta_EDSS >= 0.5,
      TRUE                        ~ NA
    ),
    clinical_activity =
      relapse_T6  == "Yes" |
      relapse_T12 == "Yes" |
      relapse_T18 == "Yes" |
      relapse_T24 == "Yes" |
      relapse_T30 == "Yes" |
      relapse_T36 == "Yes",
    
    radiological_activity =
      lesioni_T6  == "Yes" |
      lesioni_T12 == "Yes" |
      lesioni_T18 == "Yes" |
      lesioni_T24 == "Yes" |
      lesioni_T30 == "Yes" |
      lesioni_T36 == "Yes",
    
    status_NEDA3 = if_else(
      clinical_activity | radiological_activity | EDSS_progression,
      1, 0
    )
  )

# dataframe per cox — SOLO VARIABILI PRESCELTE PRIMA
df_cox <- merged_data %>%
  select(
    time_followup,
    status_NEDA3,
    lesion_burden,
    gadolinium_contrast,
    bone_marrow_lesions,
    subtentorial_lesions,
    gc_treatment) %>%
  na.omit()


# modello cox ridotto 
cox_neda3 <- coxph(
  Surv(time_followup, status_NEDA3) ~ 
    gc_treatment + lesion_burden + gadolinium_contrast +
    bone_marrow_lesions + subtentorial_lesions,
  data = df_cox)

summary(cox_neda3)

# forest 
df_forest <- tidy(cox_neda3, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    HR = exp(estimate),
    lower = exp(conf.low),
    upper = exp(conf.high),
    term = factor(term, levels = rev(term))
  )

ggplot(df_forest, aes(x = term, y = HR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = "Forest plot – Cox NEDA-3",
    x = "",
    y = "Hazard Ratio (HR)"
  )



### GC TREATMENT SOLO mav NAIVE
merged_data_naive <- merged_data %>% 
  filter(naive == "Yes")

merged_data_naive <- merged_data_naive %>%
  rowwise() %>%
  mutate(
    time_followup = case_when(
      !is.na(EDSS_T36) ~ 36,
      !is.na(EDSS_T30) ~ 30,
      !is.na(EDSS_T24) ~ 24,
      !is.na(EDSS_T18) ~ 18,
      !is.na(EDSS_T12) ~ 12,
      !is.na(EDSS_T6)  ~ 6,
      TRUE ~ 0)) %>%
  ungroup()


# NEDA-3
merged_data_naive <- merged_data_naive %>%
  rowwise() %>%
  mutate(
    EDSS_finale = last(na.omit(c(
      EDSS_T0, EDSS_T6, EDSS_T12, EDSS_T18,
      EDSS_T24, EDSS_T30, EDSS_T36)))) %>%
  ungroup() %>%
  mutate(
    delta_EDSS = as.numeric(as.character(EDSS_finale)) -
      as.numeric(as.character(EDSS_T0)),
    EDSS_progression = case_when(EDSS_T0 == 0  ~ delta_EDSS >= 1.5,
      EDSS_T0 >= 1 & EDSS_T0 <= 5 ~ delta_EDSS >= 1,
      EDSS_T0 > 5                 ~ delta_EDSS >= 0.5,
      TRUE                        ~ NA
    ),
    clinical_activity =
      relapse_T6  == "Yes" |
      relapse_T12 == "Yes" |
      relapse_T18 == "Yes" |
      relapse_T24 == "Yes" |
      relapse_T30 == "Yes" |
      relapse_T36 == "Yes",
    
    radiological_activity =
      lesioni_T6  == "Yes" |
      lesioni_T12 == "Yes" |
      lesioni_T18 == "Yes" |
      lesioni_T24 == "Yes" |
      lesioni_T30 == "Yes" |
      lesioni_T36 == "Yes",
    
    status_NEDA3 = if_else(
      clinical_activity | radiological_activity | EDSS_progression,
      1, 0
    )
  )

# dataframe per cox — SOLO VARIABILI PRESCELTE PRIMA
df_cox <- merged_data_naive %>%
  select(
    time_followup,
    status_NEDA3,
    lesion_burden,
    gadolinium_contrast,
    bone_marrow_lesions,
    subtentorial_lesions,
    gc_treatment) %>%
  na.omit()


# modello cox ridotto
cox_neda3 <- coxph(
  Surv(time_followup, status_NEDA3) ~ 
    gc_treatment + lesion_burden + gadolinium_contrast +
    bone_marrow_lesions + subtentorial_lesions,
  data = df_cox)

summary(cox_neda3)


# forest 
df_forest <- tidy(cox_neda3, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    HR = exp(estimate),
    lower = exp(conf.low),
    upper = exp(conf.high),
    term = factor(term, levels = rev(term)))

ggplot(df_forest, aes(x = term, y = HR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = "Forest plot – Cox NEDA-3",
    x = "",
    y = "Hazard Ratio (HR)"
  )




#################gc_at_least#######################################
# GC AT LEAST SOLO mav
# creo colonna gc at least con livelli corretti 

merged_data <- merged_data %>%
  mutate(cortisone_any = rowSums(select(., starts_with("cortisone_T")) != 0, na.rm = TRUE) > 0,
      gc_at_least = case_when(
      gc_treatment == "Negative" & !cortisone_any ~ 0,
      gc_treatment == "Positive" & !cortisone_any ~ 1,
      gc_treatment == "Negative" & cortisone_any ~ 2,
      gc_treatment == "Positive" & cortisone_any ~ 3,
      TRUE ~ NA_real_),
        gc_at_least = factor(gc_at_least, levels = c(0,1,2,3))) %>%
  select(-cortisone_any)

class(merged_data$gc_at_least) #deve essere factor

merged_data <- merged_data %>%
  mutate(
    gadolinium_contrast = case_when(
      gadolinium_contrast %in% c("Yes", "1") ~ 1,
      gadolinium_contrast %in% c("No", "0")  ~ 0,
      TRUE ~ NA_real_),
    
    lesion_burden = case_when(
      lesion_burden %in% c("Yes", "1") ~ 1,
      lesion_burden %in% c("No", "0")  ~ 0,
      TRUE ~ NA_real_),
    
    bone_marrow_lesions = case_when(
      bone_marrow_lesions %in% c("Yes", "1") ~ 1,
      bone_marrow_lesions %in% c("No", "0")  ~ 0,
      TRUE ~ NA_real_),
    
    subtentorial_lesions = case_when(
      subtentorial_lesions %in% c("Yes", "1") ~ 1,
      subtentorial_lesions %in% c("No", "0")  ~ 0,
      TRUE ~ NA_real_),
    
    gadolinium_contrast = factor(gadolinium_contrast, levels = c(0, 1)),
    lesion_burden = factor(lesion_burden, levels = c(0, 1)),
    bone_marrow_lesions = factor(bone_marrow_lesions, levels = c(0, 1)),
    subtentorial_lesions = factor(subtentorial_lesions, levels = c(0, 1)),
    gc_at_least = factor(gc_at_least))

merged_data <- merged_data %>%
  rowwise() %>%
  mutate(
    time_followup = case_when(
      !is.na(EDSS_T36) ~ 36,
      !is.na(EDSS_T30) ~ 30,
      !is.na(EDSS_T24) ~ 24,
      !is.na(EDSS_T18) ~ 18,
      !is.na(EDSS_T12) ~ 12,
      !is.na(EDSS_T6)  ~ 6,
      TRUE ~ 0)) %>%
  ungroup()


# NEDA-3
merged_data <- merged_data %>%
  rowwise() %>%
  mutate(
    EDSS_finale = last(na.omit(c(
      EDSS_T0, EDSS_T6, EDSS_T12, EDSS_T18,
      EDSS_T24, EDSS_T30, EDSS_T36)))) %>%
  ungroup() %>%
  mutate(EDSS_T0_num = as.numeric(as.character(EDSS_T0)),
      delta_EDSS = as.numeric(as.character(EDSS_finale)) - EDSS_T0_num,
      EDSS_progression = case_when(
      EDSS_T0_num == 0 ~ delta_EDSS >= 1.5,
      EDSS_T0_num >= 1 & EDSS_T0_num <= 5 ~ delta_EDSS >= 1,
      EDSS_T0_num > 5 ~ delta_EDSS >= 0.5,
      TRUE ~ NA),
      clinical_activity =
      relapse_T6  == "Yes" |
      relapse_T12 == "Yes" |
      relapse_T18 == "Yes" |
      relapse_T24 == "Yes" |
      relapse_T30 == "Yes" |
      relapse_T36 == "Yes",
    
    radiological_activity =
      lesioni_T6  == "Yes" |
      lesioni_T12 == "Yes" |
      lesioni_T18 == "Yes" |
      lesioni_T24 == "Yes" |
      lesioni_T30 == "Yes" |
      lesioni_T36 == "Yes",
    status_NEDA3 = if_else(
      clinical_activity | radiological_activity | EDSS_progression,
      1, 0
    )
  )

# dataframe per cox — SOLO VARIABILI PRESCELTE PRIMA
df_cox <- merged_data %>%
  select(
    time_followup,
    status_NEDA3,
    lesion_burden,
    gadolinium_contrast,
    bone_marrow_lesions,
    subtentorial_lesions,
    gc_at_least) %>%
  na.omit()


# modello cox ridotto
cox_neda3 <- coxph(
  Surv(time_followup, status_NEDA3) ~ 
    gc_at_least + lesion_burden + gadolinium_contrast +
    bone_marrow_lesions + subtentorial_lesions,
  data = df_cox)

summary(cox_neda3)

# forest 
df_forest <- tidy(cox_neda3, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    HR = exp(estimate),
    lower = exp(conf.low),
    upper = exp(conf.high),
    term = factor(term, levels = rev(term)))

ggplot(df_forest, aes(x = term, y = HR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = "Forest plot – Cox NEDA-3",
    x = "",
    y = "Hazard Ratio (HR)"
  )



#### GC AT LEAST PER mav NAIVE
merged_data <- merged_data %>%
  mutate(
    cortisone_any = rowSums(select(., starts_with("cortisone_T")) != 0, na.rm = TRUE) > 0,
      gc_at_least = case_when(
      gc_treatment == "Negative" & !cortisone_any ~ 0,
      gc_treatment == "Positive" & !cortisone_any ~ 1,
      gc_treatment == "Negative" & cortisone_any ~ 2,
      gc_treatment == "Positive" & cortisone_any ~ 3,
      TRUE ~ NA_real_),
      gc_at_least = factor(gc_at_least, levels = c(0,1,2,3))) %>%
  select(-cortisone_any)

class(merged_data$gc_at_least) #deve essere factor

merged_data <- merged_data %>%
  mutate(
    gadolinium_contrast = case_when(
      gadolinium_contrast %in% c("Yes", "1") ~ 1,
      gadolinium_contrast %in% c("No", "0")  ~ 0,
      TRUE ~ NA_real_),
    lesion_burden = case_when(
      lesion_burden %in% c("Yes", "1") ~ 1,
      lesion_burden %in% c("No", "0")  ~ 0,
      TRUE ~ NA_real_),
    bone_marrow_lesions = case_when(
      bone_marrow_lesions %in% c("Yes", "1") ~ 1,
      bone_marrow_lesions %in% c("No", "0")  ~ 0,
      TRUE ~ NA_real_),
    subtentorial_lesions = case_when(
      subtentorial_lesions %in% c("Yes", "1") ~ 1,
      subtentorial_lesions %in% c("No", "0")  ~ 0,
      TRUE ~ NA_real_),
    gadolinium_contrast = factor(gadolinium_contrast, levels = c(0, 1)),
    lesion_burden = factor(lesion_burden, levels = c(0, 1)),
    bone_marrow_lesions = factor(bone_marrow_lesions, levels = c(0, 1)),
    subtentorial_lesions = factor(subtentorial_lesions, levels = c(0, 1)),
    gc_at_least = factor(gc_at_least))

merged_data <- merged_data %>%
  rowwise() %>%
  mutate(
    time_followup = case_when(
      !is.na(EDSS_T36) ~ 36,
      !is.na(EDSS_T30) ~ 30,
      !is.na(EDSS_T24) ~ 24,
      !is.na(EDSS_T18) ~ 18,
      !is.na(EDSS_T12) ~ 12,
      !is.na(EDSS_T6)  ~ 6,
      TRUE ~ 0)) %>%
  ungroup()

# NEDA-3
merged_data <- merged_data %>%
  rowwise() %>%
  mutate(EDSS_finale = last(na.omit(c(
      EDSS_T0, EDSS_T6, EDSS_T12, EDSS_T18,
      EDSS_T24, EDSS_T30, EDSS_T36)))) %>%
  ungroup() %>%
  mutate(EDSS_T0_num = as.numeric(as.character(EDSS_T0)),
    delta_EDSS = as.numeric(as.character(EDSS_finale)) - EDSS_T0_num,
    EDSS_progression = case_when(
      EDSS_T0_num == 0 ~ delta_EDSS >= 1.5,
      EDSS_T0_num >= 1 & EDSS_T0_num <= 5 ~ delta_EDSS >= 1,
      EDSS_T0_num > 5 ~ delta_EDSS >= 0.5,
      TRUE ~ NA),
    clinical_activity =
      relapse_T6  == "Yes" |
      relapse_T12 == "Yes" |
      relapse_T18 == "Yes" |
      relapse_T24 == "Yes" |
      relapse_T30 == "Yes" |
      relapse_T36 == "Yes",
    radiological_activity =
      lesioni_T6  == "Yes" |
      lesioni_T12 == "Yes" |
      lesioni_T18 == "Yes" |
      lesioni_T24 == "Yes" |
      lesioni_T30 == "Yes" |
      lesioni_T36 == "Yes",
    status_NEDA3 = if_else(
      clinical_activity | radiological_activity | EDSS_progression,
      1, 0
    )
  )


merged_data_naive <- merged_data %>% 
  filter(naive == "Yes")

# dataframe per cox — SOLO VARIABILI PRESCELTE
df_cox <- merged_data_naive %>%
  select(
    time_followup,
    status_NEDA3,
    lesion_burden,
    gadolinium_contrast,
    bone_marrow_lesions,
    subtentorial_lesions,
    gc_at_least) %>%
  na.omit()


# modello cox ridotto
cox_neda3 <- coxph(
  Surv(time_followup, status_NEDA3) ~ 
    gc_at_least + lesion_burden + gadolinium_contrast +
    bone_marrow_lesions + subtentorial_lesions,
  data = df_cox)

summary(cox_neda3)

# forest 
df_forest <- tidy(cox_neda3, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    HR = exp(estimate),
    lower = exp(conf.low),
    upper = exp(conf.high),
    term = factor(term, levels = rev(term)))

ggplot(df_forest, aes(x = term, y = HR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = "Forest plot – Cox NEDA-3",
    x = "",
    y = "Hazard Ratio (HR)"
  )












###########################
#####prova
############
merged_data <- merged_data %>%
  mutate(
    cortisone_any = rowSums(select(., starts_with("cortisone_T")) != 0, na.rm = TRUE) > 0,
    
    gc_at_least = case_when(
      gc_treatment == "Negative" & !cortisone_any ~ 0,
      gc_treatment == "Positive" & !cortisone_any ~ 1,
      gc_treatment == "Negative" & cortisone_any ~ 2,
      gc_treatment == "Positive" & cortisone_any ~ 3,
      TRUE ~ NA_real_
    ),
    
    gc_at_least = factor(gc_at_least, levels = c(0,1,2,3))
  ) %>%
  select(-cortisone_any)

class(merged_data$gc_at_least) #deve essere factor

merged_data <- merged_data %>%
  mutate(
    gadolinium_contrast = case_when(
      gadolinium_contrast %in% c("Yes", "1") ~ 1,
      gadolinium_contrast %in% c("No", "0")  ~ 0,
      TRUE ~ NA_real_
    ),
    
    lesion_burden = case_when(
      lesion_burden %in% c("Yes", "1") ~ 1,
      lesion_burden %in% c("No", "0")  ~ 0,
      TRUE ~ NA_real_
    ),
    
    bone_marrow_lesions = case_when(
      bone_marrow_lesions %in% c("Yes", "1") ~ 1,
      bone_marrow_lesions %in% c("No", "0")  ~ 0,
      TRUE ~ NA_real_
    ),
    
    subtentorial_lesions = case_when(
      subtentorial_lesions %in% c("Yes", "1") ~ 1,
      subtentorial_lesions %in% c("No", "0")  ~ 0,
      TRUE ~ NA_real_
    ),
    
    gadolinium_contrast = factor(gadolinium_contrast, levels = c(0, 1)),
    lesion_burden = factor(lesion_burden, levels = c(0, 1)),
    bone_marrow_lesions = factor(bone_marrow_lesions, levels = c(0, 1)),
    subtentorial_lesions = factor(subtentorial_lesions, levels = c(0, 1)),
    gc_at_least = factor(gc_at_least),
    age_at_t0 = as.numeric(age_at_t0),
    sex = factor(sex)
  )

merged_data <- merged_data %>%
  rowwise() %>%
  mutate(
    time_followup = case_when(
      !is.na(EDSS_T36) ~ 36,
      !is.na(EDSS_T30) ~ 30,
      !is.na(EDSS_T24) ~ 24,
      !is.na(EDSS_T18) ~ 18,
      !is.na(EDSS_T12) ~ 12,
      !is.na(EDSS_T6)  ~ 6,
      TRUE ~ 0
    )
  ) %>%
  ungroup()


# NEDA-3
merged_data <- merged_data %>%
  rowwise() %>%
  mutate(
    EDSS_finale = last(na.omit(c(
      EDSS_T0, EDSS_T6, EDSS_T12, EDSS_T18,
      EDSS_T24, EDSS_T30, EDSS_T36
    )))
  ) %>%
  ungroup() %>%
  mutate(
    EDSS_T0_num = as.numeric(as.character(EDSS_T0)),
    
    delta_EDSS = as.numeric(as.character(EDSS_finale)) - EDSS_T0_num,
    
    EDSS_progression = case_when(
      EDSS_T0_num == 0 ~ delta_EDSS >= 1.5,
      EDSS_T0_num >= 1 & EDSS_T0_num <= 5 ~ delta_EDSS >= 1,
      EDSS_T0_num > 5 ~ delta_EDSS >= 0.5,
      TRUE ~ NA
    ),
    
    clinical_activity =
      relapse_T6  == "Yes" |
      relapse_T12 == "Yes" |
      relapse_T18 == "Yes" |
      relapse_T24 == "Yes" |
      relapse_T30 == "Yes" |
      relapse_T36 == "Yes",
    
    radiological_activity =
      lesioni_T6  == "Yes" |
      lesioni_T12 == "Yes" |
      lesioni_T18 == "Yes" |
      lesioni_T24 == "Yes" |
      lesioni_T30 == "Yes" |
      lesioni_T36 == "Yes",
    
    status_NEDA3 = if_else(
      clinical_activity | radiological_activity | EDSS_progression,
      1, 0
    )
  )

# dataframe per cox
df_cox <- merged_data %>%
  select(
    time_followup,
    status_NEDA3,
    lesion_burden,
    gadolinium_contrast,
    bone_marrow_lesions,
    subtentorial_lesions,
    gc_at_least,
    sex,
    age_at_t0
  ) %>%
  na.omit()


# modello cox ridotto 
cox_neda3 <- coxph(
  Surv(time_followup, status_NEDA3) ~ 
    gc_at_least + lesion_burden + gadolinium_contrast +
    bone_marrow_lesions + subtentorial_lesions + sex + age_at_t0,
  data = df_cox
)

summary(cox_neda3)


# forest 

df_forest <- tidy(cox_neda3, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    HR = exp(estimate),
    lower = exp(conf.low),
    upper = exp(conf.high),
    term = factor(term, levels = rev(term))
  )

ggplot(df_forest, aes(x = term, y = HR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = "Forest plot – Cox NEDA-3",
    x = "",
    y = "Hazard Ratio (HR)"
  )

summary(cox_neda3)
table(df_cox$status_NEDA3)
table(df_cox$gc_at_least, df_cox$status_NEDA3)







###########KAPLAN MAIER DOPO COX
#subtentorial x tesi #tutti mav stratificati
##########
merged_data <- merged_data %>%
  mutate(subtento = ifelse( subtentorial_lesions == "Yes", 1, 0),
    sub_group = case_when(
      subtento == 1 ~ "MAV subtentorial +",
      subtento == 0 ~ "MAV subtentorial -"
    )
  )

merged_data$therapy_group <- factor(
  merged_data$sub_group,
  levels = c("MAV subtentorial +", "MAV subtentorial -"))

# Ultimo EDSS disponibile
merged_data <- merged_data %>%
  rowwise() %>%
  mutate(EDSS_last = {
      x <- c(EDSS_T36, EDSS_T30, EDSS_T24,
             EDSS_T18, EDSS_T12, EDSS_T6, EDSS_T0)
      x[which(!is.na(x))[1]]
    }
  ) %>%
  ungroup()


# Peggioramento EDSS
merged_data <- merged_data %>%
  mutate(EDSS_worsening = case_when(
      EDSS_T0 == 0 ~ (EDSS_last - EDSS_T0) >= 1.5,
      EDSS_T0 >= 1 & EDSS_T0 <= 5 ~ (EDSS_last - EDSS_T0) >= 1,
      EDSS_T0 > 5 ~ (EDSS_last - EDSS_T0) >= 0.5,
      TRUE ~ FALSE
    )
  )

# Eventi NEDA
merged_data <- merged_data %>%
  mutate(
    event_6_C  = EDSS_worsening | relapse_T6  == "Yes" | lesioni_T6  == "Yes",
    event_12_C = EDSS_worsening | relapse_T12 == "Yes" | lesioni_T12 == "Yes",
    event_18_C = EDSS_worsening | relapse_T18 == "Yes" | lesioni_T18 == "Yes",
    event_24_C = EDSS_worsening | relapse_T24 == "Yes" | lesioni_T24 == "Yes",
    event_30_C = EDSS_worsening | relapse_T30 == "Yes" | lesioni_T30 == "Yes",
    event_36_C = EDSS_worsening | relapse_T36 == "Yes" | lesioni_T36 == "Yes")


# Tempo al primo evento
merged_data <- merged_data %>%
  mutate(
    event_time_C = case_when(
      event_6_C  ~ 6,
      event_12_C ~ 12,
      event_18_C ~ 18,
      event_24_C ~ 24,
      event_30_C ~ 30,
      event_36_C ~ 36,
      TRUE ~ NA_real_),
    event_status_C = ifelse(is.na(event_time_C), 0, 1),
    time_C = ifelse(is.na(event_time_C), 36, event_time_C))


# Kaplan-Meier
surv_obj <- Surv(time = merged_data$time_C, event = merged_data$event_status_C)

fit <- survfit(surv_obj ~ sub_group, data = merged_data)

ggsurvplot(
  fit,
  data = merged_data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  title= "NEDA-3 survival curve for MAV patients stratified by subtentorial lesions",
  legend.title = "Groups",
  legend.labs= c("MAV subtentorial lesions -", "MAV subtentorial lesions +"),
  xlab = "Months",
  ylab = "NEDA-3 loss probability"
)



##################################################
#NEDA 3 SURVIVAL GENERALE MAV E TEC
#######################################
#GRAFICO DA INSERIRE IN TESI!!!!!!!!!!!!!!!!!!!!!!!!!
###################################################
merged$GROUP <- factor(
  merged$GROUP,
  levels = c("MAV", "TEC"))


# Ultimo EDSS disponibile
merged <- merged %>%
  rowwise() %>%
  mutate(EDSS_last = {
      x <- c(EDSS_T36, EDSS_T30, EDSS_T24,
             EDSS_T18, EDSS_T12, EDSS_T6, EDSS_T0)
      x[which(!is.na(x))[1]]
    }
  ) %>%
  ungroup()


# Peggioramento EDSS
merged <- merged %>%
  mutate(EDSS_worsening = case_when(
      EDSS_T0 == 0 ~ (EDSS_last - EDSS_T0) >= 1.5,
      EDSS_T0 >= 1 & EDSS_T0 <= 5 ~ (EDSS_last - EDSS_T0) >= 1,
      EDSS_T0 > 5 ~ (EDSS_last - EDSS_T0) >= 0.5,
      TRUE ~ FALSE
    )
  )


# Eventi NEDA
merged <- merged %>%
  mutate(
    event_6_C  = EDSS_worsening | relapse_T6  == 1 | lesioni_T6  == 1,
    event_12_C = EDSS_worsening | relapse_T12 == 1 | lesioni_T12 == 1,
    event_18_C = EDSS_worsening | relapse_T18 == 1 | lesioni_T18 == 1,
    event_24_C = EDSS_worsening | relapse_T24 == 1 | lesioni_T24 == 1,
    event_30_C = EDSS_worsening | relapse_T30 == 1 | lesioni_T30 == 1,
    event_36_C = EDSS_worsening | relapse_T36 == 1 | lesioni_T36 == 1)


# Tempo al primo evento
merged <- merged %>%
  mutate(
    event_time_C = case_when(
      event_6_C  ~ 6,
      event_12_C ~ 12,
      event_18_C ~ 18,
      event_24_C ~ 24,
      event_30_C ~ 30,
      event_36_C ~ 36,
      TRUE ~ NA_real_),
    event_status_C = ifelse(is.na(event_time_C), 0, 1),
    time_C = ifelse(is.na(event_time_C), 36, event_time_C))

surv_obj4 <- Surv(time = merged$time_C, event = merged$event_status_C)
fit4 <- survfit(surv_obj4 ~ GROUP, data = merged)

ggsurvplot(
  fit4,
  data = merged,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  title = "NEDA-3 free survival curves",
  legend.title = "Groups",
  legend.labs = c("MAV", "TEC"),
  xlab = "Months",
  ylab = "NEDA-3 loss probability"
)



#################SOLO NAIVE !!!
#NEDA 3 SURVIVAL
#GRAFICO DA INSERIRE IN TESI! (ALLA FINE NO)
###################################################
merged$GROUP <- factor(
  merged$GROUP,
  levels = c("MAV", "TEC"))


# Ultimo EDSS disponibile

merged_naive <- merged %>%
  filter(naive != "no")

merged_naive <- merged_naive %>%
  rowwise() %>%
  mutate(EDSS_last = {
      x <- c(EDSS_T36, EDSS_T30, EDSS_T24,
             EDSS_T18, EDSS_T12, EDSS_T6, EDSS_T0)
      x[which(!is.na(x))[1]]
    }
  ) %>%
  ungroup()


# Peggioramento EDSS
merged_naive <- merged_naive %>%
  mutate(
    EDSS_worsening = case_when(
      EDSS_T0 == 0 ~ (EDSS_last - EDSS_T0) >= 1.5,
      EDSS_T0 >= 1 & EDSS_T0 <= 5 ~ (EDSS_last - EDSS_T0) >= 1,
      EDSS_T0 > 5 ~ (EDSS_last - EDSS_T0) >= 0.5,
      TRUE ~ FALSE
    )
  )


# Eventi NEDA
merged_naive <- merged_naive %>%
  mutate(
    event_6_C  = EDSS_worsening | relapse_T6  == 1 | lesioni_T6  == 1,
    event_12_C = EDSS_worsening | relapse_T12 == 1 | lesioni_T12 == 1,
    event_18_C = EDSS_worsening | relapse_T18 == 1 | lesioni_T18 == 1,
    event_24_C = EDSS_worsening | relapse_T24 == 1 | lesioni_T24 == 1,
    event_30_C = EDSS_worsening | relapse_T30 == 1 | lesioni_T30 == 1,
    event_36_C = EDSS_worsening | relapse_T36 == 1 | lesioni_T36 == 1
  )

# Tempo al primo evento
merged_naive <- merged_naive %>%
  mutate(
    event_time_C = case_when(
      event_6_C  ~ 6,
      event_12_C ~ 12,
      event_18_C ~ 18,
      event_24_C ~ 24,
      event_30_C ~ 30,
      event_36_C ~ 36,
      TRUE ~ NA_real_),
    event_status_C = ifelse(is.na(event_time_C), 0, 1),
    time_C = ifelse(is.na(event_time_C), 36, event_time_C))

surv_obj4 <- Surv(time = merged_naive$time_C, event = merged_naive$event_status_C)
fit4 <- survfit(surv_obj4 ~ GROUP, data = merged_naive)


ggsurvplot(
  fit4,
  data = merged_naive,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  title = "NEDA-3 free survival curves - Naive patients",
  legend.title = "Groups",
  legend.labs = c("MAV", "TEC"),
  xlab = "Months",
  ylab = "NEDA-3 loss probability"
)

