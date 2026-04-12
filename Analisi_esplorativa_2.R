# ANALISI METADATA INIZIALE
# merged_data = read.csv("C:/Users/gialu/Desktop/tesi_magistrale/Git.hub/Thesis_crepaldi/Input/merged_data.csv")
# merged_data <- read.csv("Input/merged_data.csv") # DT: dato che abbiamo la repo sincronizzata, puoi fare riferimento al path relativo, in modo che lo script funzioni indipendentemente da dove lo apri
getwd()
setwd("C:/Users/gialu/Desktop/tesi_magistrale/Git.hub/Thesis_crepaldi")
merged_data <- read.csv("Input/merged_data.csv")

str(merged_data)

# metto a posto la tipologia di variabile di nuovo:
# le date
cols_date <- c("Data_di_nascita", "Data_esordio", "Data_diagnosi", "Data_prelievo_T0", "data_referto_T6", "data_referto_T12", "data_referto_T18", "data_referto_T24", "data_referto_T30", "data_referto_T36", "data_cambio_terapia", "sample_collection_date", "glc_date_before_t0")
merged_data[cols_date] <- lapply(merged_data[cols_date], as.Date)

# categorico
cols_factor <- c("lesioni_T6", "relapse_T6", "lesioni_T12", "relapse_T12", "lesioni_T18", "relapse_T18", "lesioni_T24", "relapse_T24", "lesioni_T30", "relapse_T30", "lesioni_T36", "relapse_T36", "cambio_terapia", "lesion_burden", "bone_marrow_lesions", "subtentorial_lesions", "gadolinium_contrast", "full_recovery")
merged_data[cols_factor] <- lapply(merged_data[cols_factor], factor, levels = c("No", "Yes"))
cols_factor2 <- c("naive", "smoking_habit", "antibiotic_use", "spike_in")
merged_data[cols_factor2] <- lapply(merged_data[cols_factor2], factor, levels = c("Yes", "No"))
cols_factor3 <- c("glucocorticoid_treatment", "gc_treatment", "prognosis_at_onset")
merged_data[cols_factor3] <- lapply(merged_data[cols_factor3], factor, levels = c("Positive", "Negative"))
merged_data$physical_activity <- factor(merged_data$physical_activity, levels = c("active", "sedentary", "partially_active"))
merged_data$simplified_batch <- factor(merged_data$simplified_batch, levels = c("1-3", "4-5"))
merged_data$disease_activity <- factor(merged_data$disease_activity, levels = c("higher_activity", "lower_activity"))
merged_data$sex <- factor(merged_data$sex, levels = c("M", "F"))
# DT: rendi factor anche EDSS, CENTRO, sequencing batch, clinical presentation

merged_data$sequencing_batch <- factor(merged_data$sequencing_batch, levels = c("1", "2", "3"))
merged_data$CENTRO <- factor(merged_data$CENTRO, levels = c("CLERICO", "NAPOLI", "MOLINETTE", "CRESM"))
merged_data$clinical_presentation <- factor(merged_data$clinical_presentation)
merged_data <- merged_data %>%
  mutate(
    across(starts_with("EDSS"), as.factor)
  )

# numerici
cols_numeric <- c("MEDAS", "cortisone_T6", "cortisone_T12", "cortisone_T18", "cortisone_T24", "cortisone_T30", "cortisone_T36", "age_at_t0", "bmi", "pyr_mds", "dna_quantification", "days_gc_collection", "percent_th17", "treg_il10_pos", "treg_cd39_pos", "age", "age_at_diagnosis")
merged_data[cols_numeric] <- lapply(merged_data[cols_numeric], as.numeric)

#####Inizio analisi####
# analisi esplorativa min,max, media, mediana di alcuni dati
summary(merged_data[sapply(merged_data, is.numeric)])

lapply(merged_data[sapply(merged_data, is.factor)], table, useNA = "always")

# analisi delle possibili correlazioni tra variabili:

# EDSS vs smoking_habit
t.test(EDSS_T0 ~ smoking_habit, data = merged_data) # differenza non significativa
t.test(EDSS_T6 ~ smoking_habit, data = merged_data) # differenza non significativa
t.test(EDSS_T12 ~ smoking_habit, data = merged_data) # differenza non significativa
t.test(EDSS_T18 ~ smoking_habit, data = merged_data) # differenza non significativa
t.test(EDSS_T24 ~ smoking_habit, data = merged_data) # differenza non significativa
t.test(EDSS_T30 ~ smoking_habit, data = merged_data) # differenza non significativa
t.test(EDSS_T36 ~ smoking_habit, data = merged_data) # non ci sono abbastanza osservazioni

library(ggplot2)
ggplot(merged_data, aes(x = smoking_habit, y = EDSS_T0, fill = smoking_habit)) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal()
ggplot(merged_data, aes(x = smoking_habit, y = EDSS_T6, fill = smoking_habit)) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal()
ggplot(merged_data, aes(x = smoking_habit, y = EDSS_T12, fill = smoking_habit)) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal()
ggplot(merged_data, aes(x = smoking_habit, y = EDSS_T18, fill = smoking_habit)) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal()
ggplot(merged_data, aes(x = smoking_habit, y = EDSS_T24, fill = smoking_habit)) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal()
ggplot(merged_data, aes(x = smoking_habit, y = EDSS_T30, fill = smoking_habit)) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal()
ggplot(merged_data, aes(x = smoking_habit, y = EDSS_T36, fill = smoking_habit)) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal()

# EDSS vs sex
ggplot(subset(merged_data, !is.na(sex) & !is.na(EDSS_T0)), aes(x = sex, y = EDSS_T0, fill = sex)) +
  geom_boxplot() +
  theme_minimal()
ggplot(subset(merged_data, !is.na(sex) & !is.na(EDSS_T6)), aes(x = sex, y = EDSS_T6, fill = sex)) +
  geom_boxplot() +
  theme_minimal()
ggplot(subset(merged_data, !is.na(sex) & !is.na(EDSS_T12)), aes(x = sex, y = EDSS_T12, fill = sex)) +
  geom_boxplot() +
  theme_minimal()
ggplot(subset(merged_data, !is.na(sex) & !is.na(EDSS_T18)), aes(x = sex, y = EDSS_T18, fill = sex)) +
  geom_boxplot() +
  theme_minimal()
ggplot(subset(merged_data, !is.na(sex) & !is.na(EDSS_T24)), aes(x = sex, y = EDSS_T24, fill = sex)) +
  geom_boxplot() +
  theme_minimal()
ggplot(subset(merged_data, !is.na(sex) & !is.na(EDSS_T30)), aes(x = sex, y = EDSS_T30, fill = sex)) +
  geom_boxplot() +
  theme_minimal()
ggplot(subset(merged_data, !is.na(sex) & !is.na(EDSS_T36)), aes(x = sex, y = EDSS_T36, fill = sex)) +
  geom_boxplot() +
  theme_minimal()

# EDSS vs physical_activity
ggplot(subset(merged_data, !is.na(physical_activity) & !is.na(EDSS_T0)), aes(x = physical_activity, y = EDSS_T0, fill = physical_activity)) +
  geom_boxplot() +
  theme_minimal()
ggplot(subset(merged_data, !is.na(physical_activity) & !is.na(EDSS_T6)), aes(x = physical_activity, y = EDSS_T6, fill = physical_activity)) +
  geom_boxplot() +
  theme_minimal()
ggplot(subset(merged_data, !is.na(physical_activity) & !is.na(EDSS_T12)), aes(x = physical_activity, y = EDSS_T12, fill = physical_activity)) +
  geom_boxplot() +
  theme_minimal()
ggplot(subset(merged_data, !is.na(physical_activity) & !is.na(EDSS_T18)), aes(x = physical_activity, y = EDSS_T18, fill = physical_activity)) +
  geom_boxplot() +
  theme_minimal()
ggplot(subset(merged_data, !is.na(physical_activity) & !is.na(EDSS_T24)), aes(x = physical_activity, y = EDSS_T24, fill = physical_activity)) +
  geom_boxplot() +
  theme_minimal()
ggplot(subset(merged_data, !is.na(physical_activity) & !is.na(EDSS_T30)), aes(x = physical_activity, y = EDSS_T30, fill = physical_activity)) +
  geom_boxplot() +
  theme_minimal()
ggplot(subset(merged_data, !is.na(physical_activity) & !is.na(EDSS_T36)), aes(x = physical_activity, y = EDSS_T36, fill = physical_activity)) +
  geom_boxplot() +
  theme_minimal()

# EDSS vs age_at_diagnosis
cor.test(merged_data$age_at_diagnosis, merged_data$EDSS_T0) # correlazione non significativa
cor.test(merged_data$age_at_diagnosis, merged_data$EDSS_T6) # correlazione non significativa
cor.test(merged_data$age_at_diagnosis, merged_data$EDSS_T12) # correlazione non significativa
cor.test(merged_data$age_at_diagnosis, merged_data$EDSS_T18) # correlazione non significativa
cor.test(merged_data$age_at_diagnosis, merged_data$EDSS_T24) # correlazione positiva SIGNIFICATIVA
cor.test(merged_data$age_at_diagnosis, merged_data$EDSS_T30) # correlazione non significativa
cor.test(merged_data$age_at_diagnosis, merged_data$EDSS_T36) # correlazione non significativa

# linee di regressione EDSS vs age_at_diagnosis
ggplot(merged_data, aes(x = age_at_diagnosis, y = EDSS_T0)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()
ggplot(merged_data, aes(x = age_at_diagnosis, y = EDSS_T6)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()
ggplot(merged_data, aes(x = age_at_diagnosis, y = EDSS_T12)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()
ggplot(merged_data, aes(x = age_at_diagnosis, y = EDSS_T18)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()
ggplot(merged_data, aes(x = age_at_diagnosis, y = EDSS_T24)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()
ggplot(merged_data, aes(x = age_at_diagnosis, y = EDSS_T30)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()
ggplot(merged_data, aes(x = age_at_diagnosis, y = EDSS_T36)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()


# EDSS vs bmi
cor.test(merged_data$bmi, merged_data$EDSS_T0) # correlazione non significativa
cor.test(merged_data$bmi, merged_data$EDSS_T6) # correlazione non significativa
cor.test(merged_data$bmi, merged_data$EDSS_T12) # correlazione non significativa
cor.test(merged_data$bmi, merged_data$EDSS_T18) # correlazione non significativa
cor.test(merged_data$bmi, merged_data$EDSS_T24) # correlazione non significativa
cor.test(merged_data$bmi, merged_data$EDSS_T30) # correlazione non significativa
cor.test(merged_data$bmi, merged_data$EDSS_T36) # correlazione non significativa

# linee di regressione bmi vs EDSS
ggplot(merged_data, aes(x = bmi, y = EDSS_T0)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()
ggplot(merged_data, aes(x = bmi, y = EDSS_T6)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()
ggplot(merged_data, aes(x = bmi, y = EDSS_T12)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()
ggplot(merged_data, aes(x = bmi, y = EDSS_T18)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()
ggplot(merged_data, aes(x = bmi, y = EDSS_T24)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()
ggplot(merged_data, aes(x = bmi, y = EDSS_T30)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()
ggplot(merged_data, aes(x = bmi, y = EDSS_T36)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()


# EDSS vs MEDAS
cor.test(merged_data$MEDAS, merged_data$EDSS_T0) # correlazione non significativa
cor.test(merged_data$MEDAS, merged_data$EDSS_T6) # correlazione non significativa
cor.test(merged_data$MEDAS, merged_data$EDSS_T12) # correlazione non significativa
cor.test(merged_data$MEDAS, merged_data$EDSS_T18) # correlazione non significativa
cor.test(merged_data$MEDAS, merged_data$EDSS_T24) # correlazione non significativa
cor.test(merged_data$MEDAS, merged_data$EDSS_T30) # correlazione non significativa
cor.test(merged_data$MEDAS, merged_data$EDSS_T36) # non ci sono abbastanza osservazioni

# linee di regressione EDSS vs MEDAS
ggplot(merged_data, aes(x = EDSS_T0, y = MEDAS)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()
ggplot(merged_data, aes(x = EDSS_T6, y = MEDAS)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()
ggplot(merged_data, aes(x = EDSS_T12, y = MEDAS)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()
ggplot(merged_data, aes(x = EDSS_T18, y = MEDAS)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()
ggplot(merged_data, aes(x = EDSS_T24, y = MEDAS)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()
ggplot(merged_data, aes(x = EDSS_T30, y = MEDAS)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()
ggplot(merged_data, aes(x = EDSS_T36, y = MEDAS)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()

# age_at_diagnosis vs relapse
t.test(age_at_diagnosis ~ relapse_T6, data = merged_data) # correlazione non significativa
t.test(age_at_diagnosis ~ relapse_T24, data = merged_data) # correlazione non significativa
t.test(age_at_diagnosis ~ relapse_T30, data = merged_data) # significativo

# age_at_diagnosis vs lesions
t.test(age_at_diagnosis ~ lesioni_T6, data = merged_data) # correlazione non significativa
t.test(age_at_diagnosis ~ lesioni_T12, data = merged_data) # correlazione non significativa
t.test(age_at_diagnosis ~ lesioni_T18, data = merged_data) # correlazione non significativa
t.test(age_at_diagnosis ~ lesioni_T24, data = merged_data) # correlazione non significativa
t.test(age_at_diagnosis ~ lesioni_T30, data = merged_data) # correlazione non significativa
t.test(age_at_diagnosis ~ lesioni_T36, data = merged_data) # non ci sono abbastanza osservazioni

# bmi vs relapse
t.test(bmi ~ relapse_T6, data = merged_data) # correlazione non significativa
t.test(bmi ~ relapse_T24, data = merged_data) # correlazione non significativa
t.test(bmi ~ relapse_T30, data = merged_data) # correlazione non significativa

# bmi vs new lesions
t.test(bmi ~ lesioni_T6, data = merged_data) # correlazione non significativa
t.test(bmi ~ lesioni_T12, data = merged_data) # correlazione non significativa
t.test(bmi ~ lesioni_T18, data = merged_data) # correlazione non significativa
t.test(bmi ~ lesioni_T24, data = merged_data) # correlazione non significativa
t.test(bmi ~ lesioni_T30, data = merged_data) # correlazione non significativa


# distribuzione di singole variabili con grafico

# age_at_diagnosis
ggplot(merged_data, aes(x = age_at_diagnosis)) +
  geom_histogram(binwidth = 0.5, fill = "#4C78A8", color = "white", linewidth = 0.4) +
  scale_x_continuous(breaks = seq(15, 60, 10)) +
  theme_minimal() +
  labs(title = "Distribution age_at_diagnosis", x = "age_at_diagnosis", y = "Frequency")

# EDSS
ggplot(merged_data, aes(x = EDSS_T0)) +
  geom_histogram(binwidth = 0.5, fill = "#4C78A8", color = "white", linewidth = 0.4) +
  scale_x_continuous(breaks = seq(0, 10, 0.5)) +
  theme_minimal() +
  labs(title = "Distribution EDSS_T0", x = "EDSS_T0", y = "Frequency")
ggplot(merged_data, aes(x = EDSS_T6)) +
  geom_histogram(binwidth = 0.5, fill = "#4C78A8", color = "white", linewidth = 0.4) +
  scale_x_continuous(breaks = seq(0, 10, 0.5)) +
  theme_minimal() +
  labs(title = "Distribution EDSS_T6", x = "EDSS_T6", y = "Frequency")
ggplot(merged_data, aes(x = EDSS_T12)) +
  geom_histogram(binwidth = 0.5, fill = "#4C78A8", color = "white", linewidth = 0.4) +
  scale_x_continuous(breaks = seq(0, 10, 0.5)) +
  theme_minimal() +
  labs(title = "Distribution EDSS_T12", x = "EDSS_T12", y = "Frequency")
ggplot(merged_data, aes(x = EDSS_T18)) +
  geom_histogram(binwidth = 0.5, fill = "#4C78A8", color = "white", linewidth = 0.4) +
  scale_x_continuous(breaks = seq(0, 10, 0.5)) +
  theme_minimal() +
  labs(title = "Distribution EDSS_T18", x = "EDSS_T18", y = "Frequency")
ggplot(merged_data, aes(x = EDSS_T24)) +
  geom_histogram(binwidth = 0.5, fill = "#4C78A8", color = "white", linewidth = 0.4) +
  scale_x_continuous(breaks = seq(0, 10, 0.5)) +
  theme_minimal() +
  labs(title = "Distribution EDSS_T24", x = "EDSS_T24", y = "Frequency")
ggplot(merged_data, aes(x = EDSS_T30)) +
  geom_histogram(binwidth = 0.5, fill = "#4C78A8", color = "white", linewidth = 0.4) +
  scale_x_continuous(breaks = seq(0, 10, 0.5)) +
  theme_minimal() +
  labs(title = "Distribution EDSS_T30", x = "EDSS_T30", y = "Frequency")
ggplot(merged_data, aes(x = EDSS_T36)) +
  geom_histogram(binwidth = 0.5, fill = "#4C78A8", color = "white", linewidth = 0.4) +
  scale_x_continuous(breaks = seq(0, 10, 0.5)) +
  theme_minimal() +
  labs(title = "Distribution EDSS_T36", x = "EDSS_T36", y = "Frequency")


# bmi
ggplot(merged_data, aes(x = bmi)) +
  geom_histogram(binwidth = 0.5, fill = "#4C78A8", color = "white", linewidth = 0.4) +
  scale_x_continuous(breaks = seq(15, 35, 2)) +
  theme_minimal()
# DT: per una variabile numerica continua come il BMI è meglio usare un boxplot per visualizzare la distribuzione; in alternativa puoi discretizzare il bmi in, ad esempio, tre classi (sottopeso, normopeso, sovrappeso/obeso) e vedere la frequenza con un grafico a barre
ggplot(merged_data, aes(y = bmi)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplot_BMI",y = "BMI",x = "")



# MEDAS
ggplot(merged_data, aes(x = MEDAS)) +
  geom_histogram(binwidth = 0.5, fill = "#4C78A8", color = "white", linewidth = 0.4) +
  scale_x_continuous(breaks = seq(20, 100, 5)) +
  theme_minimal()

# grafico andamento di alcune variabili
# EDSS nel tempo
plot(c(0, 6, 12, 18, 24, 30, 36), c(mean(merged_data$EDSS_T0, na.rm = TRUE), mean(merged_data$EDSS_T6, na.rm = TRUE), mean(merged_data$EDSS_T12, na.rm = TRUE), mean(merged_data$EDSS_T18, na.rm = TRUE), mean(merged_data$EDSS_T24, na.rm = TRUE), mean(merged_data$EDSS_T30, na.rm = TRUE), mean(merged_data$EDSS_T36, na.rm = TRUE)), type = "b", xlab = "Tempo", ylab = "EDSS_medio", main = "Andamento medio EDSS nel tempo")

# grafico che associa EDSS nei vari time point e smoking habit
library(dplyr)
library(tidyr)

edss_long <- merged_data %>%
  select(ID, smoking_habit, EDSS_T0, EDSS_T6, EDSS_T12, EDSS_T24, EDSS_T36) %>%
  pivot_longer(cols = starts_with("EDSS"), names_to = "time", values_to = "EDSS") %>%
  mutate(time = factor(time, levels = c("EDSS_T0", "EDSS_T6", "EDSS_T12", "EDSS_T24", "EDSS_T36")))

ggplot(edss_long, aes(x = time, y = EDSS, color = smoking_habit, group = smoking_habit)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  theme_minimal()
# DT: qui puoi migliorare il grafico modificando i livelli di Time con il numero del mese.

# grafico che associa EDSS nei time point e sex
edss_sex <- merged_data %>%
  select(ID, sex, EDSS_T0, EDSS_T6, EDSS_T12, EDSS_T24, EDSS_T36) %>%
  pivot_longer(cols = starts_with("EDSS"), names_to = "time", values_to = "EDSS") %>%
  mutate(time = factor(time, levels = c("EDSS_T0", "EDSS_T6", "EDSS_T12", "EDSS_T24", "EDSS_T36")))

ggplot(edss_sex, aes(x = time, y = EDSS, color = sex, group = sex)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  theme_minimal()

# grafico che associa EDSS nei time point e naive
edss_naive <- merged_data %>%
  select(ID, naive, EDSS_T0, EDSS_T6, EDSS_T12, EDSS_T24, EDSS_T36) %>%
  pivot_longer(cols = starts_with("EDSS"), names_to = "time", values_to = "EDSS") %>%
  mutate(time = factor(time, levels = c("EDSS_T0", "EDSS_T6", "EDSS_T12", "EDSS_T24", "EDSS_T36")))

ggplot(edss_naive, aes(x = time, y = EDSS, color = naive, group = naive)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  theme_minimal()

# grafico che associa EDSS nei time point e physical_activity
edss_activity <- merged_data %>%
  select(ID, physical_activity, EDSS_T0, EDSS_T6, EDSS_T12, EDSS_T24, EDSS_T36) %>%
  pivot_longer(cols = starts_with("EDSS"), names_to = "time", values_to = "EDSS") %>%
  mutate(time = factor(time, levels = c("EDSS_T0", "EDSS_T6", "EDSS_T12", "EDSS_T24", "EDSS_T36")))

ggplot(edss_activity, aes(x = time, y = EDSS, color = physical_activity, group = physical_activity)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  theme_minimal()
