# ANALISI METADATA INIZIALE
# merged_data = read.csv("C:/Users/gialu/Desktop/tesi_magistrale/Git.hub/Thesis_crepaldi/Input/merged_data.csv")
# merged_data <- read.csv("Input/merged_data.csv") # DT: dato che abbiamo la repo sincronizzata, puoi fare riferimento al path relativo, in modo che lo script funzioni indipendentemente da dove lo apri
source("Libraries.R")
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
    across(
      starts_with("EDSS"),
      ~ ordered(as.numeric(.x), levels = seq(1, 10, by = 0.5))
    )
  )

# numerici
cols_numeric <- c("MEDAS", "cortisone_T6", "cortisone_T12", "cortisone_T18", "cortisone_T24", "cortisone_T30", "cortisone_T36", "age_at_t0", "bmi", "dna_quantification", "days_gc_collection", "percent_th17", "treg_il10_pos", "treg_cd39_pos", "age")
merged_data[cols_numeric] <- lapply(merged_data[cols_numeric], as.numeric)

#####Inizio analisi####
# analisi esplorativa min,max, media, mediana di alcuni dati
summary(merged_data[sapply(merged_data, is.numeric)])

lapply(merged_data[sapply(merged_data, is.factor)], table, useNA = "always")

#analisi NA globale
na_count <- sapply(merged_data, function(x) sum(is.na(x)))
na_percent <- sapply(merged_data, function(x) mean(is.na(x)) * 100)

na_summary <- data.frame(variable = names(na_count), n_NA = na_count, perc_NA = round(na_percent, 2))
na_summary[order(-na_summary$perc_NA), ]

gg_miss_var(merged_data) + theme_minimal() + labs(title = "Missing values per variable")

# analisi delle possibili correlazioni tra variabili:

# EDSS vs smoking_habit
wilcox.test(as.numeric(as.character(EDSS_T0)) ~ smoking_habit, data = merged_data, exact = FALSE) # differenza non significativa
wilcox.test(as.numeric(as.character(EDSS_T6)) ~ smoking_habit, data = merged_data, exact = FALSE) # differenza non significativa
wilcox.test(as.numeric(as.character(EDSS_T12)) ~ smoking_habit, data = merged_data, exact = FALSE) # differenza non significativa
wilcox.test(as.numeric(as.character(EDSS_T18)) ~ smoking_habit, data = merged_data, exact = FALSE) # differenza non significativa
wilcox.test(as.numeric(as.character(EDSS_T24)) ~ smoking_habit, data = merged_data, exact = FALSE) # differenza non significativa
wilcox.test(as.numeric(as.character(EDSS_T30)) ~ smoking_habit, data = merged_data, exact = FALSE) # differenza non significativa

ggplot(subset(merged_data, !is.na(smoking_habit) & !is.na(EDSS_T0)), aes(x = smoking_habit, y = as.numeric(as.character(EDSS_T0)), fill = smoking_habit)) + geom_boxplot(alpha = 0.8) + theme_minimal() + labs(title = "EDSS vs smoking_habit", x = "smoking_habit", y = "EDSS_T0")
ggplot(subset(merged_data, !is.na(smoking_habit) & !is.na(EDSS_T6)), aes(x = smoking_habit, y = as.numeric(as.character(EDSS_T6)), fill = smoking_habit)) + geom_boxplot(alpha = 0.8) + theme_minimal() + labs(title = "EDSS vs smoking_habit", x = "smoking_habit", y = "EDSS_T6")
ggplot(subset(merged_data, !is.na(smoking_habit) & !is.na(EDSS_T12)), aes(x = smoking_habit, y = as.numeric(as.character(EDSS_T12)), fill = smoking_habit)) + geom_boxplot(alpha = 0.8) + theme_minimal() + labs(title = "EDSS vs smoking_habit", x = "smoking_habit", y = "EDSS_T12")
ggplot(subset(merged_data, !is.na(smoking_habit) & !is.na(EDSS_T18)), aes(x = smoking_habit, y = as.numeric(as.character(EDSS_T18)), fill = smoking_habit)) + geom_boxplot(alpha = 0.8) + theme_minimal() + labs(title = "EDSS vs smoking_habit", x = "smoking_habit", y = "EDSS_T18")
ggplot(subset(merged_data, !is.na(smoking_habit) & !is.na(EDSS_T24)), aes(x = smoking_habit, y = as.numeric(as.character(EDSS_T24)), fill = smoking_habit)) + geom_boxplot(alpha = 0.8) + theme_minimal() + labs(title = "EDSS vs smoking_habit", x = "smoking_habit", y = "EDSS_T24")
ggplot(subset(merged_data, !is.na(smoking_habit) & !is.na(EDSS_T30)), aes(x = smoking_habit, y = as.numeric(as.character(EDSS_T30)), fill = smoking_habit)) + geom_boxplot(alpha = 0.8) + theme_minimal() + labs(title = "EDSS vs smoking_habit", x = "smoking_habit", y = "EDSS_T30")


# EDSS vs sex
wilcox.test(as.numeric(as.character(EDSS_T0)) ~ sex, data = merged_data, exact = FALSE) # differenza non significativa
wilcox.test(as.numeric(as.character(EDSS_T12)) ~ sex, data = merged_data, exact = FALSE) # differenza non significativa


ggplot(subset(merged_data, !is.na(sex) & !is.na(EDSS_T0)), aes(x = sex, as.numeric(as.character(EDSS_T0)), fill = sex)) +geom_boxplot() + theme_minimal() + labs(title = "EDSS vs Sex", x = "Sex", y = "EDSS_T0")
ggplot(subset(merged_data, !is.na(sex) & !is.na(EDSS_T6)), aes(x = sex, as.numeric(as.character(EDSS_T6)), fill = sex)) +geom_boxplot() + theme_minimal() + labs(title = "EDSS vs Sex", x = "Sex", y = "EDSS_T6")
ggplot(subset(merged_data, !is.na(sex) & !is.na(EDSS_T12)), aes(x = sex, as.numeric(as.character(EDSS_T12)), fill = sex)) +geom_boxplot() + theme_minimal() + labs(title = "EDSS vs Sex", x = "Sex", y = "EDSS_T12")
ggplot(subset(merged_data, !is.na(sex) & !is.na(EDSS_T18)), aes(x = sex, as.numeric(as.character(EDSS_T18)), fill = sex)) +geom_boxplot() + theme_minimal() + labs(title = "EDSS vs Sex", x = "Sex", y = "EDSS_T18")
ggplot(subset(merged_data, !is.na(sex) & !is.na(EDSS_T24)), aes(x = sex, as.numeric(as.character(EDSS_T24)), fill = sex)) +geom_boxplot() + theme_minimal() + labs(title = "EDSS vs Sex", x = "Sex", y = "EDSS_T24")
ggplot(subset(merged_data, !is.na(sex) & !is.na(EDSS_T30)), aes(x = sex, as.numeric(as.character(EDSS_T30)), fill = sex)) +geom_boxplot() + theme_minimal() + labs(title = "EDSS vs Sex", x = "Sex", y = "EDSS_T30")


# EDSS vs physical_activity
ggplot(subset(merged_data, !is.na(physical_activity) & !is.na(EDSS_T0)), aes(x = physical_activity, y = as.numeric(as.character(EDSS_T0)), fill = physical_activity)) + geom_boxplot() + theme_minimal() + labs(title = "EDSS vs Physical_activity", x = "Physical_activity", y = "EDSS_T0")
ggplot(subset(merged_data, !is.na(physical_activity) & !is.na(EDSS_T6)), aes(x = physical_activity, y = as.numeric(as.character(EDSS_T6)), fill = physical_activity)) + geom_boxplot() + theme_minimal() + labs(title = "EDSS vs Physical_activity", x = "Physical_activity", y = "EDSS_T6")
ggplot(subset(merged_data, !is.na(physical_activity) & !is.na(EDSS_T12)), aes(x = physical_activity, y = as.numeric(as.character(EDSS_T12)), fill = physical_activity)) + geom_boxplot() + theme_minimal() + labs(title = "EDSS vs Physical_activity", x = "Physical_activity", y = "EDSS_T12")
ggplot(subset(merged_data, !is.na(physical_activity) & !is.na(EDSS_T18)), aes(x = physical_activity, y = as.numeric(as.character(EDSS_T18)), fill = physical_activity)) + geom_boxplot() + theme_minimal() + labs(title = "EDSS vs Physical_activity", x = "Physical_activity", y = "EDSS_T18")
ggplot(subset(merged_data, !is.na(physical_activity) & !is.na(EDSS_T24)), aes(x = physical_activity, y = as.numeric(as.character(EDSS_T24)), fill = physical_activity)) + geom_boxplot() + theme_minimal() + labs(title = "EDSS vs Physical_activity", x = "Physical_activity", y = "EDSS_T24")
ggplot(subset(merged_data, !is.na(physical_activity) & !is.na(EDSS_T30)), aes(x = physical_activity, y = as.numeric(as.character(EDSS_T30)), fill = physical_activity)) + geom_boxplot() + theme_minimal() + labs(title = "EDSS vs Physical_activity", x = "Physical_activity", y = "EDSS_T30")


# EDSS vs age_at_t0
cor.test(merged_data$age_at_t0, as.numeric(merged_data$EDSS_T0), method = "spearman", exact = FALSE) # non significativo
cor.test(merged_data$age_at_t0, as.numeric(merged_data$EDSS_T6), method = "spearman", exact = FALSE) # non significativo
cor.test(merged_data$age_at_t0, as.numeric(merged_data$EDSS_T12), method = "spearman", exact = FALSE) # non significativo
cor.test(merged_data$age_at_t0, as.numeric(merged_data$EDSS_T18), method = "spearman", exact = FALSE) # non significativo
cor.test(merged_data$age_at_t0, as.numeric(merged_data$EDSS_T24), method = "spearman", exact = FALSE) #SIGNIFICATIVO
cor.test(merged_data$age_at_t0, as.numeric(merged_data$EDSS_T30), method = "spearman", exact = FALSE) #SIGNIFICATIVO

# linee di regressione EDSS vs age_at_t0
ggplot(merged_data, aes(x = age_at_t0, y = as.numeric(as.character(EDSS_T0)))) + geom_point(size = 3, alpha = 0.8) + geom_smooth(method = "lm", se = TRUE) + theme_minimal() + labs(title = "EDSS vs Age at T0", x = "Age_at_T0", y = "EDSS_T0")
ggplot(merged_data, aes(x = age_at_t0, y = as.numeric(as.character(EDSS_T12)))) + geom_point(size = 3, alpha = 0.8) + geom_smooth(method = "lm", se = TRUE) + theme_minimal() + labs(title = "EDSS vs Age at T0", x = "Age_at_T0", y = "EDSS_T12")
ggplot(merged_data, aes(x = age_at_t0, y = as.numeric(as.character(EDSS_T24)))) + geom_point(size = 3, alpha = 0.8) + geom_smooth(method = "lm", se = TRUE) + theme_minimal() + labs(title = "EDSS vs Age at T0", x = "Age_at_T0", y = "EDSS_T24")


# EDSS vs bmi
cor.test(merged_data$bmi, as.numeric(merged_data$EDSS_T0), method = "spearman", exact = FALSE) # non significativo
cor.test(merged_data$bmi, as.numeric(merged_data$EDSS_T6), method = "spearman", exact = FALSE) # non significativo
cor.test(merged_data$bmi, as.numeric(merged_data$EDSS_T12), method = "spearman", exact = FALSE) # non significativo
cor.test(merged_data$bmi, as.numeric(merged_data$EDSS_T24), method = "spearman", exact = FALSE) # non significativo
cor.test(merged_data$bmi, as.numeric(merged_data$EDSS_T30), method = "spearman", exact = FALSE) # non significativo

# EDSS vs MEDAS
cor.test(merged_data$MEDAS, as.numeric(merged_data$EDSS_T0), method = "spearman", exact = FALSE) # non significativo
cor.test(merged_data$MEDAS, as.numeric(merged_data$EDSS_T6), method = "spearman", exact = FALSE) # non significativo
cor.test(merged_data$MEDAS, as.numeric(merged_data$EDSS_T12), method = "spearman", exact = FALSE) # non significativo
cor.test(merged_data$MEDAS, as.numeric(merged_data$EDSS_T18), method = "spearman", exact = FALSE) # non significativo
cor.test(merged_data$MEDAS, as.numeric(merged_data$EDSS_T24), method = "spearman", exact = FALSE) # non significativo
cor.test(merged_data$MEDAS, as.numeric(merged_data$EDSS_T30), method = "spearman", exact = FALSE) # non significativo

# age_at_t0 vs relapse
t.test(age_at_t0 ~ lesioni_T6, data = merged_data) # correlazione SIGNIFICATIVA (di poco)
t.test(age_at_t0 ~ lesioni_T12, data = merged_data) # correlazione SIGNIFICATIVA
t.test(age_at_t0 ~ lesioni_T18, data = merged_data) # correlazione non significativa
t.test(age_at_t0 ~ lesioni_T24, data = merged_data) # correlazione non significativa
t.test(age_at_t0 ~ lesioni_T30, data = merged_data) # correlazione non significativa

# age_at_t0 vs lesions
t.test(age_at_t0 ~ lesioni_T6, data = merged_data) # correlazione SIGNIFICATIVA (di poco)
t.test(age_at_t0 ~ lesioni_T12, data = merged_data) # correlazione SIGNIFICATIVA
t.test(age_at_t0 ~ lesioni_T18, data = merged_data) # correlazione non significativa
t.test(age_at_t0 ~ lesioni_T24, data = merged_data) # correlazione non significativa
t.test(age_at_t0 ~ lesioni_T30, data = merged_data) # correlazione non significativa

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

# age_at_t0
ggplot(merged_data, aes(x = age_at_t0)) +
  geom_histogram(binwidth = 0.5, fill = "#4C78A8", color = "white", linewidth = 0.4) +
  scale_x_continuous(breaks = seq(15, 60, 10)) +
  theme_minimal() +
  labs(title = "Distribution age_at_T0", x = "age_at_t0", y = "Frequency")

# EDSS (T30 E T36 HANNO TROPPI NA)
ggplot(merged_data, aes(x = EDSS_T0)) + geom_bar(fill = "#4C78A8", color = "white", linewidth = 0.4) + theme_minimal() + labs(title = "Distribution EDSS_T0", x = "EDSS_T0", y = "Frequency")
ggplot(merged_data, aes(x = EDSS_T6)) + geom_bar(fill = "#4C78A8", color = "white", linewidth = 0.4) + theme_minimal() + labs(title = "Distribution EDSS_T6", x = "EDSS_T6", y = "Frequency")
ggplot(merged_data, aes(x = EDSS_T12)) + geom_bar(fill = "#4C78A8", color = "white", linewidth = 0.4) + theme_minimal() + labs(title = "Distribution EDSS_T12", x = "EDSS_T12", y = "Frequency")
ggplot(merged_data, aes(x = EDSS_T18)) + geom_bar(fill = "#4C78A8", color = "white", linewidth = 0.4) + theme_minimal() + labs(title = "Distribution EDSS_T18", x = "EDSS_T18", y = "Frequency")
ggplot(merged_data, aes(x = EDSS_T24)) + geom_bar(fill = "#4C78A8", color = "white", linewidth = 0.4) + theme_minimal() + labs(title = "Distribution EDSS_T24", x = "EDSS_T24", y = "Frequency")
ggplot(merged_data, aes(x = EDSS_T30)) + geom_bar(fill = "#4C78A8", color = "white", linewidth = 0.4) + theme_minimal() + labs(title = "Distribution EDSS_T30", x = "EDSS_T30", y = "Frequency")
ggplot(merged_data, aes(x = EDSS_T36)) + geom_bar(fill = "#4C78A8", color = "white", linewidth = 0.4) + theme_minimal() + labs(title = "Distribution EDSS_T36", x = "EDSS_T36", y = "Frequency")

# bmi
# DT: per una variabile numerica continua come il BMI è meglio usare un boxplot per visualizzare la distribuzione; in alternativa puoi discretizzare il bmi in, ad esempio, tre classi (sottopeso, normopeso, sovrappeso/obeso) e vedere la frequenza con un grafico a barre
ggplot(merged_data, aes(x = "", y = bmi)) + geom_boxplot(width = 0.25, outlier.shape = NA) + geom_jitter(width = 0.08, alpha = 0.5) + theme_minimal() +
  labs(title = "Distribution of BMI",x = "",y = "BMI")
#oppure discretizzando in categorie
merged_data$bmi_class <- factor(case_when(
    merged_data$bmi < 18.5 ~ "Underweight",
    merged_data$bmi >= 18.5 & merged_data$bmi < 25 ~ "Normal_weight",
    merged_data$bmi >= 25 & merged_data$bmi < 30 ~ "Overweight",
    merged_data$bmi >= 30 ~ "Obese",
    TRUE ~ NA_character_),
  levels = c("Underweight", "Normal_weight", "Overweight", "Obese"),ordered = TRUE)

ggplot(merged_data, aes(x = bmi_class)) + geom_bar(fill = "#4C78A8", color = "white") + theme_minimal()

#stesso grafico ma riempiti per il sesso
ggplot(merged_data, aes(x = bmi_class, fill = sex)) + geom_bar(position = "fill") + theme_minimal()

# MEDAS
ggplot(merged_data, aes(x = MEDAS)) + geom_histogram(bins = 10, fill = "#4C78A8", color = "white") + theme_minimal() + labs(title="MEDAS distribution")

#oppure come per il bmi discretizzo in categorie il MEDAS
summary(merged_data$MEDAS)
merged_data$MEDAS_class <- factor(case_when(
  merged_data$MEDAS < 40 ~ "Low_adherence",
  merged_data$MEDAS >= 40 & merged_data$MEDAS < 60 ~ "Moderate_adherence",
  merged_data$MEDAS >= 60 ~ "High_adherence",
  TRUE ~ NA_character_),
  levels = c("Low_adherence", "Moderate_adherence", "High_adherence"),ordered = TRUE)

ggplot(merged_data, aes(x = MEDAS_class)) + geom_bar(fill = "#4C78A8", color = "white") + theme_minimal() + labs(title = "MEDAS adherence categories", x = "Adherence level", y = "Frequency")

#preparo sottogruppo di dati per i grafici di andamento EDSS
# DT: qui puoi migliorare il grafico modificando i livelli di Time con il numero del mese.
#EDSS NEL TEMPO (non considero T36 perchè troppe poche osservazioni)
edss_long <- merged_data %>%
  select(ID, smoking_habit, sex, naive, physical_activity, EDSS_T0, EDSS_T6, EDSS_T12, EDSS_T18, EDSS_T24, EDSS_T30) %>%
  pivot_longer(cols = starts_with("EDSS"), names_to = "time", values_to = "EDSS_val") %>%
  mutate(EDSS_val = as.numeric(as.character(EDSS_val)), sex = trimws(as.character(sex)), time = factor(time, levels = c("EDSS_T0","EDSS_T6","EDSS_T12","EDSS_T18","EDSS_T24","EDSS_T30"))) %>%
  filter(!is.na(EDSS_val))

# 1 ANDAMENTO GENERALE
ggplot(edss_long, aes(x = time, y = EDSS_val, group = 1)) +
  stat_summary(fun = median, geom = "line", linewidth = 1.2, na.rm = TRUE) +
  stat_summary(fun = median, geom = "point", size = 3, na.rm = TRUE) +
  theme_minimal() + labs(title="EDSS during time", x = "Time", y = "Median EDSS")

# 2 ANDAMENTO PER SESSO 
ggplot(edss_long, aes(x = time, y = EDSS_val, color = sex, group = sex)) +
  stat_summary(fun = median, geom = "line", linewidth = 1.2, na.rm = TRUE) +
  stat_summary(fun = median, geom = "point", size = 3, na.rm = TRUE) +
  theme_minimal() + labs(title="EDSS during time by Sex", x = "Time", y = "Median EDSS")

# 3 ANDAMENTO PER NAIVE
ggplot(edss_long, aes(x = time, y = EDSS_val, color = naive, group = naive)) +
  stat_summary(fun = median, geom = "line", linewidth = 1.2, na.rm = TRUE) +
  stat_summary(fun = median, geom = "point", size = 3, na.rm = TRUE) +
  theme_minimal() + labs(title="EDSS during time by Naive", x = "Time", y = "Median EDSS")

# 4 ANDAMENTO PER SMOKING
ggplot(edss_long, aes(x = time, y = EDSS_val, color = smoking_habit, group = smoking_habit)) +
  stat_summary(fun = median, geom = "line", linewidth = 1.2, na.rm = TRUE) +
  stat_summary(fun = median, geom = "point", size = 3, na.rm = TRUE) +
  theme_minimal() + labs(title="EDSS during time by Smoking habit", x = "Time", y = "Median EDSS")

# 5 ANDAMENTO PER PHYSICAL ACTIVITY
#questo sembra avere un trend di aumento dell'EDSS per i pazienti sedentari rispetto agli attivi
ggplot(edss_long, aes(x = time, y = EDSS_val, color = physical_activity, group = physical_activity)) +
  stat_summary(fun = median, geom = "line", linewidth = 1.2, na.rm = TRUE) +
  stat_summary(fun = median, geom = "point", size = 3, na.rm = TRUE) +
  theme_minimal() + labs(title="EDSS during time by Physical Activity", x = "Time", y = "Median EDSS")

#analisi lesioni
#lesion burden
wilcox.test(as.numeric(as.character(EDSS_T0)) ~ lesion_burden, data = merged_data, exact = FALSE) # differenza non significativa
wilcox.test(as.numeric(as.character(EDSS_T6)) ~ lesion_burden, data = merged_data, exact = FALSE) # differenza non significativa
wilcox.test(as.numeric(as.character(EDSS_T12)) ~ lesion_burden, data = merged_data, exact = FALSE) # differenza non significativa
wilcox.test(as.numeric(as.character(EDSS_T24)) ~ lesion_burden, data = merged_data, exact = FALSE) # differenza non significativa

merged_data %>%
  filter(!is.na(EDSS_T0), !is.na(lesion_burden)) %>%
  mutate(EDSS_T0 = as.factor(EDSS_T0), lesion_burden = factor(lesion_burden, labels = c("No", "Yes"))) %>%
  ggplot(aes(x = lesion_burden, fill = EDSS_T0)) +
  geom_bar(position = "fill") + scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(option = "mako") + theme_minimal() + labs(title = "EDSS_T0 by lesion burden", x = "Lesion burden", y = "Percentage", fill = "EDSS")

#Bone marrow lesions
wilcox.test(as.numeric(as.character(EDSS_T0)) ~ bone_marrow_lesions, data = merged_data, exact = FALSE) # differenza non significativa
wilcox.test(as.numeric(as.character(EDSS_T6)) ~ bone_marrow_lesions, data = merged_data, exact = FALSE) # differenza non significativa
wilcox.test(as.numeric(as.character(EDSS_T12)) ~ bone_marrow_lesions, data = merged_data, exact = FALSE) # differenza non significativa
wilcox.test(as.numeric(as.character(EDSS_T24)) ~ bone_marrow_lesions, data = merged_data, exact = FALSE) # differenza non significativa

merged_data %>%
  filter(!is.na(EDSS_T0), !is.na(bone_marrow_lesions)) %>%
  mutate(EDSS_T0 = as.factor(EDSS_T0), bone_marrow_lesions = factor(bone_marrow_lesions, labels = c("No", "Yes"))) %>%
  ggplot(aes(x = bone_marrow_lesions, fill = EDSS_T0)) +
  geom_bar(position = "fill") + scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(option = "mako") + theme_minimal() + labs(title = "EDSS_T0 by bone marrow lesions", x = "Bone marrow lesions", y = "Percentage", fill = "EDSS")

#subtentorial lesions
wilcox.test(as.numeric(as.character(EDSS_T0)) ~ subtentorial_lesions, data = merged_data, exact = FALSE) # differenza non significativa
wilcox.test(as.numeric(as.character(EDSS_T6)) ~ subtentorial_lesions, data = merged_data, exact = FALSE) # differenza non significativa
wilcox.test(as.numeric(as.character(EDSS_T12)) ~ subtentorial_lesions, data = merged_data, exact = FALSE) # differenza non significativa
wilcox.test(as.numeric(as.character(EDSS_T24)) ~ subtentorial_lesions, data = merged_data, exact = FALSE) # differenza non significativa


merged_data %>%
  filter(!is.na(EDSS_T0), !is.na(subtentorial_lesions)) %>%
  mutate(EDSS_T0 = as.factor(EDSS_T0),subtentorial_lesions = factor(subtentorial_lesions, labels = c("No", "Yes"))) %>%
  ggplot(aes(x = subtentorial_lesions, fill = EDSS_T0)) +
  geom_bar(position = "fill") + scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(option = "mako") + theme_minimal() + labs(title = "EDSS_T0 by subtentorial_lesions", x = "Subtentorial_lesions", y = "Percentage", fill = "EDSS")

#gadolinium contrast
wilcox.test(as.numeric(as.character(EDSS_T0)) ~ gadolinium_contrast, data = merged_data, exact = FALSE) # differenza non significativa
wilcox.test(as.numeric(as.character(EDSS_T6)) ~ gadolinium_contrast, data = merged_data, exact = FALSE) # differenza non significativa
wilcox.test(as.numeric(as.character(EDSS_T12)) ~ gadolinium_contrast, data = merged_data, exact = FALSE) # differenza non significativa
wilcox.test(as.numeric(as.character(EDSS_T24)) ~ gadolinium_contrast, data = merged_data, exact = FALSE) # differenza non significativa

merged_data %>%
  filter(!is.na(EDSS_T0), !is.na(gadolinium_contrast)) %>%
  mutate( EDSS_T0 = as.factor(EDSS_T0), gadolinium_contrast = factor(gadolinium_contrast, labels = c("No", "Yes"))) %>%
  ggplot(aes(x = gadolinium_contrast, fill = EDSS_T0)) +
  geom_bar(position = "fill") + scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(option = "mako") + theme_minimal() + labs(title = "EDSS_T0 by gadolinium_contrast", x = "Gadolinium_contrast", y = "Percentage", fill = "EDSS")

##relazioni tra categoriche##
#relapse vs lesioni
chisq.test(table(merged_data$relapse_T6, merged_data$lesioni_T6))
chisq.test(table(merged_data$relapse_T18, merged_data$lesioni_T18))
chisq.test(table(merged_data$relapse_T24, merged_data$lesioni_T24)) #significativo
chisq.test(table(merged_data$relapse_T30, merged_data$lesioni_T30))

ggplot(merged_data, aes(x = relapse_T6, fill = lesioni_T6)) + geom_bar(position = "fill") + theme_minimal() + labs(y = "Proportion")
ggplot(merged_data, aes(x = relapse_T24, fill = lesioni_T24)) + geom_bar(position = "fill") + theme_minimal() + labs(y = "Proportion")

# relapse vs gadolinium #no correlzione
chisq.test(table(merged_data$relapse_T6, merged_data$gadolinium_contrast))
chisq.test(table(merged_data$relapse_T30, merged_data$gadolinium_contrast))

#naive vs relapse #no correlazione
chisq.test(table(merged_data$naive, merged_data$relapse_T6))
chisq.test(table(merged_data$naive, merged_data$relapse_T18))
chisq.test(table(merged_data$naive, merged_data$relapse_T24))


#glucocorticoid treatment before t0 e edss t0 #SIGNIFICATIVO
wilcox.test(as.numeric(as.character(EDSS_T0)) ~ gc_treatment, data = merged_data)

ggplot(merged_data, aes(x = gc_treatment, y = as.numeric(as.character(EDSS_T0)))) + geom_boxplot() + theme_minimal()






