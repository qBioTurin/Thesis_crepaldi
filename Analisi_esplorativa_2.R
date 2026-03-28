#ANALISI METADATA INIZIALE
merged_data = read.csv("C:/Users/gialu/Desktop/tesi_magistrale/Git.hub/Thesis_crepaldi/Input/merged_data.csv")
str(merged_data)

#metto a posto la tipologia di variabile di nuovo:
#le date
cols_date <- c("Data_di_nascita", "Data_esordio", "Data_diagnosi", "Data_prelievo_T0", "data_referto_T6", "data_referto_T12", "data_referto_T18", "data_referto_T24", "data_referto_T30", "data_referto_T36", "data_cambio_terapia", "sample_collection_date", "glc_date_before_t0")
merged_data[cols_date] <- lapply(merged_data[cols_date], as.Date)

#categorico
cols_factor <- c("lesioni_T6", "relapse_T6", "lesioni_T12", "relapse_T12", "lesioni_T18","relapse_T18" ,"lesioni_T24","relapse_T24", "lesioni_T30", "relapse_T30", "lesioni_T36", "relapse_T36","cambio_terapia", "lesion_burden", "bone_marrow_lesions", "subtentorial_lesions", "gadolinium_contrast", "full_recovery")
merged_data[cols_factor] = lapply(merged_data[cols_factor], factor, levels = c(0, 1))
cols_factor2 <- c("naive", "smoking_habit", "antibiotic_use", "spike_in")
merged_data[cols_factor2] = lapply(merged_data[cols_factor2], factor, levels = c("yes", "no"))
cols_factor3 <- c("glucocorticoid_treatment", "gc_treatment", "prognosis_at_onset")
merged_data[cols_factor3] = lapply(merged_data[cols_factor3], factor, levels = c("positive", "negative"))
merged_data$physical_activity= factor(merged_data$physical_activity, levels = c( "active", "sedentary", "partially_active"))
merged_data$simplified_batch <- factor(merged_data$simplified_batch, levels = c("1-3", "4-5"))
merged_data$disease_activity= factor(merged_data$disease_activity, levels = c("higher_activity", "lower_activity"))
merged_data$sex = factor(merged_data$sex, levels=c("M", "F"))

#numerici
cols_numeric = c("MEDAS", "EDSS_T0", "EDSS_T6", "cortisone_T6", "EDSS_T12", "cortisone_T12", "EDSS_T18", "cortisone_T18", "EDSS_T24", "cortisone_T24", "EDSS_T30", "cortisone_T30", "EDSS_T36", "cortisone_T36", "age_at_t0", "bmi", "pyr_mds", "sequencing_batch", "dna_quantification", "days_gc_collection", "percent_th17", "treg_il10_pos", "treg_cd39_pos", "age", "age_at_diagnosis")
merged_data[cols_numeric] = lapply(merged_data[cols_numeric], as.numeric)

#analisi esplorativa min,max, media, mediana di alcuni dati
summary(merged_data[sapply(merged_data, is.numeric)])

lapply(merged_data[sapply(merged_data, is.factor)],table, useNA = "always")


#analisi delle possibili correlazioni tra variabili:

#EDSS vs smoking_habit
t.test(EDSS_T0 ~ smoking_habit, data = merged_data) #differenza non significativa
t.test(EDSS_T6 ~ smoking_habit, data = merged_data) #differenza non significativa
t.test(EDSS_T12 ~ smoking_habit, data = merged_data) #differenza non significativa
t.test(EDSS_T18 ~ smoking_habit, data = merged_data) #differenza non significativa
t.test(EDSS_T24 ~ smoking_habit, data = merged_data) #differenza non significativa
t.test(EDSS_T30 ~ smoking_habit, data = merged_data) #differenza non significativa
t.test(EDSS_T36 ~ smoking_habit, data = merged_data) #non ci sono abbastanza osservazioni

#EDSS vs age_at_diagnosis
cor.test(merged_data$age_at_diagnosis, merged_data$EDSS_T0) #correlazione non significativa
cor.test(merged_data$age_at_diagnosis, merged_data$EDSS_T6) #correlazione non significativa
cor.test(merged_data$age_at_diagnosis, merged_data$EDSS_T12) #correlazione non significativa
cor.test(merged_data$age_at_diagnosis, merged_data$EDSS_T18) #correlazione non significativa
cor.test(merged_data$age_at_diagnosis, merged_data$EDSS_T24) #correlazione positiva SIGNIFICATIVA
cor.test(merged_data$age_at_diagnosis, merged_data$EDSS_T30) #correlazione non significativa
cor.test(merged_data$age_at_diagnosis, merged_data$EDSS_T36) #correlazione non significativa

#EDSS vs bmi
cor.test(merged_data$bmi, merged_data$EDSS_T0) #correlazione non significativa
cor.test(merged_data$bmi, merged_data$EDSS_T6) #correlazione non significativa
cor.test(merged_data$bmi, merged_data$EDSS_T12) #correlazione non significativa
cor.test(merged_data$bmi, merged_data$EDSS_T18) #correlazione non significativa
cor.test(merged_data$bmi, merged_data$EDSS_T24) #correlazione non significativa
cor.test(merged_data$bmi, merged_data$EDSS_T30) #correlazione non significativa
cor.test(merged_data$bmi, merged_data$EDSS_T36) #correlazione non significativa

#EDSS vs MEDAS
cor.test(merged_data$MEDAS, merged_data$EDSS_T0) #correlazione non significativa
cor.test(merged_data$MEDAS, merged_data$EDSS_T6) #correlazione non significativa
cor.test(merged_data$MEDAS, merged_data$EDSS_T12) #correlazione non significativa
cor.test(merged_data$MEDAS, merged_data$EDSS_T18) #correlazione non significativa
cor.test(merged_data$MEDAS, merged_data$EDSS_T24) #correlazione non significativa
cor.test(merged_data$MEDAS, merged_data$EDSS_T30) #correlazione non significativa
cor.test(merged_data$MEDAS, merged_data$EDSS_T36) #non ci sono abbastanza osservazioni

#age_at_diagnosis vs relapse
t.test(age_at_diagnosis ~ relapse_T6, data = merged_data) #correlazione non significativa
t.test(age_at_diagnosis ~ relapse_T24, data = merged_data) #correlazione non significativa
t.test(age_at_diagnosis ~ relapse_T30, data = merged_data) #significativo


