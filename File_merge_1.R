source("Libraries.R")
metadata_pazienti_lesioni <- read_excel("Input/metadati_pazienti_Alessandro.xlsx")
metadata_paper <- read.csv("Input/metadata_paper.csv")
dati_filtrati_pazienti <- metadata_pazienti_lesioni[startsWith(metadata_pazienti_lesioni$ID, "MAV"), ]
# metto le date come date anno-mese-giorno
str(dati_filtrati_pazienti)
# queste sono le date in chr
dati_filtrati_pazienti$data_prelievo_T0 <- as.numeric(dati_filtrati_pazienti$data_prelievo_T0)
dati_filtrati_pazienti$data_prelievo_T0 <- as.Date(dati_filtrati_pazienti$data_prelievo_T0, origin = "1899-12-30")
dati_filtrati_pazienti$data_referto_T12 <- as.numeric(dati_filtrati_pazienti$data_referto_T12)
dati_filtrati_pazienti$data_referto_T12 <- as.Date(dati_filtrati_pazienti$data_referto_T12, origin = "1899-12-30")
dati_filtrati_pazienti$data_referto_T24 <- as.numeric(dati_filtrati_pazienti$data_referto_T24)
dati_filtrati_pazienti$data_referto_T24 <- as.Date(dati_filtrati_pazienti$data_referto_T24, origin = "1899-12-30")
dati_filtrati_pazienti$data_cambio_terapia <- as.numeric(dati_filtrati_pazienti$data_cambio_terapia)
dati_filtrati_pazienti$data_cambio_terapia <- as.Date(dati_filtrati_pazienti$data_cambio_terapia, origin = "1899-12-30")

# queste sono le date in posixct
dati_filtrati_pazienti$`data di nascita` <- as.Date(dati_filtrati_pazienti$`data di nascita`)
dati_filtrati_pazienti$data_referto_T6 <- as.Date(dati_filtrati_pazienti$data_referto_T6)
dati_filtrati_pazienti$data_referto_T18 <- as.Date(dati_filtrati_pazienti$data_referto_T18)
dati_filtrati_pazienti$data_referto_T30 <- as.Date(dati_filtrati_pazienti$data_referto_T30)
dati_filtrati_pazienti$data_referto_T36 <- as.Date(dati_filtrati_pazienti$data_referto_T36)
# date numeriche
dati_filtrati_pazienti$`DATA ESORDIO` <- as.Date(dati_filtrati_pazienti$`DATA ESORDIO`, origin = "1899-12-30")
dati_filtrati_pazienti$data_diagnosi <- as.Date(dati_filtrati_pazienti$data_diagnosi, origin = "1899-12-30")

# modifico le variabili non corrette
str(dati_filtrati_pazienti)

dati_filtrati_pazienti$EDSS_T24 <- as.numeric(dati_filtrati_pazienti$EDSS_T24)
dati_filtrati_pazienti$EDSS_T12 <- as.numeric(dati_filtrati_pazienti$EDSS_T12)

cols <- c("lesioni_T6", "relapse_T6", "lesioni_T12", "relapse_T12", "lesioni T18", "lesioni T24", "relapse T24", "lesioni T30", "relapse T30", "lesioni T36", "relapse T36", "relapse T18", "cambio_terapia")
dati_filtrati_pazienti[cols] <- lapply(dati_filtrati_pazienti[cols], factor, levels = c(0, 1))


# trasformo ? in NA
# sostituire "?" con NA in tutte le colonne
dati_filtrati_pazienti$EDSS_T24 <- as.numeric(
  ifelse(
    as.character(dati_filtrati_pazienti$EDSS_T24) == "?", # se è "?"
    NA, # sostituisci con NA
    dati_filtrati_pazienti$EDSS_T24 # altrimenti tieni il valore
  )
)
summary(dati_filtrati_pazienti$EDSS_T24)

dati_filtrati_pazienti$EDSS_T12 <- as.numeric(
  ifelse(
    as.character(dati_filtrati_pazienti$EDSS_T12) == "?", # se è "?"
    NA, # sostituisci con NA
    dati_filtrati_pazienti$EDSS_T12 # altrimenti tieni il valore
  )
)

# cambio nome alle colonne
names(dati_filtrati_pazienti)[names(dati_filtrati_pazienti) == "data di nascita"] <- "Data_di_nascita"
names(dati_filtrati_pazienti)[names(dati_filtrati_pazienti) == "DATA ESORDIO"] <- "Data_esordio"
names(dati_filtrati_pazienti)[names(dati_filtrati_pazienti) == "data_diagnosi"] <- "Data_diagnosi"
names(dati_filtrati_pazienti)[names(dati_filtrati_pazienti) == "data_prelievo_T0"] <- "Data_prelievo_T0"
# sostituisco tutti gli spazi dei nomi con underscore
names(dati_filtrati_pazienti) <- gsub(" ", "_", names(dati_filtrati_pazienti))
names(dati_filtrati_pazienti)[names(dati_filtrati_pazienti) == "nuova_TERAPIA"] <- "nuova_terapia"

## DT: questo codice va bene. In alternativa puoi usare la funzione clean_names() del pacchetto janitor, che rende i nomi delle colonne minuscole e separate da underscore

#-----------------------------------------
# modifico metadata_paper
names(metadata_paper)[names(metadata_paper) == "id"] <- "ID"
names(metadata_paper)[names(metadata_paper) == "disease_actvity"] <- "disease_activity"
# estraggo solo MAV
metadata_MAV_paper <- metadata_paper[startsWith(metadata_paper$ID, "MAV"), ]
str(metadata_MAV_paper)
# trasformo simplified_batch in livelli corretti
metadata_MAV_paper$simplified_batch[metadata_MAV_paper$simplified_batch == "3-Jan"] <- "1-3"
metadata_MAV_paper$simplified_batch[metadata_MAV_paper$simplified_batch == "5-Apr"] <- "4-5"
metadata_MAV_paper$simplified_batch <- factor(metadata_MAV_paper$simplified_batch, levels = c("1-3", "4-5"))
metadata_MAV_paper$sex <- factor(metadata_MAV_paper$sex, levels = c("M", "F"))
metadata_MAV_paper$physical_activity <- factor(metadata_MAV_paper$physical_activity, levels = c("active", "sedentary", "partially_active"))
metadata_MAV_paper$disease_activity <- factor(metadata_MAV_paper$disease_activity, levels = c("higher_activity", "lower_activity"))

## DT: il codice va bene, ma si può scrivere in maniera più compatta usando tidyverse. ti faccio un esempio sotto:
# metadata_MAV_paper <- metadata_MAV_paper %>%
#   mutate(
#     across(c(age, pyr_mds, sequencing_batch, days_gc_collection),
#            as.numeric),

#     across(c(naive, smoking_habit, antibiotic_use, spike_in),
#            ~ factor(.x,
#                     levels = c("yes", "no"),
#                     labels = c("Yes", "No"))),

#     across(c(glucocorticoid_treatment, gc_treatment, prognosis_at_onset),
#            ~ factor(.x,
#                     levels = c("positive", "negative"),
#                     labels = c("Positive", "Negative"))),

#     across(c(lesion_burden, bone_marrow_lesions, subtentorial_lesions,
#              gadolinium_contrast, full_recovery,
#              new_lesions_t12, relapse_t12,
#              new_lesions_t24, relapse_t24,
#              therapy_change),
#            ~ factor(.x,
#                     levels = c("0", "1"),
#                     labels = c("No", "Yes")))
#   )
# qui ho aggiunto anche levels e labels. I levels sono i livelli che vuoi dare al fattore, i labels sono i nomi che vuoi dare ai livelli.
# ad esempio, se hai una colonna con "1" e "0", puoi trasformarla in un fattore con labels "Yes" e "No" in modo che appaia così nei grafici e nelle tabelle.


#uso codice suggerito da Dora
metadata_MAV_paper <- metadata_MAV_paper %>%
  mutate(
    across(c(age, sequencing_batch, days_gc_collection),
       as.numeric),
    across(c(naive, smoking_habit, antibiotic_use, spike_in),
     ~ factor(.x,
       levels = c("yes", "no"),
       labels = c("Yes", "No"))),
    across(c(glucocorticoid_treatment, gc_treatment, prognosis_at_onset),
         ~ factor(.x,
           levels = c("positive", "negative"),
           labels = c("Positive", "Negative"))),
    across(c(lesion_burden, bone_marrow_lesions, subtentorial_lesions,
        gadolinium_contrast, full_recovery,
        new_lesions_t12, relapse_t12,
        new_lesions_t24, relapse_t24,
        therapy_change),
       ~ factor(.x,
            levels = c("0", "1"),
            labels = c("No", "Yes")))
 )




str(metadata_MAV_paper)

# cambio date con aa/m/gg
metadata_MAV_paper$sample_collection_date <- as.Date(metadata_MAV_paper$sample_collection_date, format = "%m/%d/%y")

# elimino colonne doppie
metadata_MAV_paper <- metadata_MAV_paper %>%
  select(-c(
    X, category, edss, edss_t0, edss_t12, edss_t24,
    onset_date, diagnosis_date,
    new_lesions_t12, relapse_t12,
    new_lesions_t24, relapse_t24,
    medas_percent, therapy_change, pyr_mds,center
  ))

names(metadata_MAV_paper)

# cambio nome ID in modo che coincidano
metadata_MAV_paper$ID <- sub("T0$", "", metadata_MAV_paper$ID)


# merge dei dataset
merge_metadata <- merge(dati_filtrati_pazienti, metadata_MAV_paper, by = "ID", all.x = TRUE)
names(merge_metadata)
str(merge_metadata)

binary_cols_01 <- c(
  "lesioni_T6", "relapse_T6", "lesioni_T12", "relapse_T12",
  "lesioni_T18", "lesioni_T24", "relapse_T24",
  "lesioni_T30", "relapse_T30",
  "lesioni_T36", "relapse_T36",
  "relapse_T18", "cambio_terapia"
)

merge_metadata[binary_cols_01] <- lapply(
  merge_metadata[binary_cols_01],
  function(x) factor(as.character(x),
                     levels = c("0", "1"),
                     labels = c("No", "Yes"))
)

# cambio glc_date
merge_metadata$glc_date_before_t0 <- as.Date(merge_metadata$glc_date_before_t0, format = "%m/%d/%y")

# modifico MAV07 come non naive
merge_metadata[merge_metadata$ID == "MAV07", "naive"] <- "No"
merge_metadata[merge_metadata$ID == "MAV27", "naive"] <- "No"

# modifico la colonna therapy completandola poichè ho solo MAV
merge_metadata$therapy <- "cladribine"

# completo la colonna sex
merge_metadata$ID[is.na(merge_metadata$sex)] # per sapere chi sono gli NA

merge_metadata[merge_metadata$ID %in% c("MAV07", "MAV12", "MAV15", "MAV16", "MAV20", "MAV23", "MAV24", "MAV25", "MAV31", "MAV34"), "sex"] <- "F"
merge_metadata[merge_metadata$ID %in% c("MAV08", "MAV22"), "sex"] <- "M"
merge_metadata[, c("ID", "sex")] # per controllare

# completo la colonna age (at t0)
merge_metadata$ID[is.na(merge_metadata$age)]
merge_metadata[merge_metadata$ID == "MAV07", "age"] <- 42
merge_metadata[merge_metadata$ID == "MAV08", "age"] <- 36
merge_metadata[merge_metadata$ID == "MAV12", "age"] <- 36
merge_metadata[merge_metadata$ID == "MAV15", "age"] <- 48
merge_metadata[merge_metadata$ID == "MAV16", "age"] <- 55
merge_metadata[merge_metadata$ID == "MAV20", "age"] <- 49
merge_metadata[merge_metadata$ID == "MAV22", "age"] <- 39
merge_metadata[merge_metadata$ID == "MAV23", "age"] <- 52
merge_metadata[merge_metadata$ID == "MAV24", "age"] <- 35
merge_metadata[merge_metadata$ID == "MAV25", "age"] <- 42

# cambio nome alla colonna age in age_at_t0
names(merge_metadata)[names(merge_metadata) == "age"] <- "age_at_t0"
merge_metadata$age_at_t0 <- as.numeric(merge_metadata$age_at_t0)

# completo la colonna naive
merge_metadata$ID[is.na(merge_metadata$naive)]
merge_metadata[merge_metadata$ID %in% c("MAV08", "MAV12", "MAV15", "MAV20", "MAV24"), "naive"] <- "No"
merge_metadata[merge_metadata$ID %in% c("MAV16", "MAV22", "MAV23", "MAV25"), "naive"] <- "Yes"


# metto la data_t0 di MAV22 che ho trovato
merge_metadata[merge_metadata$ID == "MAV22", "Data_prelievo_T0"] <- as.Date("2021-11-01")

# modifico il medas in modo che sia numerico (I tre trasformati in NA sono recuperabili dalla tabella originale, es. contattare madre)
merge_metadata$MEDAS[merge_metadata$ID %in% c("MAV10", "MAV19", "MAV21")] <- NA
merge_metadata$MEDAS <- as.numeric(merge_metadata$MEDAS)


# AGGIUNGO DATI TROVATI CON ALESSANDRO IN DATA 02/04
merge_metadata[merge_metadata$ID == "MAV16", c("bmi", "smoking_habit", "physical_activity")] <- list(20.1, "Yes", "partially_active")
merge_metadata[merge_metadata$ID == "MAV12", c("bmi", "smoking_habit", "physical_activity")] <- list(20.3, "No", "partially_active")
merge_metadata[merge_metadata$ID == "MAV15", c("bmi", "smoking_habit", "physical_activity")] <- list(32.0, "Yes", "sedentary")
merge_metadata[merge_metadata$ID == "MAV25", c("bmi", "smoking_habit", "physical_activity")] <- list(25.7, "No", "sedentary")
merge_metadata[merge_metadata$ID == "MAV22", "MEDAS"] <- 42
merge_metadata[merge_metadata$ID == "MAV10", "MEDAS"] <- 35
merge_metadata[merge_metadata$ID == "MAV21", "MEDAS"] <- 35
merge_metadata[merge_metadata$ID == "MAV23", c("bmi", "smoking_habit", "physical_activity")] <- list(20.43, "No", "active")


# MAV DI ALTRI CENTRI
merge_metadata[merge_metadata$ID == "MAV31", c("Data_prelievo_T0", "cambio_terapia", "EDSS_T0", "EDSS_T12", "EDSS_T24", "Data_esordio", "Data_diagnosi", "gc_treatment", "smoking_habit")] <- list(as.Date("2021-07-01"), "No", 2.5, 2.5, 2.5, as.Date("2017-08-07"), as.Date("2017-12-15"), "Negative", "No")
merge_metadata[merge_metadata$ID == "MAV31", c("lesion_burden", "bone_marrow_lesions", "subtentorial_lesions", "gadolinium_contrast", "clinical_presentation", "full_recovery", "lesioni_T12", "relapse_T12", "lesioni_T24", "relapse_T24", "Data_di_nascita", "naive", "previous_therapy", "age_at_t0")] <- list("Yes", "Yes", "Yes", "No", "sensoriale", "No", "Yes", "No", "No", "No", as.Date("1975-10-07"), "No", "Rebif", 46)

merge_metadata[merge_metadata$ID == "MAV32", c("Data_prelievo_T0", "cambio_terapia", "EDSS_T0", "EDSS_T12", "EDSS_T24", "Data_esordio", "Data_diagnosi", "gc_treatment", "smoking_habit")] <- list(as.Date("2021-08-23"), "No", 4, 4, 4, as.Date("2015-07-01"), as.Date("2016-07-01"), "Negative", "Yes")
merge_metadata[merge_metadata$ID == "MAV32", c("lesion_burden", "bone_marrow_lesions", "subtentorial_lesions", "gadolinium_contrast", "clinical_presentation", "full_recovery", "lesioni_T12", "relapse_T12", "lesioni_T24", "relapse_T24", "Data_di_nascita", "naive", "previous_therapy", "age_at_t0")] <- list("Yes", "Yes", "No", "No", "tronco", "Yes", "No", "No", "No", "No", as.Date("1976-01-26"), "No", "Lemtrada", 45)

merge_metadata[merge_metadata$ID == "MAV33", c("Data_prelievo_T0", "cambio_terapia", "EDSS_T0", "EDSS_T12", "EDSS_T24", "Data_esordio", "Data_diagnosi", "gc_treatment", "smoking_habit")] <- list(as.Date("2021-09-06"), "No", 2.5, 2.5, 2.5, as.Date("2021-02-15"), as.Date("2021-08-20"), "Negative", "No")
merge_metadata[merge_metadata$ID == "MAV33", c("lesion_burden", "bone_marrow_lesions", "subtentorial_lesions", "gadolinium_contrast", "clinical_presentation", "full_recovery", "lesioni_T12", "relapse_T12", "lesioni_T24", "relapse_T24", "Data_di_nascita", "naive", "age_at_t0")] <- list("Yes", "Yes", "No", "No", "sensoriale", "Yes", "No", "No", "No", "No", as.Date("1998-11-09"), "Yes", 22)

merge_metadata[merge_metadata$ID == "MAV34", c("Data_prelievo_T0", "EDSS_T0", "EDSS_T12", "EDSS_T24", "Data_esordio", "Data_diagnosi", "gc_treatment")] <- list(as.Date("2022-04-19"), 2, 2, 2, as.Date("2008-01-01"), as.Date("2019-11-01"), "Positive")
merge_metadata[merge_metadata$ID == "MAV34", c("lesion_burden", "bone_marrow_lesions", "gadolinium_contrast", "clinical_presentation", "full_recovery", "lesioni_T12", "relapse_T12", "lesioni_T24", "relapse_T24", "Data_di_nascita", "age_at_t0")] <- list("No", "Yes", "No", "norb/sensoriale", "Yes", "No", "No", "No", "No", as.Date("1987-09-22"), 34)


# creo la colonna age calcolata al giorno d'oggi
merge_metadata$age <- floor(time_length(interval(as.Date(merge_metadata$Data_di_nascita), Sys.Date()), "years"))


# salvo il merge come csv
write.csv(merge_metadata, "C:/Users/gialu/Desktop/tesi_magistrale/Git.hub/Thesis_crepaldi/Input/merged_data.csv", row.names = FALSE)

str(merge_metadata)
