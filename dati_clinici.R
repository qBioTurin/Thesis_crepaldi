source("Libraries.R")
merged <- read.csv2("C:\\Users\\gialu\\Desktop\\merged.csv")

merged$age_at_T0 <- floor(
  as.numeric(interval(merged$data_nascita, merged$data_prelievo_T0) / years(1))
)



