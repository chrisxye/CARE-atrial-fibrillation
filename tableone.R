library(tidyverse)
library(tidylog)
library(tableone)
library(MatchIt)
library(survival)
library(lubridate)
library(ggpubr)
library(SCCS)
library(data.table)
library(readxl)

load("C:/Users/LabPC14CSMPR/Desktop/Chris/HOPE/booster dose and COVID-19 mortality in multimorbidity/RDS/DX.RData")
final_cohort <- read_rds("AF result/final_cohort_stroke.rds") %>% select(patient_pssn, vaccine.brand.1st, Age, sex, death_date_ymd, death_diag_cd) %>%
  mutate(index.date = date("2021-02-23"))
# final_cohort <- read_rds("AF result/final_cohort_bleeding.rds") %>% select(patient_pssn, vaccine.brand.1st, Age, sex, death_date_ymd, death_diag_cd) %>%
#   mutate(index.date = date("2021-02-23"))

adj <- c("dx.htn", "dx.cancer", "dx.crf", "dx.respdz", "dx.dm", "dx.dementia", "rx.ras", "rx.bb", "rx.ccb", "rx.diuretic", "rx.nitrates", "rx.lipid", "rx.insulin", "rx.dm", "rx.oac", "rx.apt")

covar_dz <- function(ct, name, icd) {
  hx1820 <- dx_clean[grepl(icd, code, ignore.case=T), unique(patient_pssn)]
  hx21 <- dx_latest_index[date < index.date & grepl(icd, code, ignore.case=T), unique(patient_pssn)]
  hx <- unique(c(hx1820, hx21))
  ct[, c(name) := as.numeric(patient_pssn %in% hx)]
  cat(length(hx),"\n")
  return(ct)
}

drg <- function(name, dc, db) (!is.na(drug_codes[Name==name,Item_Code]) & grepl(drug_codes[Name==name,Item_Code], dc, ignore.case=T)) | (!is.na(drug_codes[Name==name,BNF]) & grepl(drug_codes[Name==name,BNF], db, ignore.case=T))

covar_drug <- function(ct, name) {  # 365 days before index
  hx21 <- drug_latest_index[presc_start < index.date & presc_end >= (index.date-365) & drg(name, item_cd, bnfno_p), unique(patient_pssn)]
  hx <- unique(c(hx21))
  ct[, c(name) := as.numeric(patient_pssn %in% hx)]
  cat(length(hx),"\n")
  return(ct)
}


final_cohort <- setDT(final_cohort)
dx_codes <- setDT(read_excel("AF result/codes.xlsx", sheet="codes"))[!is.na(Name)] %>% filter(Name %in% adj)
dx_latest_index <- merge(dx_latest, final_cohort[, .(patient_pssn, index.date)], by="patient_pssn")
for(i in 1:nrow(dx_codes)) {
  cat(dx_codes[i,Name], "...")
  final_cohort <- covar_dz(final_cohort, dx_codes[i,Name], dx_codes[i,Regex])
}

drug_latest <- read_rds("C:/Users/LabPC14CSMPR/Desktop/Chris/HOPE/booster dose and COVID-19 mortality in multimorbidity/RDS/RX_latest.RDS")
drug_latest[, presc_start := as_date(presc_start_date_ymd)]
drug_latest[, presc_end := presc_start + presc_duration_day]
drug_latest_index <- merge(drug_latest, final_cohort[, .(patient_pssn, index.date)], by="patient_pssn")
drug_latest_index[, index.date := as.Date(index.date)]
drug_codes <- setDT(read_excel("AF result/codes.xlsx", sheet="drugs"))[!is.na(Name)] %>% filter(Name %in% adj)
for(i in 1:nrow(drug_codes)) {
  cat(drug_codes[i,Name], "...")
  final_cohort <- covar_drug(final_cohort, drug_codes[i,Name])
}

final_cohort <- final_cohort %>% 
  mutate(type = if_else(vaccine.brand.1st == "BioNTech/Fosun", "BNT162b2", "")) %>% 
  mutate(type = if_else(vaccine.brand.1st == "Sinovac", "CoronaVac", type)) %>% 
  mutate(type = if_else(is.na(vaccine.brand.1st), "unvaccinated", type))

final_cohort[,c(seq(8, 23))] <- lapply(final_cohort[,c(seq(8, 23))], factor)

vars <- c("Age", "sex", colnames(select(final_cohort, starts_with("dx"), starts_with("rx"))))
tableone <- CreateTableOne(vars = vars, strata = c("type"), data = final_cohort)
print(tableone, quote = F, noSpaces = F)
# write.csv(print(tableone, quote = F, noSpaces = F), file = "AF result/tableone_SCCS_stroke.csv")
# write.csv(print(tableone, quote = F, noSpaces = F), file = "AF result/tableone_SCCS_bleeding.csv")
