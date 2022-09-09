library(tidyverse)
library(data.table)
library(readxl)
library(fasttime)
library(lubridate)

# -1. Helpers ----

drg <- function(name, dc, db) (!is.na(drug_codes[Name==name,Item_Code]) & grepl(drug_codes[Name==name,Item_Code], dc, ignore.case=T)) | (!is.na(drug_codes[Name==name,BNF]) & grepl(drug_codes[Name==name,BNF], db, ignore.case=T))

covar_dz <- function(ct, name, icd) {
  hx1820 <- dx_clean_[grepl(icd, code, ignore.case=T), unique(patient_pssn)]
  hx21 <- dx_latest_index[date < index.date & grepl(icd, code, ignore.case=T), unique(patient_pssn)]
  hx <- unique(c(hx1820, hx21))
  ct[, c(name) := as.numeric(patient_pssn %in% hx)]
  cat(length(hx),"\n")
  return(ct)
}

bl_dz <- function(ct, dz_def) {
  for(i in 1:nrow(dz_def)) {
    cat(dz_def[i,Name], "...")
    ct <- covar_dz(ct, dz_def[i,Name], dz_def[i,Regex])
    gc()
  }
  return(ct)
}

covar_drug <- function(ct, name, within=90) {
  hx21 <- drug_latest_index[presc_start < index.date & presc_end >= (index.date-within) & drg(name, item_cd, bnfno_p), unique(patient_pssn)]
  hx <- unique(c(hx21))
  coln <- if(within==90) name else paste0(name, ".", within)
  ct[, c(coln) := as.numeric(patient_pssn %in% hx)]
  cat(length(hx),"\n")
  return(ct)
}

gp.match <- function(cohort, case.var, match.var, seed) {
  setorderv(cohort, c(case.var, match.var, "PseudoID"))
  cases <- cohort[get(case.var)==1]
  ctrls <- cohort[get(case.var)==0]
  gps <- unique(cases[, match.var, with=F])
  .Random.seed <- seed
  
  res <- list()
  for(i in 1:nrow(gps)) {
    cat("\n", gps[i,Age], gps[i,sex])
    gp1 <- cases[Age==gps[i,Age] & sex==gps[i,sex]]
    gp0 <- ctrls[Age==gps[i,Age] & sex==gps[i,sex]]
    K = floor(nrow(gp0) / nrow(gp1)); cat(" -> 1 :",K)
    gp0_rand <- sample(gp0$PseudoID, nrow(gp0), replace=F)
    
    mat <- matrix(if(K==0) numeric() else gp0_rand[1:(nrow(gp1)*K)], nrow=nrow(gp1), ncol=K)
    mat <- cbind(mat, c(if(K==0) gp0_rand else gp0_rand[-c(1:(nrow(gp1)*K))], rep(NA, nrow(gp1)*(K+1)-nrow(gp0))))
    gpm <- cbind(gp1[, .(PseudoID)], mat)
    res[[i]] <- gpm
  }
  
  return(res)
}

output_matches <- function(md, ct) {
  md_sizes <- unlist(lapply(md, nrow))
  ms <- rbindlist(lapply(1:length(md), function(i) melt(cbind(match=(sum(md_sizes[0:(i-1)]) + (1:nrow(md[[i]]))), md[[i]][, c("PseudoID", grep("^V[0-9]+", names(md[[i]]), value=T)), with=F]), id.vars=c("match"), value.name="PseudoID")))
  ct.m <- merge(ct, ms[, .(PseudoID, match, match_order=as.character(variable))], by="PseudoID")
  return(ct.m)
}

pp.match <- function(dz_pop, nodz_pop, K, seed) {
  dz_cases <- dz_pop[case==1]
  nodz_cases <- nodz_pop[case==1]
  nodz_cases_matched <- pp.match_(dz_cases, nodz_cases, c("Age","sex"), K, seed)
  cat("\n", nrow(dz_cases), ":", nrow(nodz_cases), "|", length(nodz_cases_matched), "-", nrow(nodz_pop[case==0]))
  rbind(nodz_cases[PseudoID %in% nodz_cases_matched], nodz_pop[case==0])
}

pp.match_ <- function(ct_ref, ct_tar, match.var, K, seed) {
  setorderv(ct_ref, c(match.var, "PseudoID"))
  setorderv(ct_tar, c(match.var, "PseudoID"))
  gps <- unique(ct_ref[, match.var, with=F])
  info <- merge(ct_ref[, .N, keyby=match.var], ct_tar[, .N, keyby=match.var], by=match.var, all=T)[, K := N.y/N.x]; print(summary(info$K))
  
  res <- character()
  .Random.seed <- seed
  
  for(i in 1:nrow(gps)) {
    cat("\n", gps[i,Age], gps[i,sex])
    gp_ref <- ct_ref[Reduce(`&`,lapply(match.var, function(v) ct_ref[[v]]==gps[i][[v]]))]
    gp_tar <- ct_tar[Reduce(`&`,lapply(match.var, function(v) ct_tar[[v]]==gps[i][[v]]))]
    gp_tar_rand <- sample(gp_tar$PseudoID, pmin(K*nrow(gp_ref), nrow(gp_tar)), replace=F)
    res <- c(res, gp_tar_rand)
    cat("->", nrow(gp_ref), ":", length(gp_tar_rand), "|", nrow(gp_tar))
  }
  
  return(res)
}

# 0. Helpers ----

table1 <- function(cohort, strata) {
  cohort_ <- copy(cohort)
  cohort_[, (grep("^(dx|px|ip|rx|hx)",names(cohort_),value=T)) := lapply(.SD, as.logical), .SDcols=grep("^(dx|px|ip|rx|hx)",names(cohort_),value=T)]
  
  tb1_def <- setDT(read_excel("codes.xlsx", sheet="table1"))
  t1 <- CreateTableOne(tb1_def[!is.na(Name), Name], c(strata), data=cohort_)
  print(t1, smd=T, test=F)
  return(t1)
}

as.data.table.TableOne <- function(t) {
  tb1_def <- setDT(read_excel("codes.xlsx", sheet="table1"))[!is.na(Name), .(Name, Description)]
  tb1_def <- rbind(list(Name="n", Description="n"), tb1_def)
  t <- as.data.frame(print(t, test=F, dropEqual=T, noSpaces=T))
  varlabels <- rownames(t)
  t$Name = sub("^([a-zA-Z0-9._]+).*$", "\\1", varlabels)
  t <- merge(as.data.table(t), tb1_def, by="Name", all.x=T, sort=F)
  t$Description = paste0(t$Description, sub("^([a-zA-Z0-9._]+)", "", varlabels))
  t[!is.na(Description), Name:=Description][, Description:=NULL]
  return(t)
}

as.data.table.clogit <- function(r, n=3) {
  t <- cbind(OR=exp(r$coefficients), exp(confint(r)))
  t <- as.data.table(t, keep.rownames=T)
  #t$rn <- gsub("vaccinated.brand", "Vaccinated: ", t$rn)
  t[, `OR (95% CI)` := paste0(format(round(OR,n),nsmall=n,trim=T), " (", format(round(`2.5 %`,n),nsmall=n,trim=T), " - ", format(round(pmin(`97.5 %`,9999),n),nsmall=n,trim=T), ")")]
  return(t[, .(` `=rn, `OR (95% CI)`)])
}



# 0. Load data ----

# cohort
cohort <- read_rds("C:/Users/LabPC14CSMPR/Desktop/Chris/HOPE/booster dose and COVID-19 mortality in multimorbidity/RDS/4.cohort_full.RDS")
cohort <- cohort %>% filter(!is.na(patient_pssn)) %>% as_tibble()
cohort <- cohort %>% rename(vaccine.brand.1st = `Vaccine Brand.1st`, vaccine.brand.2nd = `Vaccine Brand.2nd`, 
                            vaccine.brand.3rd = `Vaccine Brand.3rd`, vaccine.brand.4th = `Vaccine Brand.4th`, 
                            date.of.vaccination.1st = `Date of vaccination.1st`, date.of.vaccination.2nd = `Date of vaccination.2nd`, 
                            date.of.vaccination.3rd = `Date of vaccination.3rd`, date.of.vaccination.4th = `Date of vaccination.4th`)

# DX
load("C:/Users/LabPC14CSMPR/Desktop/Chris/HOPE/booster dose and COVID-19 mortality in multimorbidity/RDS/DX.RData")
# AF
dx_af <- rbind(dx_clean, dx_latest) %>% as_tibble() %>% filter(str_detect(code, "^427\\.3")) %>% mutate(date = as_date(date))
dx_stroke <- rbind(dx_clean, dx_latest) %>% as_tibble() %>% filter(str_detect(code, "^433\\.[012389]1|^434|^436|^437\\.[01]")) %>% 
  filter(str_detect(Source, "1.IP")) %>% 
  mutate(date = as_date(date))
dx_se <- rbind(dx_clean, dx_latest) %>% as_tibble() %>% filter(str_detect(code, "^44[45]")) %>% 
  filter(str_detect(Source, "1.IP")) %>% 
  mutate(date = as_date(date))
# dx_ICH <- rbind(dx_clean, dx_latest) %>% as_tibble() %>% filter(str_detect(code, "^43[012]")) %>% 
#   filter(str_detect(Source, "1.IP")) %>% 
#   mutate(date = as_date(date))
# dx_GIB <- rbind(dx_clean, dx_latest) %>% as_tibble() %>% filter(str_detect(code, "^53[1234]\\.[0246]|^535\\.[01234567]1|^562\\.[01][23]|^569\\.3|^569\\.8[56]|^578\\.[019]")) %>% 
#   filter(str_detect(Source, "1.IP")) %>% 
#   mutate(date = as_date(date))
# dx_OB <- rbind(dx_clean, dx_latest) %>% as_tibble() %>% filter(str_detect(code, "^423\\.0|^459\\.0|^593\\.81|^599\\.7|^623\\.8|^626\\.[26]|^719\\.1|^784\\.[78]|^786\\.3")) %>% 
#   filter(str_detect(Source, "1.IP")) %>% 
#   mutate(date = as_date(date))

# # Lab COVID
# lab_all_covid <- readRDS("C:/Users/LabPC14CSMPR/Desktop/Chris/HOPE/booster dose and COVID-19 mortality in multimorbidity/RDS/LAB_ALL_COVID.RDS")

# # Attendance latest
# load("C:/Users/LabPC14CSMPR/Desktop/Chris/HOPE/booster dose and COVID-19 mortality in multimorbidity/RDS/A10X.RData")
# admin_latest[, date := as.Date(fastPOSIXct(date))]

cohort_AF <- cohort %>% filter(patient_pssn %in% filter(dx_af, date < date("2021-02-23")|is.na(date))$patient_pssn)
dx_case <- rbind(dx_stroke, dx_se) %>% 
  filter(date >= date("2021-02-23")) %>%
  filter(patient_pssn %in% cohort_AF$patient_pssn) %>% 
  group_by(patient_pssn) %>% # not recurrent event to be consider
  arrange(date) %>%
  slice(1) %>%
  ungroup()

case <- cohort_AF %>% 
  filter(patient_pssn %in% dx_case$patient_pssn) %>% 
  filter(!patient_pssn %in% filter(rbind(dx_stroke, dx_se), date < date("2021-02-23"))$patient_pssn) %>% 
  mutate(stroke.date = dx_case$date[match(patient_pssn, dx_case$patient_pssn)])

control <- cohort_AF %>% 
  filter(!patient_pssn %in% filter(rbind(dx_stroke, dx_se), date < date("2021-02-23"))$patient_pssn) %>% 
  filter(!patient_pssn %in% case$patient_pssn) %>% 
  mutate(stroke.date = NA)

cohort_cc <- rbind(mutate(case, case = 1), mutate(control, case = 0))
cohort_cc$case %>% table()
cohort_cc <- as.data.table(cohort_cc)
# Pseudo index date
dist.match <- merge(cohort_cc[case==1, .N, keyby=.(Age, sex)], cohort_cc[case==0, .N, keyby=.(Age, sex)], by=c("Age","sex"), all.x=T)[, Ratio:=N.y/N.x]
set.seed(475)
seed <- .Random.seed
matched <- gp.match(cohort_cc, "case", c("Age","sex"), seed)
cohort.matched <- output_matches(matched, cohort_cc)
cohort.matched$case %>% table()

case_ind <- cohort.matched[match_order=="PseudoID", .(index.date=stroke.date), keyby=match]
cohort.matched <- merge(cohort.matched, case_ind, by="match")
cohort.matched[, index.date := as.Date(index.date)]

cohort_cc$case %>% table()
# 4. Baseline characteristics ----
cohort_cc <- data.table(cohort.matched)
cohort_cc[, Age_5yb := floor(Age/5)]

# * History of comobidities (2018 onwards) ---
dx_clean_ <- dx_clean[patient_pssn %in% cohort_cc$patient_pssn]
dx_latest_index <- merge(dx_latest, cohort_cc[, .(patient_pssn, index.date)], by="patient_pssn")

dxpx_codes <- setDT(read_excel("codes.xlsx", sheet="codes"))[!is.na(Name)]
cohort_cc <- bl_dz(cohort_cc, dxpx_codes)
cohort_cc[, score.cci := (dx.mi+dx.chf+dx.pvd+dx.cbd+dx.copd+dx.dementia+dx.paralysis+(dx.dm_com0&!dx.dm_com1)+dx.dm_com1*2+dx.crf*2+(dx.liver_mild&!dx.liver_modsev)+dx.liver_modsev*3+dx.ulcers+dx.ra+dx.aids*6+dx.cancer*2+dx.cancer_mets*6)]
cohort_cc[score.cci==0, score.cci.trunc := "0"]
cohort_cc[score.cci==1 | score.cci==2, score.cci.trunc := "1-2"]
cohort_cc[score.cci==3 | score.cci==4, score.cci.trunc := "3-4"]
cohort_cc[score.cci>=5, score.cci.trunc := ">=5"]

# * Medication use within 90 days (unless specified) ---
# RX
gc()
drug_latest <- readRDS("C:/Users/LabPC14CSMPR/Desktop/Chris/HOPE/booster dose and COVID-19 mortality in multimorbidity/RDS/RX_latest.RDS")
drug_latest[, presc_start := as.Date(fastPOSIXct(presc_start_date_ymd))]
drug_latest[, presc_end := presc_start + presc_duration_day]
gc()
drug_latest_index <- merge(drug_latest, cohort_cc[, .(patient_pssn, index.date)], by="patient_pssn")
gc()
drug_codes <- setDT(read_excel("codes.xlsx", sheet="drugs"))[!is.na(Name)]
for(i in 1:nrow(drug_codes)) {
  cat(drug_codes[i,Name], "...")
  cohort_cc <- covar_drug(cohort_cc, drug_codes[i,Name])
  gc()
}
rm(i)

# 6. Analysis ----
cohort.matched.bl <- cohort_cc %>% as.data.table() %>% mutate(date.of.vaccination.1st=ymd_hms(date.of.vaccination.1st)) %>% 
  mutate(date.of.vaccination.1st = as_date(date.of.vaccination.1st))

# count only on/before event
cohort.matched.bl[date.of.vaccination.1st > index.date, `:=`(date.of.vaccination.1st=NA, vaccine.brand.1st=NA)]

#cohort.matched.bl[`Vaccine Brand.1st`!=`Vaccine Brand.2nd`, .N]  # =0
cohort.matched.bl[, vacc.1x_biontech := as.numeric(!is.na(date.of.vaccination.1st) & vaccine.brand.1st=="BioNTech/Fosun")]
cohort.matched.bl[, vacc.1x_sinovac := as.numeric(!is.na(date.of.vaccination.1st) & vaccine.brand.1st=="Sinovac")]


# factor vaccination status
table(cohort.matched.bl[, rowSums(.SD), .SDcols=grep("^vacc[.]", names(cohort.matched.bl), value=T)], useNA="always")
cohort.matched.bl$id <- seq_along(cohort.matched.bl$match)
cohort.matched.bl <- merge(cohort.matched.bl, unique(melt(cohort.matched.bl[, c("id", "match", grep("^vacc[.]", names(cohort.matched.bl), value=T)), with=F], c("id","match"))[value==1, .(id, match, vacc_status=sub("vacc.","",variable))]), by=c("id","match"), all.x=T)
table(cohort.matched.bl$vacc_status, useNA="always")

# Duration since nth dose of vaccine
cohort.matched.bl[, time.since.1st_dose := as.numeric(index.date-date.of.vaccination.1st)]
# cohort.matched.bl[, time.since.recent_dose := as.numeric(index.date - pmax(date.of.vaccination.1st, na.rm=T))]
# cohort.matched.bl[, .(min=min(time.since.recent_dose), max=max(time.since.recent_dose)), keyby=vacc_status]

# Days since study start
cohort.matched.bl[, calendar_days := as.numeric(index.date - min(index.date))]



# 2. Descriptives ----
library(tableone)

get_tb1 <- function(ct) {
  ts <- list(table1(ct, c("case")), table1(ct, c("vacc_status")), table1(ct[case==1], c("vacc_status")), table1(ct[case==0], c("vacc_status")))
  names(ts) <- paste0(c("case-ctrl","vacc_status","case-vacc","ctrl-vacc"), " (", deparse(substitute(ct)), ")")
  return(sapply(ts, as.data.table, simplify=F))
}

get_tb2 <- function(ct) {
  ts <- list(ct[, .(Case=sum(case==1), Control=sum(case==0)), keyby=vacc_status])
  names(ts) <- paste0(c("events"), " (", deparse(substitute(ct)), ")")
  return(ts)
}



# 3. Conditional logit ----
library(survival)
#data <- copy(cohort.matched.bl)

#f0 <- "case ~ vacc.1x_biontech + vacc.1x_sinovac + vacc.2x_biontech + vacc.2x_sinovac + vacc.3x_biontech + vacc.3x_sinovac + vacc.3x_bbs + vacc.3x_ssb"
f0_ <- c("1x_biontech", "1x_sinovac")

adj <- c("dx.cancer", "dx.crf", "dx.respdz", "dx.dm", "dx.dementia", "rx.ras", "rx.bb", "rx.ccb", "rx.diuretic", "rx.nitrates", "rx.lipid", "rx.insulin", "rx.dm", "rx.oac", "rx.apt")

clean_adj <- function(adj, data) {
  data.table(var=adj, N=unlist(lapply(adj, function(a) data[get(a)==1, .N])))[N>0, var]
}


get_res <- function(f,a,d) {
  res <- list(clogit(as.formula(paste(f, "+ strata(match)")), d),
              clogit(as.formula(paste(f, "+", paste(a, collapse=" + "), "+ strata(match)")), d))
  names(res) <- c(paste0(deparse(substitute(f)), " (", deparse(substitute(d)), ")"),
                  paste0(deparse(substitute(f)), ".", deparse(substitute(a)), " (", deparse(substitute(d)), ")"))
  res.tab <- sapply(res, as.data.table, simplify=F)
  return(list(model=res, table=res.tab))
}

#r0 <- get_res(f0, adj, data)
#r1 <- get_res(f0, adj1, data)

library(writexl)
#write_xlsx(c(get_tb1(data), get_tb2(data), r0$table, r1$table), paste0("res.", OUTCOME, SA_NAME, ".xlsx"))

#save(r0, r1, file=paste0("res.", OUTCOME, SA_NAME, ".RData"))



# Timing since vaccination

# if 1-dose-only
data <- copy(cohort.matched.bl[(vacc.1x_biontech==1|vacc.1x_sinovac==1) | is.na(vacc_status)])
# data <- data %>% filter(sex == "F")
data <- data %>% filter(sex == "M")

data$case %>% table()

evt.timing <- list()
res.timing <- list()

for(p in list(c(0,13), c(14,27))) {
  cat(paste0(p[1],"-",p[2]), "...\n")
  dt <- copy(data[is.na(vacc_status) | (time.since.1st_dose >= p[1] & time.since.1st_dose <= p[2])])
  dt <- dt[match %in% dt[case==1, match]]
  dt <- dt[!match %in% dt[, sum(case==0), by=match][V1==0, match]]
  ft <- paste("case ~", paste(paste0("vacc.", intersect(f0_, dt[,.N, keyby=vacc_status][N>0, vacc_status])), collapse=" + "))
  adj_ <- clean_adj(adj, dt)
  r0 <- get_res(ft, adj_, dt)
  res.timing[[paste0(p[1],"-",p[2])]] <- list(model=c(r0$model), table=c(r0$table))
  evt.timing[[paste0(p[1],"-",p[2])]] <- get_tb2(dt)
}
rm(p, dt, ft, r0)

merge_ <- function(L, by=NULL) Reduce(function(x,y) merge(x, y, all=T, sort=F), L)
res.timing.summary <- c(
  list(events = merge_(lapply(names(res.timing), function(n) setnames(evt.timing[[n]]$`events (dt)`[c(NA,f0_), .(vacc_status, paste(Case,"/",Control))], c("vacc_status", paste(n,"days")))), by="vacc_status")),
  lapply(list("crude"="ft (dt)", "adj"="ft.adj_ (dt)"), function(x) merge_(lapply(names(res.timing), function(n) setnames(res.timing[[n]]$table[[x]], c("vacc_status", paste(n,"days")))), by="vacc_status")) )

tab.timing <- sapply(list("tab.adj"="adj"), function(a) {
  t <- lapply(list("1x_biontech", "1x_sinovac"), 
              function(n) merge_(list(data.table::transpose(res.timing.summary$events, keep.names="Days", make.names="vacc_status")[, .(Days=sub(" days","",Days), Case=tstrsplit(get(n)," / ")[[1]], Control=tstrsplit(get(n)," / ")[[2]])], 
                                      data.table::transpose(res.timing.summary$crude, keep.names="Days", make.names="vacc_status")[, .(Days=sub(" days","",Days), `Crude OR (95% CI)`=gsub(" - ","-",get(paste0("vacc.",n))))], 
                                      data.table::transpose(res.timing.summary[[a]], keep.names="Days", make.names="vacc_status")[, .(Days=sub(" days","",Days), `Adjusted OR (95% CI)`=gsub(" - ","-",get(paste0("vacc.",n))))])))
  rbind(c("BNT162b2", rep(list(""),ncol(t[[1]])-1)), t[[1]], c("CoronaVac", rep(list(""),ncol(t[[2]])-1)), t[[2]], fill=T)
}, USE.NAMES = T, simplify = F)

# write_xlsx(c(res.timing.summary, tab.timing), paste0("res.timing.", OUTCOME, SA_NAME, ".xlsx"))
