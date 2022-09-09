library(tidyverse)
library(tidylog)
library(tableone)
library(MatchIt)
library(survival)
library(lubridate)
library(ggpubr)
library(SCCS)

load("C:/Users/LabPC14CSMPR/Desktop/Chris/HOPE/booster dose and COVID-19 mortality in multimorbidity/RDS/DX.RData")
# RX_latest <- read_rds("C:/Users/LabPC14CSMPR/Desktop/Chris/HOPE/booster dose and COVID-19 mortality in multimorbidity/RDS/RX_latest.RDS")
# anticoagulation <- RX_latest %>% filter(str_detect(item_cd, "APIX|DABI|EDOX|RIVA|WARF"))
# write_rds(anticoagulation, "anticoagulation.rds")

cohort <- read_rds("C:/Users/LabPC14CSMPR/Desktop/Chris/HOPE/booster dose and COVID-19 mortality in multimorbidity/RDS/4.cohort_full.RDS")
cohort <- cohort %>% filter(!is.na(patient_pssn)) %>% as_tibble()
cohort <- cohort %>% rename(vaccine.brand.1st = `Vaccine Brand.1st`, vaccine.brand.2nd = `Vaccine Brand.2nd`, 
                      vaccine.brand.3rd = `Vaccine Brand.3rd`, vaccine.brand.4th = `Vaccine Brand.4th`, 
                      date.of.vaccination.1st = `Date of vaccination.1st`, date.of.vaccination.2nd = `Date of vaccination.2nd`, 
                      date.of.vaccination.3rd = `Date of vaccination.3rd`, date.of.vaccination.4th = `Date of vaccination.4th`)

dx_clean
dx_latest

# AF
dx_af <- rbind(dx_clean, dx_latest) %>% as_tibble() %>% filter(str_detect(code, "^427\\.3")) %>% mutate(date = as_date(date))
dx_stroke <- rbind(dx_clean, dx_latest) %>% as_tibble() %>% filter(str_detect(code, "^433\\.[012389]1|^434|^436|^437\\.[01]")) %>% 
  filter(str_detect(Source, "1.IP")) %>% 
  mutate(date = as_date(date))
dx_se <- rbind(dx_clean, dx_latest) %>% as_tibble() %>% filter(str_detect(code, "^44[45]")) %>% 
  filter(str_detect(Source, "1.IP")) %>% 
  mutate(date = as_date(date))
dx_ICH <- rbind(dx_clean, dx_latest) %>% as_tibble() %>% filter(str_detect(code, "^43[012]")) %>% 
  filter(str_detect(Source, "1.IP")) %>% 
  mutate(date = as_date(date))
dx_GIB <- rbind(dx_clean, dx_latest) %>% as_tibble() %>% filter(str_detect(code, "^53[1234]\\.[0246]|^535\\.[01234567]1|^562\\.[01][23]|^569\\.3|^569\\.8[56]|^578\\.[019]")) %>% 
  filter(str_detect(Source, "1.IP")) %>% 
  mutate(date = as_date(date))
dx_OB <- rbind(dx_clean, dx_latest) %>% as_tibble() %>% filter(str_detect(code, "^423\\.0|^459\\.0|^593\\.81|^599\\.7|^623\\.8|^626\\.[26]|^719\\.1|^784\\.[78]|^786\\.3")) %>% 
  filter(str_detect(Source, "1.IP")) %>% 
  mutate(date = as_date(date))

dx_af_3weeksbefore <- dx_af %>% filter(date > date("2021-02-01") & date < date("2021-02-23"))

dx_pre <- dx_af %>%
  filter(date < date("2021-02-23")|is.na(date)) %>%
  group_by(patient_pssn) %>%
  slice(1)

dx_after <- 
  dx_stroke %>% rbind(dx_se) %>%
  # dx_ICH %>% rbind(dx_GIB) %>% rbind(dx_OB) %>%
  filter(!is.na(date)) %>%
  select(patient_pssn, date) %>%
  filter(date >= date("2021-02-23")) %>%
  group_by(patient_pssn) %>% # not recurrent event to be consider
  arrange(date) %>%
  slice(1) %>%
  ungroup()

hx_stroke_se <- rbind(dx_clean, dx_latest) %>% 
  filter(str_detect(code, "^433\\.[012389]1|^434|^436|^437\\.[01]|^44[45]")) %>% 
  as_tibble() %>% 
  mutate(date = as_date(date)) %>%
  filter(date < date("2021-02-23")|is.na(date))

hx_bleeding <- rbind(dx_clean, dx_latest) %>% 
  filter(str_detect(code, "^43[012]|^53[1234]\\.[0246]|^535\\.[01234567]1|^562\\.[01][23]|^569\\.3|^569\\.8[56]|^578\\.[019]|^423\\.0|^459\\.0|^593\\.81|^599\\.7|^623\\.8|^626\\.[26]|^719\\.1|^784\\.[78]|^786\\.3")) %>% 
  as_tibble() %>% 
  mutate(date = as_date(date)) %>%
  filter(date < date("2021-02-23")|is.na(date))

# # sen5 eTable7
# dx_after <- dx_latest %>%
#   as_tibble() %>% filter(str_detect(code, "^8[0-2][0-9]")) %>% mutate(date = as_date(date)) %>%
#   filter(!is.na(date)) %>%
#   select(patient_pssn, date) %>%
#   filter(date >= date("2021-02-23")) %>%
#   group_by(patient_pssn) %>% # not recurrent event to be consider
#   arrange(date) %>%
#   slice(1) %>%
#   ungroup()

# flow chart
# cohort %>% filter(patient_pssn %in% dx_pre$patient_pssn) %>% .$vaccine.brand.1st %>% table
# cohort %>% filter(patient_pssn %in% dx_pre$patient_pssn) %>% filter(vaccine.brand.1st == "BioNTech/Fosun") %>% filter(is.na(vaccine.brand.2nd))
# cohort %>% filter(patient_pssn %in% dx_pre$patient_pssn) %>% filter(vaccine.brand.1st == "BioNTech/Fosun") %>% filter(!is.na(vaccine.brand.2nd)) %>% filter(is.na(vaccine.brand.3rd))
# cohort %>% filter(patient_pssn %in% dx_pre$patient_pssn) %>% filter(vaccine.brand.1st == "BioNTech/Fosun") %>% filter(!is.na(vaccine.brand.2nd))  %>% .$vaccine.brand.3rd %>% table
# 
# cohort %>% filter(patient_pssn %in% dx_pre$patient_pssn) %>% filter(vaccine.brand.1st == "Sinovac") %>% filter(is.na(vaccine.brand.2nd))
# cohort %>% filter(patient_pssn %in% dx_pre$patient_pssn) %>% filter(vaccine.brand.1st == "Sinovac") %>% filter(!is.na(vaccine.brand.2nd)) %>% filter(is.na(vaccine.brand.3rd))
# cohort %>% filter(patient_pssn %in% dx_pre$patient_pssn) %>% filter(vaccine.brand.1st == "Sinovac") %>% filter(!is.na(vaccine.brand.2nd))  %>% .$vaccine.brand.3rd %>% table
# 
# cohort %>% filter(patient_pssn %in% dx_pre$patient_pssn) %>% filter(patient_pssn %in% dx_after$patient_pssn) %>% .$vaccine.brand.1st %>% table
# cohort %>% filter(patient_pssn %in% dx_pre$patient_pssn) %>% filter(patient_pssn %in% dx_after$patient_pssn) %>% filter(vaccine.brand.1st == "BioNTech/Fosun") %>% filter(is.na(vaccine.brand.2nd))
# cohort %>% filter(patient_pssn %in% dx_pre$patient_pssn) %>% filter(patient_pssn %in% dx_after$patient_pssn) %>% filter(vaccine.brand.1st == "BioNTech/Fosun") %>% filter(!is.na(vaccine.brand.2nd)) %>% filter(is.na(vaccine.brand.3rd))
# cohort %>% filter(patient_pssn %in% dx_pre$patient_pssn) %>% filter(patient_pssn %in% dx_after$patient_pssn) %>% filter(vaccine.brand.1st == "BioNTech/Fosun") %>% filter(!is.na(vaccine.brand.2nd))  %>% .$vaccine.brand.3rd %>% table
# 
# cohort %>% filter(patient_pssn %in% dx_pre$patient_pssn) %>% filter(patient_pssn %in% dx_after$patient_pssn) %>% filter(vaccine.brand.1st == "Sinovac") %>% filter(is.na(vaccine.brand.2nd))
# cohort %>% filter(patient_pssn %in% dx_pre$patient_pssn) %>% filter(patient_pssn %in% dx_after$patient_pssn) %>% filter(vaccine.brand.1st == "Sinovac") %>% filter(!is.na(vaccine.brand.2nd)) %>% filter(is.na(vaccine.brand.3rd))
# cohort %>% filter(patient_pssn %in% dx_pre$patient_pssn) %>% filter(patient_pssn %in% dx_after$patient_pssn) %>% filter(vaccine.brand.1st == "Sinovac") %>% filter(!is.na(vaccine.brand.2nd))  %>% .$vaccine.brand.3rd %>% table

LAB_ALL_COVID <- readRDS("C:/Users/LabPC14CSMPR/Desktop/Chris/HOPE/booster dose and COVID-19 mortality in multimorbidity/RDS/LAB_ALL_COVID.RDS")
hx_covid <- LAB_ALL_COVID[grepl("^21[3567]$",T_NUM) & result=="detected"] %>% as_tibble() %>% 
  group_by(patient_pssn) %>% arrange(date) %>% slice(1) %>% ungroup()

df2 <- cohort %>% 
  # filter(!patient_pssn %in% dx_hf_3weeksbefore$patient_pssn) %>% # sen4 eTable6
  mutate(death_date_ymd = as_date(death_date_ymd)) %>% 
  filter(patient_pssn %in% dx_pre$patient_pssn) %>%
  filter(patient_pssn %in% dx_after$patient_pssn) %>%
  # filter(!patient_pssn %in% hx_stroke_se$patient_pssn) %>% # should be changed when changing outcome
  # filter(!patient_pssn %in% hx_bleeding$patient_pssn) %>% # should be changed when changing outcome
  mutate(eventdate = dx_after$date[match(patient_pssn, dx_after$patient_pssn)]) %>% 
  # filter(is.na(death_date_ymd)) %>% 
  # filter(!is.na(date.of.vaccination.1st)) %>%
  group_by(patient_pssn) %>% 
  slice(1) %>% 
  ungroup() %>% 
  filter((vaccine.brand.1st == "BioNTech/Fosun" & (vaccine.brand.3rd == "BioNTech/Fosun" | is.na(vaccine.brand.3rd))) | (vaccine.brand.1st == "Sinovac" & (vaccine.brand.3rd == "Sinovac" | is.na(vaccine.brand.3rd))) | is.na(vaccine.brand.1st)) %>% 
  mutate(type = if_else(vaccine.brand.1st == "BioNTech/Fosun", "BNT162b2", "")) %>% 
  mutate(type = if_else(vaccine.brand.1st == "Sinovac", "CoronaVac", type)) %>% 
  mutate(type = if_else(is.na(vaccine.brand.1st), "unvaccinated", type)) %>% 
  mutate(date.covid = hx_covid$date[match(patient_pssn, hx_covid$patient_pssn)]) %>% 
  mutate(hx.covid = if_else(is.na(date.covid)|is.na(date.of.vaccination.1st)|date.covid > date.of.vaccination.1st, 0, 1))

# final_cohort <- df2 %>% filter((vaccine.brand.1st == "BioNTech/Fosun" & (vaccine.brand.3rd == "BioNTech/Fosun" | is.na(vaccine.brand.3rd))) | (vaccine.brand.1st == "Sinovac" & (vaccine.brand.3rd == "Sinovac" | is.na(vaccine.brand.3rd))) | is.na(vaccine.brand.1st))
# write_rds(final_cohort, "HF result/final_cohort.rds")

df2 %>% filter(death_date_ymd != "") %>% filter(type == "BNT162b2") %>% .$death_diag_cd %>% table()
df2 %>% filter(death_date_ymd != "") %>% filter(type == "CoronaVac") %>% .$death_diag_cd %>% table()

df2 %>% filter(type == "BNT162b2") %>% .$hx.covid %>% table()
df2 %>% filter(type == "CoronaVac") %>% .$hx.covid %>% table()

df3 <- df2 %>% 
  # filter(patient_pssn %in% filter(anticoagulation, str_detect(item_cd, "APIX|DABI|EDOX|RIVA"))$patient_pssn) %>%
  # filter(patient_pssn %in% filter(anticoagulation, str_detect(item_cd, "WARF"))$patient_pssn) %>%
  # filter(Age >= 80) %>%
  # filter(Age < 80) %>%
  # filter(sex == "M") %>%
  filter(sex == "F") %>%
  # filter(hx.covid == 0) %>% # sen2
  mutate(nid = seq_along(patient_pssn)) %>%
  filter(vaccine.brand.1st == "BioNTech/Fosun" & (vaccine.brand.3rd == "BioNTech/Fosun" | is.na(vaccine.brand.3rd)) | is.na(vaccine.brand.1st)) %>% # main analysis: include BNT and unvaccinated
  # filter(vaccine.brand.1st == "BioNTech/Fosun" & (vaccine.brand.3rd == "BioNTech/Fosun" | is.na(vaccine.brand.3rd))) %>%
  # filter(vaccine.brand.1st == "Sinovac" & (vaccine.brand.3rd == "Sinovac" | is.na(vaccine.brand.3rd)) | is.na(vaccine.brand.1st)) %>% # main analysis: include sinovac and unvaccinated
  # filter(vaccine.brand.1st == "Sinovac" & (vaccine.brand.3rd == "Sinovac" | is.na(vaccine.brand.3rd))) %>%
  # mutate(dob = as_date(paste0(dob_y, "-01-01"))) %>% 
  mutate(dob = as_date("2021-01-01")) %>% 
  mutate(eventdate = as.integer(as_date(eventdate) - dob)) %>% 
  mutate(obs_start = as.integer(as_date("2021-02-23")-dob)) %>% 
  mutate(vaccdate1 = as.integer(as_date(date.of.vaccination.1st)-dob)) %>% 
  mutate(vaccdate2 = as.integer(as_date(date.of.vaccination.2nd)-dob)) %>% 
  mutate(vaccdate3 = as.integer(as_date(date.of.vaccination.3rd)-dob)) %>% 
  mutate(vaccdate1_p14 = as.integer(vaccdate1 + 14)) %>% 
  mutate(vaccdate2_p14 = as.integer(vaccdate2 + 14)) %>% 
  mutate(vaccdate3_p14 = as.integer(vaccdate3 + 14)) %>% 
  mutate(obs_end = as.integer(pmin(as_date("2022-03-31"), as_date(death_date_ymd), na.rm = T)-dob)) %>% 
  select(nid, patient_pssn, starts_with("obs"), starts_with("vaccdate"), starts_with("eventdate"), vaccine.brand.1st)

# stroke_pts <- df3 %>% filter(eventdate >= vaccdate1 & eventdate <= vaccdate1 + 27)
# cohort %>% filter(patient_pssn %in% stroke_pts$patient_pssn) %>% View()
# rbind(dx_clean, dx_latest) %>% filter(patient_pssn %in% stroke_pts$patient_pssn) %>% arrange(patient_pssn, date)
# stroke_pts %>% filter(patient_pssn %in% anticoagulation$patient_pssn)

# df3 <- df3 %>% mutate(eventdate = if_else(nid == 3174, 357, as.double(eventdate))) %>% mutate(eventdate = if_else(nid == 3174, 357, as.double(eventdate)))
# df3 %>% filter(nid %in% c(1295,1916,2916,6888))
# df3 %>% filter(!is.na(vaccdate3)) %>% filter(eventdate > 320)
# df3 %>% filter(nid %in% c(1334,7007,8430))
# df3 %>% filter(!is.na(vaccdate3)) %>% filter(eventdate > 320)

# ageq <- floor(quantile(df3$eventdate, seq(0.1,0.9,0.1),
#                        names=F, na.rm = T))

# ageq <- cumsum(c(90, 61, 61, 61, 61, 61))
ageq <- cumsum(c(59, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31, 28))

(eventde_result <- eventdepenexp(indiv = nid,
                                 astart = obs_start,
                                 aend = obs_end,
                                 aevent = eventdate,
                                 adrug = cbind(vaccdate1, vaccdate1_p14, vaccdate2, vaccdate2_p14, vaccdate3, vaccdate3_p14),
                                 aedrug = cbind(vaccdate1+13, vaccdate1_p14 +13, vaccdate2+13, vaccdate2_p14+13, vaccdate3+13, vaccdate3_p14+13),
                                 sameexpopar = F,
                                 agegrp=ageq,
                                 dataformat = "multi",
                                 data = df3))

# # sen2 eTable4
# (eventde_result <- standardsccs(event ~ vaccdate1+age,
#                              indiv = nid,
#                              astart = obs_start,
#                              aend = obs_end,
#                              aevent = eventdate,
#                              adrug = cbind(vaccdate1 - 14, vaccdate1, vaccdate1_p14, vaccdate2, vaccdate2_p14, vaccdate3, vaccdate3_p14),
#                              aedrug = cbind(vaccdate1 - 1, vaccdate1+13, vaccdate1_p14 +13, vaccdate2+13, vaccdate2_p14+13, vaccdate3+13, vaccdate3_p14+13),
#                              sameexpopar = F,
#                              agegrp=ageq,
#                              dataformat = "multi",
#                              data = df3))

str_c(round(eventde_result$conf.int[seq(1,6),], 2)[,1], "(",
      round(eventde_result$conf.int[seq(1,6),], 2)[,3], "-",
      round(eventde_result$conf.int[seq(1,6),], 2)[,4], ")")


mace_data <- formatdata(indiv = nid,
                        astart = obs_start,
                        aend = obs_end,
                        aevent = eventdate,
                        adrug = cbind(vaccdate1, vaccdate1_p14, vaccdate2, vaccdate2_p14, vaccdate3, vaccdate3_p14),
                        aedrug = cbind(vaccdate1+13, vaccdate1_p14 +13, vaccdate2+13, vaccdate2_p14+13, vaccdate3+13, vaccdate3_p14+13),
                        sameexpopar = F,
                        agegrp=ageq,
                        dataformat = "multi",
                        data = df3)
tapply(mace_data$event, mace_data$vaccdate1, sum)
tapply(mace_data$interval, mace_data$vaccdate1, sum)
mace_data %>% filter(vaccdate1 == 1)
mace_data %>% filter(vaccdate1 == 2)


result <- rbind(cbind(eventde_result$conf.int, eventde_result$coefficients)[, c(1,3,4,9)][seq(1,6),]) %>% 
  as_tibble() %>% 
  mutate_if(is.numeric, round, digits = 2) %>% 
  mutate_all(format, nsmall = 2)
colnames(result) <- c("exp", "lower", "upper", "p")
result <- result %>%
  mutate(IRR = str_c(exp, " (", lower, "-", upper, ")", sep = "")) %>% 
  mutate(p = if_else(p == "0.00", "<.01", p)) %>% 
  select(IRR, p)
result <- rbind(rep("", 2), result[seq(1,6),])
result <- result %>% mutate(No.events = c(tapply(mace_data$event, mace_data$vaccdate1, sum))) %>% 
  mutate(follow.up = c(tapply(mace_data$interval, mace_data$vaccdate1, sum))) %>% 
  mutate(No.events = as.integer(No.events)) %>% 
  mutate(follow.up = as.integer(follow.up)) %>% 
  mutate(absolute.rate = format(round(No.events*1000/follow.up, digits = 1), nsmall = 1)) %>% 
  select(No.events, follow.up, absolute.rate, IRR, p)
result <- rbind(rep("", 5), result[1,], rep("", 5), result[seq(2,3),], rep("", 5), result[seq(4,5),], rep("", 5), result[seq(6,7),])
# write_csv(result, "AF result/stroke biontech M sen2.csv")
# write_csv(result, "AF result/stroke sinovac M sen2.csv")
# write_csv(result, "AF result/bleeding biontech F sen2.csv")
# write_csv(result, "AF result/bleeding sinovac F sen2.csv")

df3$vaccine.brand.1st %>% table()
df3 %>% filter(is.na(vaccine.brand.1st)) %>% dim()

