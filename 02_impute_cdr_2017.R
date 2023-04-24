# 02_impute_cdr_2017.R ----

library(here)
library(rgdal)
library(tidyverse)
library(lubridate)

source(here("code", "mat_utils.R"))

load(here("data","cleaned_2017_CDR","all-transitions-daily.Rdata"))
bang_union <- readOGR(dsn = here("data","bgd_admbnda_adm4_bbs_20180410"), layer = "bgd_admbnda_adm4_bbs_20180410", stringsAsFactors = F) 
bang_union_df = as(bang_union, "data.frame")

# add upazilas to align with 2020 data
trans.sub <- trans %>%
  # change these to align wit shape file to add district
  mutate(UnionFrom = case_when(substr(UnionFrom, 1, 4) == "3039" ~ paste0("4539", substr(UnionFrom, 5, 8)),
                               substr(UnionFrom, 1, 4) == "3061" ~ paste0("4561", substr(UnionFrom, 5, 8)),
                               substr(UnionFrom, 1, 4) == "3072" ~ paste0("4572", substr(UnionFrom, 5, 8)),
                               substr(UnionFrom, 1, 4) == "3089" ~ paste0("4589", substr(UnionFrom, 5, 8)),
                               TRUE ~ as.character(UnionFrom))) %>%
  mutate(UnionTo = case_when(substr(UnionTo, 1, 4) == "3039" ~ paste0("4539", substr(UnionTo, 5, 8)),
                             substr(UnionTo, 1, 4) == "3061" ~ paste0("4561", substr(UnionTo, 5, 8)),
                             substr(UnionTo, 1, 4) == "3072" ~ paste0("4572", substr(UnionTo, 5, 8)),
                             substr(UnionTo, 1, 4) == "3089" ~ paste0("4589", substr(UnionTo, 5, 8)),
                             TRUE ~ as.character(UnionTo))) %>%
  mutate(UnionFrom = as.character(UnionFrom), UnionTo = as.character(UnionTo)) %>%
  left_join(.,bang_union_df %>% select(ADM4_PCODE, ADM3_PCODE, ADM2_EN), by = c("UnionTo" = "ADM4_PCODE")) %>%
  rename(UpaTo = ADM3_PCODE,
         DistrictTo = ADM2_EN)  %>%
  left_join(.,bang_union_df %>% select(ADM4_PCODE, ADM3_PCODE, ADM2_EN), by = c("UnionFrom" = "ADM4_PCODE")) %>%
  rename(UpaFrom = ADM3_PCODE,
         DistrictFrom = ADM2_EN)

trans.sub %>% filter(is.na(UpaTo)) %>% distinct(UnionTo) %>% pull(UnionTo)
trans.sub %>% filter(is.na(UpaFrom)) %>% distinct(UnionFrom) %>% pull(UnionFrom)

# Remove two unions "20157867" "30333449"
trans.sub <- trans.sub %>% filter(!is.na(UpaTo) & !is.na(UpaFrom))

#trans.sub = readRDS(here("data","cleaned_2017_CDR", "trans.sub.rds"))

# Group by upazila ----

# proportion trips from Dhaka to other districts
trans.sub.upa = trans.sub %>%
  select(!average) %>%
  group_by(UpaTo, UpaFrom) %>%
  summarise(across(starts_with("Day"), ~ sum(.x, na.rm=T))) %>%
  ungroup()

#trans.sub.upa = readRDS(file = here("data","cleaned_2017_CDR", "trans.sub.upa.rds"))

# Make long df ----
distinct_sets = trans.sub.upa %>% distinct(UpaTo, UpaFrom) %>% nrow()

district_upa <- read.csv(file = here("data","popcount_by_union.csv"), header=T) %>%
  mutate(upazila_code = as.character(Upa_Code), Dis_Code = as.character(Dis_Code), District = as.character(District)) %>%
  dplyr::select(upazila_code, Upazila, Dis_Code, District) %>% distinct()

district_upa$District[district_upa$District == "Cox'S Bazar"] <- "Cox's Bazar"

# updated this to remove the double entry for upa code "303334" in the unique Upazila name
district_upa <- district_upa %>% filter(Upazila != "Kaliganj Paurashava") # two entries for this upa_code "303334"

trans.sub.upa.long = trans.sub.upa %>% pivot_longer(!c(1:2), names_to = "day", values_to = "value") %>%
  mutate(date = rep(seq(as.Date("2017-04-01"), as.Date("2017-09-30"), 1), times =distinct_sets)) %>%
  mutate(weekend = ifelse(wday(as.Date(date)) == 1, 1, 0) | ifelse(wday(as.Date(date)) == 7, 1, 0)) %>%
  # change back upazila codes to align with pop data and district crosswalk
  mutate(UpaFrom = case_when(substr(UpaFrom, 1, 4) == "4539" ~ paste0("3039", substr(UpaFrom, 5, 7)),
                             substr(UpaFrom, 1, 4) == "4561" ~ paste0("3061", substr(UpaFrom, 5, 7)),
                             substr(UpaFrom, 1, 4) == "4572" ~ paste0("3072", substr(UpaFrom, 5, 7)),
                             substr(UpaFrom, 1, 4) == "4589" ~ paste0("3089", substr(UpaFrom, 5, 7)),
                             TRUE ~ as.character(UpaFrom))) %>%
  mutate(UpaTo = case_when(substr(UpaTo, 1, 4) == "4539" ~ paste0("3039", substr(UpaTo, 5, 7)),
                           substr(UpaTo, 1, 4) == "4561" ~ paste0("3061", substr(UpaTo, 5, 7)),
                           substr(UpaTo, 1, 4) == "4572" ~ paste0("3072", substr(UpaTo, 5, 7)),
                           substr(UpaTo, 1, 4) == "4589" ~ paste0("3089", substr(UpaTo, 5, 7)),
                           TRUE ~ as.character(UpaTo))) %>%
  left_join(., district_upa %>% select(District, upazila_code), by = c("UpaFrom" = "upazila_code")) %>%
  rename(DistrictFrom = District) %>% 
  left_join(., district_upa %>% select(District, upazila_code), by = c("UpaTo" = "upazila_code")) %>%
  rename(DistrictTo = District)

trans.sub.upa.long %>% filter(is.na(DistrictFrom))

#saveRDS(trans.sub.upa.long, file = here("data","cleaned_2017_CDR", "trans.sub.upa.long.rds"))
#trans.sub.upa.long = readRDS(here("data","cleaned_2017_CDR", "trans.sub.upa.long.rds"))


## missing data
missing_op1_upa = c("302602", "302605", "302606", "302609", "302610", "302611", "302624",
                        "302629", "302632", "302630", "302633", "302637", "302663", "302665",
                        "302667", "302674", "302675", "302680", "302692", "302693", "302696")

#### subscribers ####
subsc <- read.csv(here("data","cleaned_2017_CDR","PopulationCountsUnion.csv"), header = T, sep = "|") %>%
  mutate(Union = as.character(Union))

no.subsc = subsc %>% pivot_longer(!Union, names_to = "date", values_to = "value") %>%
  group_by(date) %>%
  summarise(sum = sum(value, na.rm=T))

max(no.subsc$sum)

# Imputation for upazilas in Dhaka missing in 2020 data ----
# remove Eid days June 25-June 27 and September 1-3
# # eid dates based on https://www.officeholidays.com/countries/bangladesh/2017
trans.sub.upa.impute = trans.sub.upa.long %>%
  filter(DistrictFrom == "Dhaka" & UpaFrom %in% missing_op1_upa) %>%
  filter(!(date %in% c("2017-06-25", "2017-06-26", "2017-06-27", "2017-09-01", "2017-09-02", "2017-09-03")))

# check volumes 
weekend_df_op1 = readRDS(here("out", "2020_cdr", "data_out", "weekend_df_op1.rds"))

check.2020.dhaka = weekend_df_op1 %>%
  left_join(., district_upa %>% select(District, upazila_code), by = c("UpaFrom" = "upazila_code")) %>%
  rename(DistrictFrom = District) %>%
  filter(DistrictFrom == "Dhaka" & !(UpaFrom %in% missing_op1_upa) &
           !(UpaTo %in% missing_op1_upa)) %>% 
  mutate(year = "vol.2020")# eid dates already removed in weekend_df_op1

check.2017.dhaka = trans.sub.upa.long %>%
  filter(DistrictFrom == "Dhaka" & !(UpaFrom %in% missing_op1_upa) &
           !(UpaTo %in% missing_op1_upa)) %>%
  filter(!(date %in% c("2017-06-25", "2017-06-26", "2017-06-27", "2017-09-01", "2017-09-02", "2017-09-03"))) %>%
  group_by(UpaFrom, UpaTo, weekend, DistrictFrom) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>%
  select(UpaFrom, UpaTo, weekend, value, DistrictFrom) %>%
  mutate(year = "vol.2017") %>%
  ungroup()

compare.volumes.dhaka = bind_rows(check.2020.dhaka, check.2017.dhaka) %>%
  pivot_wider(names_from = year, values_from = value) %>%
  mutate(multiplier = vol.2020/vol.2017) %>%
  filter(is.finite(multiplier) & multiplier != 0)

summary(compare.volumes.dhaka$multiplier)
median_multiplier = summary(compare.volumes.dhaka$multiplier)[3] #1.332365 

compare.volumes.dhaka %>%
  ggplot(.) +
  geom_histogram(aes(x = multiplier))

# Impute proportions out of missing Dhaka upazilas: use 2017 proportions since lack data ====
weekend.mobility.impute = trans.sub.upa.impute %>% 
  group_by(weekend, UpaFrom, UpaTo) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>%
  ungroup()

# these are added as columns to 2020 matrix
# do not need to consider different absolute values from 2020 since just using proportions
# for trip volumes in, use same multiplier for each cell in column
prop.out.weekday.2017 = get_prop_m(get_m(weekend.mobility.impute, weekday = TRUE))
prop.out.weekend.2017 = get_prop_m(get_m(weekend.mobility.impute, weekday = FALSE))

weekday.2017.m = get_m(weekend.mobility.impute, weekday = TRUE)
weekend.2017.m = get_m(weekend.mobility.impute, weekday = FALSE)

#saveRDS(prop.out.weekday.2017, file = here("data", "cleaned_2017_cdr", "prop.out.weekday.2017.rds"))
#saveRDS(prop.out.weekend.2017, file = here("data", "cleaned_2017_cdr", "prop.out.weekend.2017.rds"))
colSums(prop.out.weekday.2017)
colSums(prop.out.weekday.2017)

# save raw counts to get uncertainty intervals
saveRDS(weekday.2017.m, file = here("data", "cleaned_2017_cdr", "weekday.2017.m.rds"))
saveRDS(weekend.2017.m, file = here("data", "cleaned_2017_cdr", "weekend.2017.m.rds"))

# Impute proportions into missing Dhaka upazilas ====
# add these trip totals to rows in 2020 matrix and then recompute prop
wekeend.allupa = trans.sub.upa.long %>%
  filter(!(date %in% c("2017-06-25", "2017-06-26", "2017-06-27", "2017-09-01", "2017-09-02", "2017-09-03"))) %>%
  group_by(weekend, UpaFrom, UpaTo) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>%
  ungroup()
  
load(here("data","full_pop_age_gender.Rdata"))
pop <- full_pop %>% 
  mutate(UPAZILA_CODE = as.character(upa_codes),
         pop_over60 = m_60 + m_65 + m_70 + m_75 + m_80 + f_60 + f_65 + f_70 + f_75 + f_80) %>%
  mutate(upazila_code = case_when(substr(UPAZILA_CODE, 1, 4) == "4539" ~ paste0("3039", substr(UPAZILA_CODE, 5, 7)),
                                  substr(UPAZILA_CODE, 1, 4) == "4561" ~ paste0("3061", substr(UPAZILA_CODE, 5, 7)),
                                  substr(UPAZILA_CODE, 1, 4) == "4572" ~ paste0("3072", substr(UPAZILA_CODE, 5, 7)),
                                  substr(UPAZILA_CODE, 1, 4) == "4589" ~ paste0("3089", substr(UPAZILA_CODE, 5, 7)),
                                  TRUE ~ as.character(UPAZILA_CODE))) %>% 
  dplyr::select(upazila_code, total_pop, pop_over60) %>%
  rename(pop = total_pop) %>%
  left_join(., district_upa, by = "upazila_code")

trips.in.weekday.2017 = wekeend.allupa %>%
  filter(weekend == FALSE & UpaFrom != "557384" & UpaTo != "557384") %>% 
  filter(UpaTo %in% missing_op1_upa) %>%
  mutate(value_multiplier = value*median_multiplier) %>%
  select(UpaFrom, UpaTo, value_multiplier) %>%
  pivot_wider(names_from = UpaFrom, values_from = value_multiplier) %>%
  remove_rownames %>% column_to_rownames(var="UpaTo") %>%
  data.matrix() 

trips.in.weekend.2017 = wekeend.allupa %>%
  filter(weekend == TRUE & UpaFrom != "557384" & UpaTo != "557384") %>% 
  filter(UpaTo %in% missing_op1_upa) %>%
  mutate(value_multiplier = value*median_multiplier) %>%
  select(UpaFrom, UpaTo, value_multiplier) %>%
  pivot_wider(names_from = UpaFrom, values_from = value_multiplier) %>%
  remove_rownames %>% column_to_rownames(var="UpaTo") %>%
  data.matrix()

saveRDS(trips.in.weekday.2017, file = here("data", "cleaned_2017_cdr", "trips.in.weekday.2017.withmult.rds"))
saveRDS(trips.in.weekend.2017, file = here("data", "cleaned_2017_cdr", "trips.in.weekend.2017.withmult.rds"))

