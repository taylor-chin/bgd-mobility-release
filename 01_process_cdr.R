# 01_process_cdr.R

# Load dependencies for this file
require(here)
library(tidyverse) ## added
library(lubridate)

rm(list = ls())

source(here("code", "utils.R")) #matrix functions in here

# load pre-saved data from # 00_mask_cdr.R
full_df_mobility <- readRDS(here("out", "cleaned_cdr_2020_masked.rds"))

# Create datasets ---- 
travelto <- full_df_mobility %>% 
  filter(UpaFrom != UpaTo)  %>% 
  group_by(UpaTo, operator) %>%  
  rename(upazila_code = UpaTo) %>% 
  summarize_if(is.numeric,sum, na.rm = T)  

travelfrom <- full_df_mobility %>% 
  filter(UpaFrom != UpaTo) %>% 
  group_by(UpaFrom, operator) %>%   
  rename(upazila_code = UpaFrom) %>% 
  summarize_if(is.numeric,sum, na.rm = T) 


stay_put <- full_df_mobility %>% 
  filter(UpaFrom == UpaTo) %>% 
  group_by(UpaFrom, operator) %>%
  rename(upazila_code = UpaFrom) %>% 
  summarize_if(is.numeric,sum, na.rm = T) 

subsc_tmp <- full_df_mobility %>% 
  group_by(UpaTo, operator) %>%  
  rename(upazila_code = UpaTo) %>% 
  summarize_if(is.numeric,sum, na.rm = T) 

# missing these dates
travelto$`2020-05-09` <- travelfrom$`2020-05-09`<- stay_put$`2020-05-08`<- NA # incomplete data

travelto$`2020-05-25` <- travelfrom$`2020-05-25` <- stay_put$`2020-05-25`  <-NA # missing
travelto$`2020-05-26` <- travelfrom$`2020-05-26` <- stay_put$`2020-05-26`  <-NA # missing

travelto$`2020-06-08` <- travelfrom$`2020-06-08` <- stay_put$`2020-06-08` <- NA # missing
travelto$`2020-06-09` <- travelfrom$`2020-06-09` <- stay_put$`2020-06-09` <- NA # missing

travelto$`2020-06-15` <- travelfrom$`2020-06-15` <- stay_put$`2020-06-15` <- NA # incomplete data

travelto$`2020-08-01` <- travelfrom$`2020-08-01` <- stay_put$`2020-08-01` <- NA # incomplete data
travelto$`2020-08-04` <- travelfrom$`2020-08-04` <- stay_put$`2020-08-04` <- NA # incomplete data

subsc_tmp2 <- as.data.frame(as.matrix(subsc_tmp[,3:ncol(subsc_tmp)]) - as.matrix(travelfrom[,3:ncol(travelfrom)]))
subsc <- as.data.frame(cbind(upazila_code = subsc_tmp$upazila_code, operator = subsc_tmp$operator, subsc_tmp2))

subsc$`2020-05-09` <- NA
subsc$`2020-05-25` <- NA
subsc$`2020-05-26` <- NA
subsc$`2020-06-08` <- NA
subsc$`2020-06-09` <- NA
subsc$`2020-06-15` <- NA
subsc$`2020-08-01` <- NA
subsc$`2020-08-04` <- NA

travelto$average = rowMeans(travelto[,3:ncol(travelto)], na.rm= TRUE)
travelfrom$average = rowMeans(travelfrom[,3:ncol(travelfrom)], na.rm= TRUE)
subsc$average <- rowMeans(subsc[,3:ncol(subsc)], na.rm= TRUE)
stay_put$average <- rowMeans(stay_put[,3:ncol(stay_put)], na.rm= TRUE)

# Average number of subscribers across 3 providers - 100M ---
subsc %>% 
  select(!average) %>%
  pivot_longer(!c(upazila_code, operator), names_to = "date", values_to = "subsc") %>%
  group_by(operator, date) %>% # aggregate over upazilas
  summarise(subsc = sum(subsc)) %>%
  filter(subsc != 0) %>% # get rid of missing dates
  group_by(operator) %>%
  summarise(ave_operator_subsc = mean(subsc)) %>% # daily average per operator
  summarise(sum(ave_operator_subsc)) 

# Create mobility matrices ----
## weekend/weekday matrices ====
# eid dates based on https://publicholidays.com.bd/2020-dates/
eid.dates = c("2020-05-24","2020-05-25", "2020-05-26", "2020-07-31", "2020-08-01", "2020-08-02")

mobility_m_op1 = full_df_mobility %>%
  filter(operator == "Operator 1") %>% 
  pivot_longer(`2020-04-28`:`2020-09-05`, names_to = "date", values_to = "value") %>%
  filter(!(date %in% eid.dates)) %>%
  mutate(weekend = ifelse(wday(as.Date(date)) == 1, 1, 0) | ifelse(wday(as.Date(date)) == 7, 1, 0))

mobility_m_op2 = full_df_mobility %>%
  filter(operator == "Operator 2") %>% 
  pivot_longer(`2020-04-28`:`2020-09-05`, names_to = "date", values_to = "value") %>%
  filter(!(date %in% eid.dates)) %>%
  mutate(weekend = ifelse(wday(as.Date(date)) == 1, 1, 0) | ifelse(wday(as.Date(date)) == 7, 1, 0))

mobility_m_op3 = full_df_mobility %>%
  filter(operator == "Operator 3") %>% 
  pivot_longer(`2020-04-28`:`2020-09-05`, names_to = "date", values_to = "value") %>%
  filter(!(date %in% eid.dates)) %>%
  mutate(weekend = ifelse(wday(as.Date(date)) == 1, 1, 0) | ifelse(wday(as.Date(date)) == 7, 1, 0))

#saveRDS(mobility_m_op1, here("out", "2020_cdr", "data_out", "mobility_m_op1.rds"))
#saveRDS(mobility_m_op2, here("out", "2020_cdr", "data_out", "mobility_m_op2.rds"))
#saveRDS(mobility_m_op3, here("out", "2020_cdr", "data_out", "mobility_m_op3.rds"))

mobility_m_op1 = readRDS(here("out", "2020_cdr", "data_out", "mobility_m_op1.rds"))
mobility_m_op2 = readRDS(here("out", "2020_cdr", "data_out", "mobility_m_op2.rds"))
mobility_m_op3 = readRDS(here("out", "2020_cdr", "data_out", "mobility_m_op3.rds"))

# summarise by weekend/weekday
weekend_df_op1 = mobility_m_op1 %>%
  group_by(weekend, UpaFrom, UpaTo) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>%
  ungroup()

weekend_df_op2 = mobility_m_op2 %>%
  group_by(weekend, UpaFrom, UpaTo) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>%
  ungroup()

weekend_df_op3 = mobility_m_op3 %>%
  group_by(weekend, UpaFrom, UpaTo) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>%
  ungroup()

#saveRDS(weekend_df_op1, here("out", "2020_cdr", "data_out", "weekend_df_op1.rds"))
#saveRDS(weekend_df_op2, here("out", "2020_cdr", "data_out", "weekend_df_op2.rds"))
#saveRDS(weekend_df_op3, here("out", "2020_cdr", "data_out", "weekend_df_op3.rds"))

weekend_df_op1 = readRDS(here("out", "2020_cdr", "data_out", "weekend_df_op1.rds"))
weekend_df_op2 = readRDS(here("out", "2020_cdr", "data_out", "weekend_df_op2.rds"))
weekend_df_op3 = readRDS(here("out", "2020_cdr", "data_out", "weekend_df_op3.rds"))

# get lists of true NA upazilas ----
op1_true_na_upa = get_na_upa_list(weekend_df_op1)
op2_true_na_upa = get_na_upa_list(weekend_df_op2)
op3_true_na_upa = get_na_upa_list(weekend_df_op3)

# differentiate na 
op1_df_withna_weekday = differentiate_zero_na_df(df = weekend_df_op1, true_na_upa = op1_true_na_upa)
op2_df_withna_weekday = differentiate_zero_na_df(df = weekend_df_op2, true_na_upa = op2_true_na_upa)
op3_df_withna_weekday = differentiate_zero_na_df(df = weekend_df_op3, true_na_upa = op3_true_na_upa)

# get m with imputed missing dhaka upazilas for Operators 1 and 3 only
# these functions use files saved from 02_impute_cdr_2017
op1_m_weekday_counts = get_impute_dhaka_m(op1_df_withna_weekday, weekday=TRUE)
op1_m_weekend_counts = get_impute_dhaka_m(op1_df_withna_weekday, weekday=FALSE)
op2_m_weekday_counts = get_m(op2_df_withna_weekday, weekday = TRUE)
op2_m_weekend_counts = get_m(op2_df_withna_weekday, weekday = FALSE)
op3_m_weekday_counts = get_impute_dhaka_m(op3_df_withna_weekday, weekday=TRUE)
op3_m_weekend_counts = get_impute_dhaka_m(op3_df_withna_weekday, weekday=FALSE)
# 
#saveRDS(op1_m_weekday_counts, here("out", "2020_cdr", "matrices", "op1_m_weekday_counts.rds"))
#saveRDS(op1_m_weekend_counts, here("out", "2020_cdr", "matrices", "op1_m_weekend_counts.rds"))
#saveRDS(op2_m_weekday_counts, here("out", "2020_cdr", "matrices", "op2_m_weekday_counts.rds"))
#saveRDS(op2_m_weekend_counts, here("out", "2020_cdr", "matrices", "op2_m_weekend_counts.rds"))
#saveRDS(op3_m_weekday_counts, here("out", "2020_cdr", "matrices", "op3_m_weekday_counts.rds"))
#saveRDS(op3_m_weekend_counts, here("out", "2020_cdr", "matrices", "op3_m_weekend_counts.rds"))

op1_m_weekday_counts = readRDS(here("out", "2020_cdr", "matrices", "op1_m_weekday_counts.rds"))
op1_m_weekend_counts = readRDS(here("out", "2020_cdr", "matrices", "op1_m_weekend_counts.rds"))
op2_m_weekday_counts = readRDS(here("out", "2020_cdr", "matrices", "op2_m_weekday_counts.rds"))
op2_m_weekend_counts = readRDS(here("out", "2020_cdr", "matrices", "op2_m_weekend_counts.rds"))
op3_m_weekday_counts = readRDS(here("out", "2020_cdr", "matrices", "op3_m_weekday_counts.rds"))
op3_m_weekend_counts = readRDS(here("out", "2020_cdr", "matrices", "op3_m_weekend_counts.rds"))


# sample using binom.bayes ----
op1_weekday_beta_param = get_beta_distr_param(op1_m_weekday_counts)
op1_weekend_beta_param = get_beta_distr_param(op1_m_weekend_counts)
op2_weekday_beta_param = get_beta_distr_param(op2_m_weekday_counts)
op2_weekend_beta_param = get_beta_distr_param(op2_m_weekend_counts)
op3_weekday_beta_param = get_beta_distr_param(op3_m_weekday_counts)
op3_weekend_beta_param = get_beta_distr_param(op3_m_weekend_counts)

saveRDS(op1_weekday_beta_param, here("out", "2020_cdr", "beta_distr", "op1_weekday_beta_param.rds"))
saveRDS(op1_weekend_beta_param, here("out", "2020_cdr", "beta_distr", "op1_weekend_beta_param.rds"))
saveRDS(op2_weekday_beta_param, here("out", "2020_cdr", "beta_distr", "op2_weekday_beta_param.rds"))
saveRDS(op2_weekend_beta_param, here("out", "2020_cdr", "beta_distr", "op2_weekend_beta_param.rds"))
saveRDS(op3_weekday_beta_param, here("out", "2020_cdr", "beta_distr", "op3_weekday_beta_param.rds"))
saveRDS(op3_weekend_beta_param, here("out", "2020_cdr", "beta_distr", "op3_weekend_beta_param.rds"))


## calculate district matrices first and then beta param
op1.district.m1 = get_district_m(op1_m_weekday_counts)
op1.district.m2 = get_district_m(op1_m_weekend_counts)
op2.district.m1 = get_district_m(op2_m_weekday_counts)
op2.district.m2 = get_district_m(op2_m_weekend_counts)
op3.district.m1 = get_district_m(op3_m_weekday_counts)
op3.district.m2 = get_district_m(op3_m_weekend_counts)

op1_weekday_beta_param_district = get_beta_distr_param(op1.district.m1)
op1_weekend_beta_param_district = get_beta_distr_param(op1.district.m2)
op2_weekday_beta_param_district = get_beta_distr_param(op2.district.m1)
op2_weekend_beta_param_district = get_beta_distr_param(op2.district.m2)
op3_weekday_beta_param_district = get_beta_distr_param(op3.district.m1)
op3_weekend_beta_param_district = get_beta_distr_param(op3.district.m2)


# save no uncertainty matrices
upa_pop = pop
op1.sym.m1 = get_norm_sym_m(op1_m_weekday_counts, pop_input = upa_pop, mobility_source = "Operator 1", sym = TRUE)
op1.district.m1 = get_district_m(op1.sym.m1)
op1.sym.m2 = get_norm_sym_m(op1_m_weekend_counts, pop_input = upa_pop, mobility_source = "Operator 1", sym = TRUE)
op1.district.m2 = get_district_m(op1.sym.m2)
op2.sym.m1 = get_norm_sym_m(op2_m_weekday_counts, pop_input = upa_pop, mobility_source = "Operator 2", sym = TRUE)
op2.district.m1 = get_district_m(op2.sym.m1)
op2.sym.m2 = get_norm_sym_m(op2_m_weekend_counts, pop_input = upa_pop, mobility_source = "Operator 2", sym = TRUE)
op2.district.m2 = get_district_m(op2.sym.m2)
op3.sym.m1 = get_norm_sym_m(op3_m_weekday_counts, pop_input = upa_pop, mobility_source = "Operator 3", sym = TRUE)
op3.district.m1 = get_district_m(op3.sym.m1)
op3.sym.m2 = get_norm_sym_m(op3_m_weekend_counts, pop_input = upa_pop, mobility_source = "Operator 3", sym = TRUE)
op3.district.m2 = get_district_m(op3.sym.m2)

saveRDS(op1.district.m1, here("out", "2020_cdr", "matrices", "op1_m_weekday_district_withNA_symmetric_dhakamult.rds"))
saveRDS(op1.district.m2, here("out", "2020_cdr", "matrices", "op1_m_weekend_district_withNA_symmetric_dhakamult.rds"))
saveRDS(op2.district.m1, here("out", "2020_cdr", "matrices", "op2_m_weekday_district_withNA_symmetric_dhakamult.rds"))
saveRDS(op2.district.m2, here("out", "2020_cdr", "matrices", "op2_m_weekend_district_withNA_symmetric_dhakamult.rds"))
saveRDS(op3.district.m1, here("out", "2020_cdr", "matrices", "op3_m_weekday_district_withNA_symmetric_dhakamult.rds"))
saveRDS(op3.district.m2, here("out", "2020_cdr", "matrices", "op3_m_weekend_district_withNA_symmetric_dhakamult.rds"))

saveRDS(op1.sym.m1, here("out", "2020_cdr", "matrices", "op1_m_weekday_withNA_symmetric_dhakamult.rds"))
saveRDS(op1.sym.m2, here("out", "2020_cdr", "matrices", "op1_m_weekend_withNA_symmetric_dhakamult.rds"))
saveRDS(op2.sym.m1, here("out", "2020_cdr", "matrices", "op2_m_weekday_withNA_symmetric_dhakamult.rds"))
saveRDS(op2.sym.m2, here("out", "2020_cdr", "matrices", "op2_m_weekend_withNA_symmetric_dhakamult.rds"))
saveRDS(op3.sym.m1, here("out", "2020_cdr", "matrices", "op3_m_weekday_withNA_symmetric_dhakamult.rds"))
saveRDS(op3.sym.m2, here("out", "2020_cdr", "matrices", "op3_m_weekend_withNA_symmetric_dhakamult.rds"))


# DIAGNOSTICS ----
### . 2017 matrix to check volumes ====
### used in 02_impute_cdr_2017.R
trans.sub.upa.long = readRDS(here("data","cleaned_2017_CDR", "trans.sub.upa.long.rds"))
weekend_df_2017 = trans.sub.upa.long %>%
  mutate(weekend = ifelse(wday(as.Date(date)) == 1, 1, 0) | ifelse(wday(as.Date(date)) == 7, 1, 0)) %>%
  group_by(weekend, UpaFrom, UpaTo) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>%
  ungroup()

op1_m_weekday_2017 = get_weekend_m_2017(weekend_df_2017, weekday = TRUE)
op1_m_weekend_2017 = get_weekend_m_2017(weekend_df_2017, weekday = FALSE)

#saveRDS(op1_m_weekday_2017, file = here("out", "2020_cdr", "data_out", "op1_m_weekday_2017.rds"))
#saveRDS(op1_m_weekend_2017, file = here("out", "2020_cdr", "data_out", "op1_m_weekend_2017.rds"))
