# 08_supp_figs.R

# Load dependences ----
library(here)
library(tidyverse)
library(RColorBrewer)
library(MetBrewer)
library(cowplot)
library(gtools) # quantcut

source(here("code", "utils.R"))

# Load data ----
# rename pop data
upa_pop = pop

## . CDR and matrices ====
# Load 2020 operator 1 data for Fig S7, S10 - created in 01_compare_cdr.R
weekend_df_op1 = readRDS(here("out", "2020_cdr", "data_out", "weekend_df_op1.rds"))
op1_m_weekday_counts = readRDS(here("out", "2020_cdr", "matrices", "op1_m_weekday_counts.rds"))
op1_df = get_df_from_m(op1_m_weekday_counts)
# Load 2017 operator 1 data for Fig S7, S10 - created in 02_impute_cdr_2017.R
trans.sub.upa.long = readRDS(here("data","cleaned_2017_CDR", "trans.sub.upa.long.rds"))

# for Figs S2-3, load final matrices - created in 01_compare_cdr.R
op1_m_upa = readRDS(here("out", "2020_cdr", "matrices", "op1_m_weekday_withNA_symmetric_dhakamult.rds"))
op2_m_upa = readRDS(here("out", "2020_cdr", "matrices", "op2_m_weekday_withNA_symmetric_dhakamult.rds"))
op3_m_upa = readRDS(here("out", "2020_cdr", "matrices", "op3_m_weekday_withNA_symmetric_dhakamult.rds"))
grav_lit_upa_m = readRDS(here("out", "metapop_model", "grav_mat.upa.w.lit.params.rds"))

## For Fig 8 date range
mobility_m_op1 = readRDS(here("out", "2020_cdr", "data_out", "mobility_m_op1.rds"))
mobility_m_op2 = readRDS(here("out", "2020_cdr", "data_out", "mobility_m_op2.rds"))
mobility_m_op3 = readRDS(here("out", "2020_cdr", "data_out", "mobility_m_op3.rds"))

## . Meta data number of subscribers ====
# from 03_fb_movement.R
# for Fig S1
pop_clean <- readRDS(file = here("data", "fb-pop-data", "fb-pop-cleaned-mapped.rds"))

pop_sub = pop_clean %>% 
  select(date_time, n_baseline, n_crisis, date.time, time, date, ADM2_EN, ADM2_PCODE) %>%
  mutate(n_baseline = as.numeric(n_baseline), 
         n_crisis = as.numeric(n_crisis))

pop_summary = pop_sub %>%
  as.data.frame() %>%
  # sum across tiles in district
  group_by(ADM2_EN, date, date.time) %>%
  summarise(n_baseline = sum(n_baseline, na.rm=T), 
            n_crisis = sum(n_crisis, na.rm=T)) %>%
  # take mean across 3 time points per day
  group_by(ADM2_EN, date) %>%
  summarise(n_baseline = sum(n_baseline, na.rm=T), 
            n_crisis = sum(n_crisis, na.rm=T),           
            perc_change = (n_crisis-n_baseline)*100/n_baseline) %>%
  ungroup()

## . district matrices w/o uncertainty ====
## For Fig S9
op1_m_district = readRDS(here("out", "2020_cdr", "matrices", "op1_m_weekday_district_withNA_symmetric_dhakamult.rds"))
op2_m_district = readRDS(here("out", "2020_cdr", "matrices", "op2_m_weekday_district_withNA_symmetric_dhakamult.rds"))
op3_m_district = readRDS(here("out", "2020_cdr", "matrices", "op3_m_weekday_district_withNA_symmetric_dhakamult.rds"))
fb_m_district = readRDS(here("out", "FB-analysis", "fb_weekday_m_withNA_symmetric.rds"))

## . Demographic data ====
# for Fig S1
## read in subscribers # from 01_process_cdr.R
subsc = readRDS(here("out", "2020_cdr", "subsc_2020.rds"))

## For Fig S1
income <- read.csv(here("data", "misc", "income_by_union.csv"), header = TRUE)
income_district = income %>%
  mutate(ADM2_EN = case_when(District == "Cox'S Bazar" ~ "Cox's Bazar", 
                             TRUE ~ as.character(District))) %>%
  group_by(ADM2_EN) %>%
  summarise(income_mean = mean(income, na.rm=T)) %>%
  mutate(income_cat = quantcut(income_mean, q=4, na.rm=TRUE))

# nightlight data as proxy of urbanicity
# source: http://geo.aiddata.org/query/#!/ (VIIRS Nighttime Lights)
# https://data.ngdc.noaa.gov/instruments/remote-sensing/passive/spectrometers-radiometers/imaging/viirs/dnb_composites/v10/README_dnb_composites_v1.txt

nl <- read.csv(here("data", "misc", "bgd_nightlight_adm2.csv"), header = TRUE) %>%
  mutate(ADM2_EN = case_when(NAME_2 == "Borgona" ~ "Barguna",
                             NAME_2 == "Jhalakati" ~ "Jhalokati",
                             NAME_2 == "Bandarbon" ~ "Bandarban",
                             NAME_2 == "Brahmanbaria" ~ "Brahamanbaria",
                             NAME_2 == "Khagrachari" ~ "Khagrachhari",
                             NAME_2 == "Parbattya Chattagram" ~ "Rangamati",
                             NAME_2 == "Gopalgonj" ~ "Gopalganj",
                             NAME_2 == "Manikgonj" ~ "Manikganj",
                             NAME_2 == "Munshigonj" ~ "Munshiganj",
                             NAME_2 == "Naray Angonj" ~ "Narayanganj",
                             NAME_2 == "Jaipurhat" ~ "Joypurhat",
                             NAME_2 == "Sirajgonj" ~ "Sirajganj",
                             NAME_2 == "Gaibanda" ~ "Gaibandha",
                             NAME_2 == "Rongpur" ~ "Rangpur",
                             NAME_2 == "Hobiganj" ~ "Habiganj",
                             NAME_2 == "Moulvibazar" ~ "Maulvibazar",
                             NAME_2 == "Sun Amgonj" ~ "Sunamganj",
                             NAME_2 == "Narshingdi" ~ "Narsingdi",
                             NAME_2 == "Nasirabad" ~ "Mymensingh",
                             NAME_2 == "Choua Danga" ~ "Chuadanga",
                             NAME_2 == "Kustia" ~ "Kushtia",
                             NAME_2 == "Shatkhira" ~ "Satkhira",
                             TRUE ~ as.character(NAME_2))) %>%
  dplyr::select(ADM2_EN, nl_mean_2018 = viirs_vcmcfg_dnb_composites_v10_yearly_max.2018.mean) %>%
  mutate(nl_cat = quantcut(nl_mean_2018, q=4, na.rm=TRUE)) 

## . Simulation results ----
# For Figs S4
all.results.bind = readRDS(here("out", "metapop_model", "cluster_results", "results.district.beta.avgsim.new.rds"))
# For Fig S6
all.results.bind.upa = readRDS(here("out", "metapop_model", "cluster_results", "results.upa.beta.avgsim.rds"))

# Create figures ----
# Fig S1. Meta vs. Operator 1 user demographics ----
## process Meta 
meta_n_crisis_pop = pop_summary %>%
  filter(!is.na(ADM2_EN)) %>%
  group_by(ADM2_EN) %>%
  summarise(fb_pop_crisis = mean(n_crisis, na.rm=T))

meta_n_baseline_pop = pop_summary %>%
  filter(!is.na(ADM2_EN)) %>%
  group_by(ADM2_EN) %>%
  summarise(fb_pop_baseline = mean(n_baseline, na.rm=T)) # average across dates

meta_pop = meta_n_crisis_pop %>%
  left_join(., meta_n_baseline_pop, by = "ADM2_EN") %>%
  left_join(., nl, by = "ADM2_EN") %>%
  left_join(., income_district, by = "ADM2_EN") %>%
  select(-c(nl_cat, income_cat, fb_pop_baseline)) %>%
  rename(District = ADM2_EN,
         users = fb_pop_crisis) %>%
  mutate(source = "Meta") %>%
  filter(District != "Jhalokati")

# process CDR 
subsc_join = subsc %>%
  filter(upazila_code != "557384") %>% # remove missing upa
  filter(operator == "Operator 1") %>%
  select(upazila_code, average) %>%
  left_join(., district_upa %>% select(upazila_code, District), by = c("upazila_code")) %>%
  group_by(District) %>% 
  summarise(subsc =sum(average, na.rm=T)) %>%  # sum across upazilas
  left_join(., nl, by = c("District" = "ADM2_EN")) %>%
  left_join(., income_district, by = c("District"= "ADM2_EN")) %>%
  select(-c(nl_cat, income_cat)) %>%
  rename(users = subsc) %>%
  mutate(source = "Operator 1") %>%
  filter(District != "Jhalokati")

# join together and pivot
demog = bind_rows(meta_pop, subsc_join) %>%
  pivot_longer(c(nl_mean_2018, income_mean), names_to = "metric", values_to = "value") %>%
  mutate(metric = case_when(metric == "income_mean" ~ "Household Income",
                            metric == "nl_mean_2018" ~ "Nightlight Level")) %>%
  group_by(source, metric) %>%
  mutate(cor = round(cor(log(users), log(value)), 2))

demog.plt = ggplot(demog) +
  geom_point(aes(x = log(users), y = log(value), group = District)) +
  geom_text(aes(x = -Inf, y = Inf, hjust = -0.2, vjust = 1.2,  
                label = paste0("r=", cor))) +
  labs(x = "Log average district number of users", 
       y = "Log average district value") +
  facet_grid(metric~source, scales = "free") +
  scale_x_continuous(labels = scales::comma) +
  theme_classic() +
  theme(axis.line = element_line(),
        axis.ticks = element_line(),
        ## Strips for facet_wrap
        strip.text=element_text(size=8),
        strip.background=element_rect(fill="#f0f0f0"),
        legend.position = "none") +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

ggsave(demog.plt, file =  here("out", "supplement_figs", "s1.user.demog.png"), width =8, height = 5)

# Fig S2A. Upazila level heat map ----
upa_compare_df = bind_rows(get_df_from_m(op1_m_upa) %>% mutate(source = "Operator 1"), 
                           get_df_from_m(op2_m_upa) %>% mutate(source = "Operator 2"),
                           get_df_from_m(op3_m_upa) %>% mutate(source = "Operator 3"),
                           get_df_from_m(grav_lit_upa_m) %>% mutate(source = "Gravity Model")) %>%
  mutate(source = factor(source, levels = c("Operator 1", "Operator 2", "Operator 3", "Gravity Model")))

#latlong from utils
latlong_recon = latlong %>%
  mutate(upazila_code = case_when(substr(upazila_code, 1, 4) == "4539" ~ paste0("3039", substr(upazila_code, 5, 7)),
                                  substr(upazila_code, 1, 4) == "4561" ~ paste0("3061", substr(upazila_code, 5, 7)),
                                  substr(upazila_code, 1, 4) == "4572" ~ paste0("3072", substr(upazila_code, 5, 7)),
                                  substr(upazila_code, 1, 4) == "4589" ~ paste0("3089", substr(upazila_code, 5, 7)),
                                  TRUE ~ as.character(upazila_code))) 

df_plt = upa_compare_df %>%
  mutate(UpaTo = factor(UpaTo, levels = (latlong_recon %>% pull(upazila_code))), 
         UpaFrom = factor(UpaFrom, levels = (latlong_recon %>% pull(upazila_code)))) %>%
  mutate(log_value = log(value)) %>%
  mutate(log_value_NA = case_when(is.infinite(log_value) ~ 10000, # random high value for plottin
                                  TRUE ~ log_value))
df_plt %>% filter(log_value_NA == 10000)

min_value = floor(summary((df_plt %>% filter(!is.infinite(log_value)))$log_value)[1])
max_value = ceiling(summary((df_plt %>% filter(!is.infinite(log_value)))$log_value)[6])


heat_upa = ggplot(df_plt %>% filter(source != "Gravity Model")) +
  geom_tile(aes(x = UpaFrom, y = UpaTo, fill = log_value_NA)) +
  scale_fill_gradientn(colors = c(brewer.pal(9, "Blues"),"#999999"), na.value = "black",
                       limit = c(min_value,max_value),oob=scales::squish) +
  facet_wrap(~source, scales = "free", ncol = 3) +
  export_theme +
  theme(axis.title.x = element_text(vjust = 1),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y= element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = "Origin Upazila", y = "Destination Upazila ", fill = "Average log \nnumber of trips")

#ggsave(heat_upa, file =  here("out", "descriptive", "s4.png"), width = 12, height = 6.5)

# Fig S2B. Histogram proportion staying put ----
prop.compare = upa_compare_df %>%
  mutate(same = ifelse(UpaTo == UpaFrom, 1, 0)) %>%
  group_by(UpaFrom, source) %>% 
  mutate(sum = sum(value, na.rm=T)) %>%
  filter(same == 1) %>%
  group_by(UpaFrom, source) %>%
  mutate(prop = value/sum)

hist.prop = prop.compare %>%
  filter(source != "Gravity Model") %>%
  ggplot() +
  geom_histogram(aes(x = prop), bins = 20) +
  facet_wrap(~source) +
  export_theme +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  theme(legend.text = element_blank(),
        legend.title = element_blank()) +
  labs(y = "Number of upazilas", x = "Proportion of people staying put in upazila")

# combine plots 
combined_plt_upa = plot_grid(heat_upa, hist.prop, nrow = 2, 
                             labels = c('A', 'B'), align = 'hv', axis = 'lr') 

ggsave(combined_plt_upa, file =  here("out", "supplement_figs", "s2.upamat.png"), 
       width = 9, height = 5)

# Fig S3. Characterize zeroes in matrices ----

# take symmetric upazila level matrices
# replace true missing NA rows and columns with NAs
# replace upper tri of mat with negative number to filter later
# melt to df
upa_compare_df_dedup = bind_rows(get_df_from_m(half_mat(get_mat_with_na(op1_m_upa, mobility_source = "Operator 1"))) %>% mutate(source = "Operator 1"), 
                                 get_df_from_m(half_mat(get_mat_with_na(op2_m_upa, mobility_source = "Operator 2"))) %>% mutate(source = "Operator 2"),
                                 get_df_from_m(half_mat(get_mat_with_na(op3_m_upa, mobility_source = "Operator 3"))) %>% mutate(source = "Operator 3"))
mutate(source = factor(source, levels = c("Operator 1", "Operator 2", "Operator 3")))

# population
zero.pop = upa_compare_df_dedup %>%
  filter(value != -10000) %>%
  filter(source == "Operator 2") %>%
  mutate(valuezero = ifelse(value ==0, "Zero", "Non-zero")) %>%
  left_join(., pop %>% select(upazila_code, pop), by = c("UpaFrom" = "upazila_code")) %>%
  ggplot() +
  geom_density(aes(x = pop, fill = valuezero), adjust = 5, alpha = 0.7) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values=met.brewer("Egypt", 4)) +
  export_theme +
  labs(y = "Density", x = "Upazila population", fill = "Operator 2 value")

op1.zero = upa_compare_df_dedup %>%
  filter(value != -10000) %>% # remove mat duplicated values
  pivot_wider(names_from = source, values_from= value) %>%
  mutate(op2_zero = ifelse(`Operator 2` == 0, "Zero", "Non-zero")) %>%
  filter(!is.na(`Operator 2`) & !is.na(`Operator 1`)) %>% 
  filter(`Operator 1` < quantile(`Operator 1`, probs = 0.95)) %>%
  ggplot() +
  geom_violin(aes(x = op2_zero, y = `Operator 1`, fill = op2_zero)) +
  scale_fill_manual(values=met.brewer("Egypt", 4)) +
  export_theme +
  labs(x = "", y = "Operator 1 value", fill = "Operator 2 value")

op3.zero = upa_compare_df_dedup %>%
  filter(value != -10000) %>% # remove mat duplicated values
  pivot_wider(names_from = source, values_from= value) %>%
  mutate(op2_zero = ifelse(`Operator 2` == 0, "Zero", "Non-zero")) %>%
  filter(!is.na(`Operator 3`) & !is.na(`Operator 2`)) %>% 
  filter(`Operator 3` < quantile(`Operator 3`, probs = 0.95)) %>%
  ggplot() +
  geom_violin(aes(x = op2_zero, y = `Operator 3`, fill = op2_zero)) +
  scale_fill_manual(values=met.brewer("Egypt", 4)) +
  export_theme +
  labs(x = "", y = "Operator 3 value", fill = "Operator 2 value")

combined_op2_plt = ggpubr::ggarrange(zero.pop, op1.zero, op3.zero, ncol = 3, common.legend = TRUE,
                                     labels = c('A', 'B', 'C'), align = 'hv') 

ggsave(combined_op2_plt, file =  here("out", "supplement_figs", "s3.op2zero.png"), width = 12, height = 3)


# Simulation supplementary figures ----
# Fig S4: district level grid (symptomatic only) ----
sim.grid = sumbgd %>%
  group_by(time, source, r0, seed_city) %>%
  summarise(mean = mean(value, na.rm=TRUE),
            lower = quantile(value, probs = 0.025, na.rm=TRUE),
            upper = quantile(value, probs = 0.975, na.rm=TRUE)) %>%
  mutate(source = factor(source, levels = c("Gravity model", "Operator 1", "Operator 2", "Operator 3", "Meta"))) %>%
  ggplot(.) +
  geom_line(aes(x = time, y = mean, colour = as.factor(source)), size = 0.65) +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper, fill = as.factor(source)), 
              alpha = 0.3) +
  facet_grid(r0~seed_city, scales = "free") +
  scale_colour_manual(values=met.brewer("Java", 5)) +
  scale_fill_manual(values=met.brewer("Java", 5)) +
  theme_bw()+
  theme(legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.key.size= unit(0.5, "cm"),
        legend.margin = margin(0,0,0,0, "cm"),
        ## Strips for facet_wrap
        strip.text=element_text(size=8),
        strip.background=element_rect(fill="#f0f0f0")) +
  labs(y = "Incidence", x = "Time (days)", colour = "Mobility matrix source", 
       fill = "Mobility matrix source")

ggsave(sim.grid,  file =  here("out", "supplement_figs", "s4.sims.district.new.png"), 
       width = 10, height = 6)

# Fig S5: district level outbreaks ----
district_outbreaks = all.results.bind %>%
  filter(state == "obs" & r0 == 1.3) %>%
  mutate(source = factor(source, levels = c("Gravity model", "Operator 1", "Operator 2", "Operator 3", "Meta"))) %>%
  ggplot(.) +
  geom_line(aes(x = time, y = mean, colour = as.factor(geo_code)), size = 0.65) +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper, fill = as.factor(geo_code)), 
              alpha = 0.3) +
  facet_grid(source~seed_city, scales = "free") +
 # scale_colour_manual(values=met.brewer("Java", 5)) +
#  scale_fill_manual(values=met.brewer("Java", 5)) +
  theme_bw()+
  theme(legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.key.size= unit(0.5, "cm"),
        legend.position = "bottom",
        legend.margin = margin(0,0,0,0, "cm"),
        ## Strips for facet_wrap
        strip.text=element_text(size=8),
        strip.background=element_rect(fill="#f0f0f0")) +
  labs(y = "Incidence", x = "Time (days)", colour = "District", 
       fill = "District")

ggsave(district_outbreaks,  file =  here("out", "supplement_figs", "s5.districtoutbreaks.png"), 
       width = 11, height = 10)

## . Values cited in results ====
# % missing upazilas
upa_compare_df %>%
  group_by(source, UpaFrom) %>%
  tally(is.na(value)) %>%
  filter(n == 544) %>% # number of upazilas 
  group_by(source) %>%
  count() %>%
  mutate(prop = n/544)

# proportion position or zero travel
upa_compare_df %>%
  filter(UpaFrom != UpaTo) %>%
  group_by(source) %>%
  tally(value > 0) %>% # either ==0, is.na or > 0
  mutate(betweenupa_comb = 544*544-544) %>%
  mutate(prop = n/betweenupa_comb)

upa_compare_df %>%
  filter(UpaFrom != UpaTo)
group_by(source) %>%
  tally(value == 0) %>% # either ==0, is.na or > 0
  mutate(betweenupa_comb = 544*544-544) %>%
  mutate(prop = n/betweenupa_comb)

# proportion staying put
prop.compare %>%
  group_by(source) %>%
  summarise(mean = mean(prop, na.rm=T),
            median = median(prop, na.rm=T),
            iqr = IQR(prop, na.rm=T))


# Fig S6: upazila level grid ----
sim.grid.upa = all.results.bind.upa %>%
  filter(state == "obs") %>%
  group_by(time, source, r0, seed_city) %>% #geo_code
  #  filter(source != "Gravity model") %>%
  # sum across districts 
  summarise(mean = sum(mean), 
            lower = sum(lower), 
            upper = sum(upper)) %>% 
  mutate(source = factor(source, levels = c("Gravity model", "Operator 1", "Operator 2", "Operator 3", "Meta"))) %>%
  ggplot(.) +
  geom_line(aes(x = time, y = mean, colour = as.factor(source)), size = 0.65) +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper, fill = as.factor(source)), 
              alpha = 0.3) +
  facet_grid(r0~seed_city, scales = "free") +
  scale_colour_manual(values=met.brewer("Java", 5)) +
  scale_fill_manual(values=met.brewer("Java", 5)) +
  theme_bw()+
  theme(legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.key.size= unit(0.5, "cm"),
        legend.margin = margin(0,0,0,0, "cm"),
        ## Strips for facet_wrap
        strip.text=element_text(size=8),
        strip.background=element_rect(fill="#f0f0f0")) +
  labs(y = "Incidence", x = "Time (days)", colour = "Mobility matrix source", 
       fill = "Mobility matrix source")

ggsave(sim.grid.upa,  file =  here("out", "supplement_figs", "s5.sims.upa.png"), 
       width = 10, height = 6)

# Fig S7. Calc overall proportion correlation 2020 vs. 2017 across all upazilas ----
# remove missing upazilas
missing_op1_upa = c("302602", "302605", "302606", "302609", "302610", "302611", "302624",
                    "302629", "302632", "302630", "302633", "302637", "302663", "302665",
                    "302667", "302674", "302675", "302680", "302692", "302693", "302696")

op1.2020.props = op1_df %>%
  filter(!(UpaFrom %in% missing_op1_upa) &
           !(UpaTo %in% missing_op1_upa)) %>% 
  group_by(UpaFrom) %>%
  mutate(prop.2020 = value/sum(value, na.rm=T)) %>%
  ungroup()

trans.sub.upa.long %>%
  filter(value != 0) %>%
  group_by(UpaFrom, date) %>%
  tally() %>%
  arrange(n)

op1.2017.data = trans.sub.upa.long %>%
  filter(!(date %in% c("2017-06-25", "2017-06-26", "2017-06-27", "2017-09-01", "2017-09-02", "2017-09-03"))) %>%
  group_by(weekend, UpaFrom, UpaTo) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>%
  ungroup()

op1.2017.props = op1.2017.data %>%
  filter(weekend == FALSE) %>%
  group_by(UpaFrom) %>%
  mutate(prop.2017= value/ sum(value, na.rm=T))  # eid dates already removed in weekend_df_op1

compare.props = op1.2020.props %>%
  left_join(., op1.2017.props, by = c("UpaFrom", "UpaTo"))%>%
  filter(UpaFrom != UpaTo)

compare.all.props = ggplot(compare.props) +
  geom_point(aes(x = prop.2020, y = prop.2017), size = 0.5) +
  geom_text(aes(x = Inf, y=Inf, hjust = 1, vjust = 1,
                label = paste0("r=",round(cor(prop.2020, prop.2017, use="complete.obs"),2)))) +
  # geom_abline() +
  export_theme +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) +
  labs(y = "Proportion of trips between upazilas in 2017",
       x = "Proportion of trips between upazilas in 2020")

ggsave(compare.all.props, file =  here("out", "supplement_figs", "s6.upa.prop.compare.png"), dpi = 300, width = 6, height = 4)

#0.73 correlation between proportion from upazila to all other upazilas, comparing 2017 vs. 2020 for Operator 1
cor(compare.props$prop.2020, compare.props$prop.2017, use="complete.obs") 

# Fig S8. date range figure ----
date_range = bind_rows(mobility_m_op1, mobility_m_op2, mobility_m_op3) %>%
  group_by(operator, date) %>%
  summarise(value = sum(value, na.rm=T)) %>%
  mutate(date_data = ifelse(value > 0 , 1, 0))

date_range_plt = ggplot(date_range) +
  geom_tile(aes(x = as.Date(date), y = factor(operator), fill = factor(date_data))) +
  scale_fill_manual(values = c("white","gray60")) +
  scale_y_discrete(limits=rev) + 
  scale_x_date(date_breaks  ="2 weeks")+
  export_theme +
  labs(y = "", x = "Date", fill = "") +
  theme(axis.text.x = element_text(angle = 45, hjust=1), 
        legend.position = "none")

ggsave(date_range_plt, file =  here("out", "supplement_figs", "s7.date.range.png"), dpi = 300, width = 7, height = 4)

# Fig S9. Effective population sizes ----
sum((upa_pop %>% filter(substr(upazila_code, 1,4)!=1042))$pop) #162717525 without Jhalakoti - missing for op1
pop.totals = upa_pop %>% left_join(., district_upa, by = "upazila_code")

pop.totals %>% group_by(District) %>% summarise(sum = sum(pop)) %>%arrange(sum) %>% View()
op1.sym = get_norm_sym_m(op1_m_district, pop_input = upa_pop, mobility_source = "Operator 1", sym = T)
op2.sym = get_norm_sym_m(op2_m_district, pop_input = upa_pop, mobility_source = "Operator 2", sym = T)
op3.sym = get_norm_sym_m(op3_m_district, pop_input = upa_pop, mobility_source = "Operator 3", sym = T)

bang.pop = data.frame(area = rep("Bangladesh", 4), 
                      value = c(sum(pop.totals$pop), sum(op1.sym), sum(op2.sym), sum(op3.sym)),
                      source = c("Actual pop", "Operator 1 effective pop", "Operator 2 effective pop", "Operator 3 effective pop"))
dhaka.pop =  data.frame(area = rep("Dhaka", 4), 
                        value = c(sum((pop.totals %>% filter(District == "Dhaka"))$pop), sum(op1.sym[, "Dhaka"]), sum(op2.sym[, "Dhaka"]), sum(op3.sym[, "Dhaka"])),
                        source = c("Actual pop", "Operator 1 effective pop", "Operator 2 effective pop", "Operator 3 effective pop"))
small.pop =  data.frame(area = rep("Bandarban", 4), 
                        value = c(sum((pop.totals %>% filter(District == "Bandarban"))$pop), sum(op1.sym[, "Bandarban"]), sum(op2.sym[, "Bandarban"]), sum(op3.sym[, "Bandarban"])),
                        source = c("Actual pop", "Operator 1 effective pop", "Operator 2 effective pop", "Operator 3 effective pop"))


compare.ep = bind_rows(bang.pop, dhaka.pop, small.pop) %>%
  mutate(area = factor(area, levels = c("Bangladesh", "Dhaka", "Bandarban"))) %>%
  ggplot(.) +
  geom_col(aes(x = source, y = value, fill = source)) +
  facet_wrap(~area, scales = "free") +
  export_theme +
  theme(legend.position = "none") +
  scale_fill_manual(values=met.brewer("Egypt", 5)) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(y = "Population", x = "", fill = "Source")

ggsave(compare.ep, file =  here("out", "supplement_figs", "s8.compare.pop.png"), dpi = 300, width = 10, height = 3)

# Fig S10. Calc proportion correlation for Dhaka to all other districts ----
prop.from.dhaka.2020 = weekend_df_op1 %>%
  filter(weekend == FALSE) %>%
  left_join(., district_upa, by = c("UpaTo" = "upazila_code")) %>%
  rename(DistrictTo = District) %>%
  left_join(., district_upa, by = c("UpaFrom" = "upazila_code")) %>%
  rename(DistrictFrom = District) %>%
  filter(DistrictFrom == "Dhaka" & !(UpaFrom %in% missing_op1_upa)) %>%
  filter(DistrictTo != "Dhaka") %>%
  group_by(UpaFrom, DistrictTo) %>%
  summarise(value = sum(value, na.rm=T)) %>%
  group_by(UpaFrom) %>%
  mutate(sum_value = sum(value, na.rm=T)) %>%
  mutate(`2020` = value/sum_value) %>%
  ungroup()

prop.from.dhaka.2017 = op1.2017.data %>%
  filter(weekend == FALSE) %>%
  left_join(., district_upa, by = c("UpaTo" = "upazila_code")) %>%
  rename(DistrictTo = District) %>%
  left_join(., district_upa, by = c("UpaFrom" = "upazila_code")) %>%
  rename(DistrictFrom = District) %>%
  filter(DistrictFrom == "Dhaka" & !(UpaFrom %in% missing_op1_upa)) %>%
  filter(DistrictTo != "Dhaka") %>%
  group_by(UpaFrom, DistrictTo) %>%
  summarise(value = sum(value, na.rm=T)) %>%
  group_by(UpaFrom) %>%
  mutate(sum_value = sum(value, na.rm=T)) %>%
  mutate(`2017` = value/sum_value) %>%
  ungroup()

nrow(prop.from.dhaka.2020)
nrow(prop.from.dhaka.2017)

combined.prop = prop.from.dhaka.2020 %>% 
  left_join(., prop.from.dhaka.2017, by = c("UpaFrom", "DistrictTo")) %>%
  select(UpaFrom, DistrictTo, `2020`, `2017`)

combined.prop.plt = combined.prop %>%
  pivot_longer(!c("UpaFrom", "DistrictTo"), values_to = "prop", names_to = "year") %>%
  ggplot() +
  geom_col(aes(x = UpaFrom, y = prop, fill = DistrictTo)) +
  facet_grid(year~.) +
  export_theme +
  scale_fill_manual(values=met.brewer("Egypt", 64)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.5, 'cm')) +
  labs(x = "Dhaka upazila",
       y = "Proportion of trips from Dhaka upazila \nto non-Dhaka Districts (Operator 1)")


ggsave(combined.prop.plt, file =  here("out", "supplement_figs", "s9.prop.from.dhaka.2017.plt.png"), dpi = 300, width = 8, height = 4)

# calculate r
round(cor(combined.prop$`2020`, combined.prop$`2017`), 2)
 

