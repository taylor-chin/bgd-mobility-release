# 06_descrip_compare.R 
#
# Descriptive comparisons of mobility patterns estimated from 
# 3 CDR operators, Meta Data for Good, gravity model
# on district level
# use averages across weekdays in 2020 from each data source

# Load dependencies ----
library(here)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(scales)
library(svglite)

source(here("code", "utils.R")) #matrix functions in here

# read in data ----
## final matrices (after normalization, symmetric step) ====
## created in 01_process_cdr.R
op1_m = readRDS(here("out", "2020_cdr", "matrices", "op1_m_weekday_district_withNA_symmetric_dhakamult.rds"))
op2_m = readRDS(here("out", "2020_cdr", "matrices", "op2_m_weekday_district_withNA_symmetric_dhakamult.rds"))
op3_m = readRDS(here("out", "2020_cdr", "matrices", "op3_m_weekday_district_withNA_symmetric_dhakamult.rds"))

## read in district-level Meta data from fb_movement.R:
## created in 03_fb_movement.R
fb_m = readRDS(here("out", "FB-analysis", "fb_weekday_m_withNA_symmetric_fulldates.rds"))

# read in gravity model from 04_gravity_model.R:
## district level
grav_lit_m = readRDS(here("out", "metapop_model", "grav_mat.w.lit.params.rds"))

# to order my latitude/longitude,
# read in df with centroid lat/long of districts (created in utils.R)
cent_district_df = readRDS(here("data", "misc", "cent_district_df.rds")) %>%
  arrange(lat) # used in get_df_plot function

op1_df = get_df_plot(op1_m) %>% mutate(source = "Operator 1")
op2_df = get_df_plot(op2_m) %>% mutate(source = "Operator 2")
op3_df = get_df_plot(op3_m) %>% mutate(source = "Operator 3")
fb_df = get_df_plot(fb_m) %>% mutate(source = "Meta")
grav_df = get_df_plot(grav_lit_m) %>% mutate(source = "Gravity Model")

compare_source_df = bind_rows(op1_df, op2_df, op3_df, fb_df, 
                                    grav_df) %>%
  mutate(source = factor(source, levels = c("Operator 1", "Operator 2", "Operator 3", "Meta", 
                                           "Gravity Model")))
  

# Add distances between districts
distanceBetweenDistricts = readRDS(here("data", "misc", "distanceBetweenDistricts.rds"))
compare_source_df = compare_source_df %>%
  left_join(., distanceBetweenDistricts, by = c("DistrictTo" = "District.x", "DistrictFrom" = "District.y"))


# Fig1A. Map of normalized subscribers ----
subsc <- readRDS(here("out","2020_cdr", "subsc_2020.rds"))

missing_op1_upa = c("302602", "302605", "302606", "302609", "302610", "302611", "302624",
                        "302629", "302632", "302630", "302633", "302637", "302663", "302665",
                        "302667", "302674", "302675", "302680", "302692", "302693", "302696")

pop_subsc = subsc %>% select(upazila_code, operator, average) %>%
  filter(upazila_code != "557384") %>%
  left_join(., district_upa, by = "upazila_code") %>%
  group_by(District, operator) %>%
  summarise(sum = sum(average)) %>%
  # district level
  mutate(sum_withNA = case_when(sum == 0 ~ NA_real_,
                                TRUE ~ sum)) %>%
  mutate(operator = factor(operator, levels = c("Operator 1", "Operator 2", "Operator 3"))) %>%
  ungroup() %>%
  group_by(operator) %>%
  mutate(norm_subsc = sum_withNA/max(sum_withNA, na.rm=T)) %>%
  ungroup() %>%
  select(District, operator, norm_subsc)

# FB population data
fb_pop_clean <- readRDS(file = here("data", "fb-pop-data", "fb-pop-cleaned-mapped.rds"))

fb_pop_sub = fb_pop_clean %>% 
  select(date_time, n_baseline, n_crisis, date.time, time, date, ADM2_EN, ADM2_PCODE) %>%
  mutate(n_baseline = as.numeric(n_baseline), 
         n_crisis = as.numeric(n_crisis))

fb_pop_summary = fb_pop_sub %>%
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

eid.dates = c("2020-05-24","2020-05-25", "2020-05-26", "2020-07-31", "2020-08-01", "2020-08-02")

fb_baseline_pop = fb_pop_summary %>%
  filter(!is.na(ADM2_EN) & !(date %in% eid.dates)) %>%
  group_by(ADM2_EN) %>%
  summarise(pop_baseline = mean(n_baseline, na.rm=T)) %>%
  mutate(norm_subsc = pop_baseline/max(pop_baseline, na.rm=T)) %>%
  ungroup() %>% 
  mutate(operator = "Meta") %>%
  rename(District = ADM2_EN) %>%
  select(-pop_baseline)

subsc_plt = bind_rows(pop_subsc, fb_baseline_pop) %>%
  mutate(operator = factor(operator, levels = c("Operator 1", "Operator 2", "Operator 3", "Meta")))

tmp_map = subsc_plt %>%
  left_join(map_district, ., 
            by=c("ADM2_EN" = "District")) %>% 
  st_as_sf %>%
  st_transform(crs = 4326)

# Fig1B. Top routes ----- 
pop_district = compare_source_df %>%
  group_by(source, DistrictFrom)%>%
  summarise(origin_pop = sum(value, na.rm=T))

rank_df = compare_source_df %>%
  left_join(., pop_district, by = c("DistrictFrom", "source")) %>% 
  group_by(source, DistrictTo, DistrictFrom) %>%
  summarise(norm_value = value/origin_pop) %>%
  filter(DistrictTo != DistrictFrom) %>%
  group_by(source) %>%
  mutate(rank = rank(desc(norm_value))) %>%
  ungroup() %>%
  arrange(source, rank) %>%
  filter(rank <= 15) %>%
  rename(operator = source)

# look at overlapping routes
#rank_df %>% mutate(route_name = paste0(DistrictTo, "-", DistrictFrom)) %>%
#  filter(operator != "Gravity Model") %>%
#  count(route_name) %>% arrange(desc(n))

ranks = rank_df %>%
  # add long/lat from centroid df
  left_join(., cent_district_df %>% mutate(district_code = as.numeric(as.character(district_code))), 
            by = c("DistrictFrom" ="District")) %>%
  mutate(x = as.numeric(as.character(long)),
         y = as.numeric(as.character(lat))) %>%
  select(-c(long, lat)) %>%
  left_join(., cent_district_df %>% mutate(district_code = as.numeric(as.character(district_code))), 
            by = c("DistrictTo" ="District")) %>%
  mutate(xend = as.numeric(as.character(long)),
         yend = as.numeric(as.character(lat))) %>%
  select(-c(long, lat))

# overlay A and B maps
combined.plt = ggplot() +
  geom_sf(data = tmp_map, aes(fill = norm_subsc),lwd = 0.3) +
  geom_curve(aes(x = x, y = y, xend = xend, yend = yend,     # draw edges as arcs
                 color = as.factor(rank)),
             data =ranks %>% filter(!(operator %in% c("Gravity Model"))), 
             curvature = 0.2, size = 1) +
  facet_wrap(~operator, nrow = 1) +
  scale_fill_distiller(palette = "Blues", na.value = "grey90", direction =+1) +
  scale_colour_viridis_d(option = "rocket") +
  theme_void(base_size = 12) +
  labs(colour = "Top 15 relative routes", 
       fill = "Number of subscribers \n(normalized)") +
  theme(legend.title=element_text(size=10), 
        legend.text  = element_text(size = 10), 
        strip.text.x = element_text(margin = margin( b = 1, t = 0)),
        strip.background=element_rect(fill="#f0f0f0"),
        legend.position = "right",
        plot.margin = unit(c(0, 0, -1, 0), "cm")) +
  guides(color=guide_legend(ncol=2))

ggsave(combined.plt, file =  here("out", "descriptive", "fig1.descrip_binary.png"), 
       width = 9, height = 7)


# Fig 2. District level matrices ----
min_value = floor(summary((compare_source_df %>% filter(!is.infinite(log_value)))$log_value)[1])
max_value = ceiling(summary((compare_source_df %>% filter(!is.infinite(log_value)))$log_value)[6])

# proportion of travel within district vs between to add to matrix
within_district_prop = compare_source_df %>%
  mutate(same_district = ifelse(DistrictFrom == DistrictTo, 1, 0)) %>%
  group_by(source, DistrictFrom, same_district) %>%
  summarise(value = sum(value, na.rm=T)) %>%
  mutate(prop = value/sum(value, na.rm=T)) %>%
  filter(same_district == 1) %>%
  group_by(source) %>%
  summarise(summary_prop = round(mean(prop, na.rm=T), 2)) %>%
  # make gravity model empty for now
  mutate(summary_prop = ifelse(source == "Gravity Model", "", summary_prop))

compare_source_df_withprop = compare_source_df %>%
  left_join(., within_district_prop, by = "source")

heat = ggplot(compare_source_df_withprop) +
  geom_tile(aes(x = DistrictFrom, y = DistrictTo, fill = log_value_NA)) +
  geom_text(aes(x = -Inf, y = Inf, hjust = -0.2, vjust = 1.2, label = summary_prop)) +
  scale_fill_gradientn(colors = c(brewer.pal(9, "Blues"),"#999999"), na.value = "black",
                       limit = c(min_value,max_value+2),oob=scales::squish) +
  facet_wrap(~source, scales = "free") +
  theme_classic(base_size = 12) + 
  theme(legend.title=element_text(size=10), 
        axis.line=element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = "Origin District", y = "Destination District", fill = "Average log number of trips")

ggsave(heat, file =  here("out", "descriptive", "fig2.districtmat_withdiag_v2.png"), width = 12, height = 6.5)

# IQR values in manuscript ====
##  IQR values at upazila-level are in 08_supp_figs.R
compare_source_df %>%
  mutate(same_district = ifelse(DistrictFrom == DistrictTo, 1, 0)) %>%
  group_by(source, DistrictFrom, same_district) %>%
  summarise(value = sum(value, na.rm=T)) %>%
  mutate(prop = value/sum(value, na.rm=T)) %>%
  filter(same_district == 1) %>%
  group_by(source) %>%
  summarise(mean = mean(prop, na.rm=T),
            median = median(prop, na.rm=T),
            iqr = IQR(prop, na.rm=T))
