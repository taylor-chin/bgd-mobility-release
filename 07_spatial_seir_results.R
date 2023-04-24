library(here)
library(tidyverse)
library(cowplot)
library(stringr) #wrap labels
library(ggh4x) # facetgrid2 for independent y axis
library(MetBrewer)

# load helper functions in utils
source(here("code","utils.R"))

# load results
all.results.bind = readRDS(here("out", "metapop_model", "cluster_results", "results.district.beta.avgsim.vF.rds"))
sumbgd = readRDS(here("out", "metapop_model", "cluster_results", "results.district.beta.sumbgd.vF.rds"))

# join to geo data for heatmap 
bang2_adm <- readOGR(dsn = here("data","bgd_admbnda_adm2_bbs_20180410"), layer = "bgd_admbnda_adm2_bbs_20180410", stringsAsFactors = F)
dist_div = bang2_adm@data %>% select(ADM2_EN, ADM1_EN)
rm(bang2_adm)
all.results.bind = all.results.bind %>% left_join(., dist_div, by = c("geo_code" = "ADM2_EN"))


# SEIR results ----
# FIG 3 ----
## . panel A: epi curves ====
p1 = sumbgd %>%
  filter(r0 == 1.3) %>%
  group_by(time, source, seed_city) %>%
  summarise(mean = mean(value, na.rm=TRUE),
            lower = quantile(value, probs = 0.025, na.rm=TRUE),
            upper = quantile(value, probs = 0.975, na.rm=TRUE)) %>%
  mutate(source = factor(source, levels = c("Gravity model", "Operator 1", "Operator 2", "Operator 3", "Meta"))) %>%
  ggplot(.) +
  geom_line(aes(x = time, y = mean, colour = as.factor(source)), size = 0.65) +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper, fill = as.factor(source)), 
              alpha = 0.3) +
  facet_wrap(~seed_city, scales = "free") +
  theme_classic(base_size = 10) +
  theme(strip.text=element_text(size=8),
        strip.background=element_rect(fill="#f0f0f0")) +
  scale_colour_manual(values=met.brewer("Java", 5)) +
  scale_fill_manual(values=met.brewer("Java", 5)) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(limits = c(50,525)) +
  labs(y = "Incidence", x = "Time (days)", colour = "Mobility matrix source", 
       fill = "Mobility matrix source")

## . panel B: peak and take off ====
peak_df = all.results.bind %>%
  filter(state == "obs" & r0 == 1.3) %>%
  group_by(geo_code, seed_city, source, time) %>%
  summarise(sum = sum(mean)) %>%
  group_by(source, geo_code, seed_city) %>%
  top_n(1, sum) %>% 
  select(time, source, geo_code, seed_city, sum) %>%
  mutate(summary = "peak_time") %>%
  filter( sum != 0)

takeoff_df = all.results.bind %>%
  filter(state == "obs" & r0 == 1.3) %>%
  group_by(geo_code, seed_city, source, time) %>%
  summarise(sum = sum(mean)) %>%
  filter(sum >= 50) %>%
  group_by(source, geo_code, seed_city) %>% 
  slice(1) %>% 
  select(time, source, geo_code, seed_city, sum) %>%
  mutate(summary = "takeoff_time")

summary_df = bind_rows(peak_df, takeoff_df)

# order y-axis by distance from seed city
# created in utils.R
distanceBetweenDistricts = readRDS(here("data", "misc", "distanceBetweenDistricts.rds"))

new_order = distanceBetweenDistricts %>% 
  filter(District.x %in% c("Dhaka", "Chittagong", "Panchagarh")) %>%
  group_by(District.x) %>%
  do(data_frame(order_geo_code=levels(reorder(interaction(.$District.x, .$District.y, drop=TRUE), .$distance_mi)))) %>% 
  pull(order_geo_code)
  
p2 = summary_df %>%
  mutate(source = factor(source, levels = c("Gravity model", "Operator 1", "Operator 2", "Operator 3", "Meta"))) %>%
  mutate(order_geo_code=factor(interaction(seed_city, geo_code), levels=new_order)) %>%
  ggplot(.) +
  geom_point(aes(x = time, 
                 y = order_geo_code, colour = summary), size = 0.25) + 
  ggh4x::facet_grid2(source~seed_city, 
             labeller = label_wrap_gen(width=12), 
             scales = "free_y", independent = "y") +
  scale_x_continuous(limits = c(50,525)) +
  scale_y_discrete(breaks= new_order, labels=gsub("^.*\\.", "", new_order)) +
  scale_colour_discrete(labels=c('Outbreak peak', 'Outbreak take off')) +
  theme_bw(base_size = 10) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = 'white'),
        panel.grid.minor = element_line(color = "white"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text=element_text(size=8),
        strip.background=element_rect(fill="#f0f0f0")) +
  labs(y = "District", x = "Time (days)", colour = "")

# ggarrange all together ----
combined_plt = plot_grid(p1, p2, labels = c('A', 'B'), nrow = 2,
                         rel_heights = c(0.4, 1), align = 'hv', axis = 'lr')

ggsave(combined_plt, file =  here("out", "metapop_model", "figures", "fig3_curves_dot_vF.png"), 
       width = 10, height = 6.5)


## numbers used in manuscript ====
all.results.bind %>%
  filter(state == "obs" & r0 == 1.3) %>%
  group_by(time, source, seed_city) %>% #geo_code
  # sum across districts 
  summarise(mean = sum(mean), 
            lower = sum(lower), 
            upper = sum(upper)) %>% 
  filter(seed_city == "Chittagong") %>%
  group_by(source) %>%
  slice(which.max(mean))

peak_df %>% filter(source == "Gravity model") %>% 
  group_by(seed_city) %>%
  summarise(min(time), max(time))

peak_df %>% ungroup() %>%
  filter(source == "Gravity model") %>% 
  filter(seed_city == "Panchagarh") %>%
  slice(which.min(time), which.max(time))

peak_df %>% filter(source == "Operator 1") %>% 
  group_by(seed_city) %>%
  summarise(min(time), max(time), mean(time))

peak_df %>% filter(source == "Operator 2") %>% 
  group_by(seed_city) %>%
  summarise(min(time), max(time), mean(time))

# FIG 4 time to intro categorical plot ----
cat_intro = all.results.bind %>%
  filter(state == "obs" & r0 == 1.3) %>%
  group_by(geo_code, source, seed_city) %>%
  mutate(cumsum = cumsum(mean)) %>% 
  mutate(time_intro_cat = case_when(cumsum >= 10 & time == 30 ~ "First 30 days",
                                    cumsum >= 10 & time == 60 ~ "31-60 days",
                                    cumsum >= 10 & time == max(time) ~ "61+ days",
                                    cumsum < 10 & time == max(time) ~ "Missing Data")) %>% 
  filter(!is.na(time_intro_cat)) %>% 
  group_by(geo_code, source, seed_city) %>%
  # slice first instance to get earliest time
  slice(1) %>%
  mutate(time_intro_cat = factor(time_intro_cat, levels = c("First 30 days", "31-60 days", "61+ days", "Missing Data"))) %>%
  mutate(source = factor(source, levels = c("Gravity model", "Operator 1", "Operator 2", "Operator 3", "Meta"))) #

colors = c(met.brewer("Hiroshige", n = 10)[c(1,3,5)], "#C0C0C0")

tmp_map = cat_intro %>%
  left_join(map_district,.,
            #  by = c("upazila_code" = "geo_code")) %>%
            by=c("ADM2_EN" = "geo_code")) %>%     # change this depending on geographic area
  st_as_sf %>%
  st_transform(crs = 4326)

time_plt= ggplot() +
  geom_sf(data = map_district, size = 0.05) + 
  geom_sf(data = tmp_map %>% mutate(source = factor(source, levels = c("Gravity model", "Operator 1", "Operator 2", "Operator 3", "Meta"))), #
          aes(fill = factor(time_intro_cat)), size = 0.05) +
  scale_fill_manual(values = colors) +
  facet_grid(seed_city~source) +
  theme_void() +
  labs(fill = "Time to introduction \n(at least 10 cases)") +
  theme(legend.text  = element_text(size = 8), 
        legend.position = "bottom",
        strip.background=element_rect(fill="#f0f0f0")) 

ggsave(time_plt, file =  here("out", "metapop_model", "figures", "fig4_cat_intro_maps_vF.png"), 
       width = 9, height = 7)

# FIG 5 summary stats by r0 ----
## load in results from cluster
final.size = readRDS(here("out", "metapop_model", "cluster_results", "results.district.beta.finalsize.vF.rds"))
num.dist = readRDS(here("out", "metapop_model", "cluster_results", "results.district.beta.numdist30.vF.rds"))

p1 = final.size %>%
  mutate(source = factor(source, levels = c("Gravity model", "Operator 1", "Operator 2", "Operator 3", "Meta"))) %>% 
  mutate(seed_city = factor(seed_city, levels = c("Dhaka", "Chittagong", "Panchagarh"))) %>%
  ggplot(.) +
  geom_bar(aes(x = factor(r0), y = final_size_prop, fill = source), position=position_dodge(0.95), stat = "identity") +
  geom_errorbar(aes(x = factor(r0), ymin=final_size_prop_lower, ymax=final_size_prop_upper, fill = source),
                width=.5, position=position_dodge(0.95)) +
  facet_grid(~seed_city) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12)) + 
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_fill_manual(values=met.brewer("Java", 5)) +
  theme_classic(base_size = 10) +
  theme( strip.text=element_text(size=8),
         strip.background=element_rect(fill="#f0f0f0")) +
  labs(y = "Final size epidemic \n(proportion)", 
       x = expression(R[0]), 
       fill = "Mobility source")

## . panel B: number of districts spread in first 30 days ====
p2 = num.dist %>%
  mutate(source = factor(source, levels = c("Gravity model", "Operator 1", "Operator 2", "Operator 3", "Meta"))) %>%  #
  mutate(seed_city = factor(seed_city, levels = c("Dhaka", "Chittagong", "Panchagarh"))) %>%
  ggplot(.) +
  geom_bar(aes(x = factor(r0), y = num_dist, fill = source), position=position_dodge(0.95), stat = "identity") +
  geom_errorbar(aes(x = factor(r0), ymin=num_dist_lower, ymax=num_dist_upper, fill = source),
                width=.5, position=position_dodge(0.95)) +
  facet_grid(~seed_city) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12)) + 
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_fill_manual(values=met.brewer("Java", 5)) +
  theme_classic(base_size = 10) +
  theme( strip.text=element_text(size=8),
         strip.background=element_rect(fill="#f0f0f0")) +
  labs(y = "Number of districts with \nat least one case within 30 days", 
       x = expression(R[0]), fill = "Mobility source")

combined_plot = plot_grid(p1, p2, labels = c('A', 'B'), nrow = 2,
                          align = 'h')


ggsave(combined_plot, file = here("out", "metapop_model", "figures", "fig5_r0_summary_plt_30days.png"), 
       width = 8, height = 5)



# figure used in manuscript
num.dist %>% 
  filter(source == "Operator 2" | source == "Meta") %>% view()


