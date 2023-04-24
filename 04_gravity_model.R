library(rgdal)
library(rgeos)
library(geosphere)

source(here("code","utils.R")) # distance between upas and districts made here
# need pop, district_upa

distanceBetweenDistricts = readRDS(here("data", "misc", "distanceBetweenDistricts.rds"))

district_pop = pop %>% # from utils
  left_join(., district_upa %>% select(District, upazila_code), by = "upazila_code") %>%
  group_by(District) %>%
  summarise(pop = sum(pop)) %>%
  arrange(match(District, colnames(matrix)))

distance.df = distanceBetweenDistricts %>%
  left_join(., district_pop, by = c("District.x" = "District")) %>%
  rename(pop.origin = pop) %>%
  left_join(., district_pop, by = c("District.y" = "District")) %>%
  rename(pop.dest = pop)
  
head(distanceBetweenDistricts)

#### calculate gravity model with parameters from Wesolowski 2015 ####
#https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004267

# full data set parameters: k = -20.61, alpha/pop from = 1.22, beta/pop to = 1.22, gamma = -2.05
grav.df <- distance.df %>%
  mutate(distance_mi = as.numeric(as.character(distance_mi))) %>%
  mutate(trav = exp(-20.61)*(pop.dest^1.22 * pop.origin^1.22)/distance_mi^2.05) %>%
  arrange(-distance_mi)

grav.mat = grav.df %>%
  select(District.x, District.y, trav) %>%
  pivot_wider(names_from = District.x, values_from = trav) %>%
  remove_rownames %>% column_to_rownames(var="District.y") %>%
  data.matrix()

grav.mat.order = grav.mat[order(rownames(grav.mat)), order(colnames(grav.mat))]
grav.mat.order[is.infinite(grav.mat.order)] <- 0
colSums(grav.mat.order)

saveRDS(grav.mat.order, here("out", "metapop_model", "grav_mat.w.lit.params.rds"))

# plot
as.data.frame(grav.mat.order) %>%
  rownames_to_column(var = "District.y") %>%
  pivot_longer(!District.y, names_to = "District.x", values_to = "trav") %>% 
  mutate(log_value = log(trav)) %>%
  ggplot(.) +geom_tile(aes(x = District.x, y = District.y, fill = log_value)) +
  theme_bw(base_size = 8) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

# upazila level
distanceBetweenUpazilas = readRDS(here("data", "misc", "distanceBetweenUpazilas.rds"))

distance.df.upa = distanceBetweenUpazilas %>%
  left_join(., pop, by = c("ADM3_PCODE_home" = "upazila_code")) %>%
  rename(pop.origin = pop) %>%
  left_join(., pop, by = c("ADM3_PCODE_dest" = "upazila_code")) %>%
  rename(pop.dest = pop)

grav.df.upa <- distance.df.upa %>%
  mutate(distance_mi = as.numeric(as.character(distance_mi))) %>%
  mutate(trav = exp(-20.61)*(pop.dest^1.22 * pop.origin^1.22)/distance_mi^2.05) %>%
  arrange(-distance_mi)

grav.mat.upa = grav.df.upa %>%
  select(ADM3_PCODE_home, ADM3_PCODE_dest, trav) %>%
  pivot_wider(names_from = ADM3_PCODE_home, values_from = trav) %>%
  remove_rownames %>% column_to_rownames(var="ADM3_PCODE_dest") %>%
  data.matrix()

grav.mat.order.upa = grav.mat.upa[order(rownames(grav.mat.upa)), order(colnames(grav.mat.upa))]
grav.mat.order.upa[is.infinite(grav.mat.order.upa)] <- 0
colSums(grav.mat.order.upa)

saveRDS(grav.mat.order.upa, here("out", "metapop_model", "grav_mat.upa.w.lit.params.rds"))

