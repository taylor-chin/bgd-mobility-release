library(here)
library(tidyverse)
library(patchwork)
library(tmap)
library(sp)
library(scales)
library(rgeos)
library(lubridate)
library(binom)

source(here("code", "utils.R")) # functions to make matrix

# pop data
# Meta pop
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

# load movement tile data and clean
movement_out  <- list.files(path = here("data", "fb-movement-data"), pattern = "*.csv", full.names = T) %>%
  lapply(read_csv) %>%
  bind_rows() %>%
  mutate(date.time = as.POSIXct(date_time, format = "%Y-%m-%d %H%M", tz = "GMT")) %>% #adjusts dates
  filter(country == "BD")

attributes(movement_out$date.time)$tzone <- "Asia/Dhaka"   
movement_out <- movement_out %>% 
  mutate(time = format(date.time, format="%H%:%M:%S"), 
         date = strftime(date.time, format="%Y-%m-%d"),
         day = strftime(date.time,'%A'))
nrow(movement_out) #9944678
movement_out %>% View()
head(movement_out)

#saveRDS(movement_out, file = here("data", "fb-movement-data", "fb-movement-cleaned.rds"))

movement_out = readRDS(file = here("data", "fb-movement-data", "fb-movement-cleaned.rds"))

# check start/end lon/lat columns - these are divisions 
movement_out %>% distinct(start_lon, start_lat, start_polygon_name) %>% arrange(start_polygon_name) %>% View()


# Map tile centroids to district polygons # -----
## . split geography column into origin, destination coordinates ====
latlon <- as.data.frame(stringr::str_split_fixed(as.character(movement_out$geometry)," ",5))

# lat/long of origin coordinates
s.tmp1 <- as.data.frame(cbind(substr(as.character(latlon$V2),2,nchar(as.character(latlon$V2))), 
                              substr(as.character(latlon$V3),1,(nchar(as.character(latlon$V3))-1))))
names(s.tmp1) <- c("lon", "lat")

nrow(s.tmp1)
s.tmp1 <- s.tmp1 %>% mutate(lat = as.numeric(as.character(lat)),
                            lon = as.numeric(as.character(lon))) %>% distinct()

nrow(s.tmp1) 
e.tmp1 <- as.data.frame(cbind(substr(as.character(latlon$V4),1,nchar(as.character(latlon$V4))),
                              substr(as.character(latlon$V5),1,(nchar(as.character(latlon$V5))-1))))

# repeat for destination tiles
names(e.tmp1) <- c("lon", "lat")
nrow(e.tmp1) 

e.tmp1 <- e.tmp1 %>% mutate(lat = as.numeric(as.character(lat)),
                            lon = as.numeric(as.character(lon))) %>% distinct()

nrow(e.tmp1) 

# df of unique combinations of origin/destination coordinates
tmp1 <- rbind(s.tmp1, e.tmp1) %>% distinct()

nrow(tmp1)

## . map lat/long to districts ====
bang_adm2 <- rgdal::readOGR(dsn = here("data","bgd_admbnda_adm2_bbs_20180410"), layer = "bgd_admbnda_adm2_bbs_20180410", stringsAsFactors = F)

bang_adm2@data <- bang_adm2@data %>% mutate(district_code = case_when(ADM2_PCODE == "4539" ~ "3039",
                                                                      ADM2_PCODE == "4561" ~ "3061",
                                                                      ADM2_PCODE == "4572" ~ "3072",
                                                                      ADM2_PCODE == "4589" ~ "3089",
                                                                      TRUE ~ as.character(ADM2_PCODE)))

coords_df <- tmp1 %>% distinct() %>%
  mutate(lon = as.numeric(lon), 
         lat = as.numeric(lat))

coordinates(coords_df) <-  ~ lon + lat
proj4string(coords_df) <- proj4string(bang_adm2)
loc <- over(coords_df, bang_adm2)
loc$lat <- coords_df$lat; loc$lon <- coords_df$lon; 

#loc <- loc %>% filter(Uni_Code != "NA") %>% distinct()
nrow(loc)
head(loc)

loclist <- loc %>% 
  select(ADM2_EN, district_code, lat, lon) %>% distinct() # distinct in redundant here

nrow(loclist)

plot(bang_adm2[bang_adm2$ADM2_EN == "Dhaka", ])
plot(coords_df, add = T, col = "red")

## . Sum baseline and crisis movements in and out for all lat/lons that fall in that district ====
# repeat steps to parse out lat/lon from geometry and add as columns to "df"
tmp2 = as.data.frame(cbind(substr(as.character(latlon$V2),2,nchar(as.character(latlon$V2))), 
                          substr(as.character(latlon$V3),1,(nchar(as.character(latlon$V3))-1)),
                          substr(as.character(latlon$V4),1,nchar(as.character(latlon$V4))),
                          substr(as.character(latlon$V5),1,(nchar(as.character(latlon$V5))-1))))



names(tmp2) <- c("s.lon", "s.lat", "e.lon", "e.lat")
tmp2 <- tmp2 %>% mutate(s.lat = as.numeric(as.character(s.lat)),
                      s.lon = as.numeric(as.character(s.lon)),
                      e.lat = as.numeric(as.character(e.lat)),
                      e.lon = as.numeric(as.character(e.lon)))

nrow(tmp2) # same as movement_out

df <- cbind(movement_out, tmp2)
head(df)

# add mapped district names
trans <- df %>% left_join(., loclist, by = c("s.lon" = "lon", "s.lat" = "lat")) %>%
  filter(district_code!="NA") %>%
  left_join(., loclist, by = c("e.lon" = "lon", "e.lat" = "lat")) %>% 
  filter(district_code.y!="NA")

nrow(df) - nrow(trans) #50951 rows dropped

# summarise across time points and aggregate to district
tm.all <- trans %>% group_by(ADM2_EN.x, ADM2_EN.y, date.time, date) %>% 
  # sum movement of tiles by start/end districts
  summarise(sum_baseline = sum(n_baseline, na.rm= T), sum_crisis = sum(n_crisis, na.rm=T)) %>%
  # average across 3 timepoints per day 
  group_by(ADM2_EN.x, ADM2_EN.y, date) %>%
  summarise(baseline = mean(sum_baseline, na.rm= T), current = mean(sum_crisis, na.rm=T)) %>%
  mutate(diff = current - baseline, prop.diff = current/baseline, 
         perc_change = (current-baseline)*100/baseline) %>%
  ungroup()

head(tm.all)
nrow(tm.all)

tm.all %>% filter(ADM2_EN.x == "Dhaka" & ADM2_EN.y == "Chittagong") %>%
  ggplot(.) +
  geom_line(aes(x = date, y = baseline, group = 1 ), colour = "red") +
  geom_line(aes(x = date, y = current, group = 1), colour = "black")
  

tm.all = readRDS(file = here("data", "fb-movement-data", "fb-movement-district-summary.rds"))

# check connectivity
# all districts have within-travel estimated
tm.all %>%
  group_by(ADM2_EN.x, ADM2_EN.y) %>%
  summarise(baseline = mean(baseline, na.rm=T), current = mean(current, na.rm=T)) %>%
  ungroup() %>%
  filter(ADM2_EN.x == ADM2_EN.y) %>%
  count(current > 0) # current > 0

tm.all %>%
  group_by(ADM2_EN.x, ADM2_EN.y) %>%
  summarise(baseline = mean(baseline, na.rm=T), current = mean(current, na.rm=T)) %>%
  ungroup() %>%
  filter(ADM2_EN.x != ADM2_EN.y & current > 0 ) %>% 
  group_by(ADM2_EN.x) %>%
  tally() %>%
  arrange(n)# current > 0

tm.all %>%
  summarise(max(date), min(date))

## Eid 2020 May 23-24 and July 31-Aug 1
eid.dates = c("2020-05-24","2020-05-25", "2020-05-26", "2020-07-31", "2020-08-01", "2020-08-02")

tm.all.mat.df = tm.all %>% filter(!(date %in% eid.dates)) %>%
  mutate(weekend = ifelse(wday(as.Date(date)) == 1, 1, 0) | ifelse(wday(as.Date(date)) == 7, 1, 0)) %>%
  rename(DistrictFrom = ADM2_EN.x, DistrictTo = ADM2_EN.y) %>%
  select(DistrictFrom, DistrictTo, date, current, weekend) %>% 
  group_by(DistrictFrom, DistrictTo, weekend) %>%
  summarise(value = mean(current, na.rm = TRUE)) %>%
  ungroup()

tm.all.mat.df %>%
  group_by(DistrictFrom, weekend) %>%
  count(value == 0) %>%
  arrange(desc(n))

# check connectivity again
# chuadanga, Lalmonirhat, meherpur
tm.all.mat.df %>%
  filter(DistrictFrom == "Chuadunga") %>% 
  filter(value == 0) %>% View()
  group_by(DistrictFrom, weekend) %>%
  tally() %>%
  arrange(desc(n))

# create matrix ----
## rerun below sections twice changing line 223 for weekend vs. weekday
raw_m = tm.all.mat.df %>% 
  filter(weekend == FALSE) %>% 
  select(DistrictFrom, DistrictTo, value) %>%
  pivot_wider(names_from = DistrictFrom, values_from = value) %>%
  remove_rownames %>% column_to_rownames(var="DistrictTo") %>%
  data.matrix()

raw_m = raw_m[order(rownames(raw_m)), order(colnames(raw_m))]
  
saveRDS(raw_m, here("out", "FB-analysis", "fb_weekday_m_raw_fulldates.rds"))

rowSums(raw_m != 0)
colSums(raw_m != 0)

op1_m_weekday_counts[rowSums(op1_m_weekday_counts != 0, na.rm=T)==544,]
rowSums(op1_m_weekday_counts == 0, na.rm=T)

## without uncertainty ====
raw_m[is.na(raw_m)] <- 0 

prop_m = prop.table(raw_m, 2)
colSums(prop_m) # check 1s
prop_m = prop_m[order(rownames(prop_m)), order(colnames(prop_m))]
prop_m[is.nan(prop_m)] <- 0 # 

# normalize to District pops
pop_order = pop %>% 
  left_join(., district_upa %>% select(District, upazila_code), by = "upazila_code") %>%
  group_by(District) %>%
  summarise(pop = sum(pop)) %>%
  arrange(match(District, colnames(prop_m)))

m_norm = prop_m%*%diag(pop_order$pop)
colnames(m_norm) = rownames(m_norm)
print(colSums(m_norm, na.rm=T)) # check pop totals


sym_norm_m = (m_norm + t(m_norm))/2

saveRDS(sym_norm_m, here("out", "FB-analysis", "fb_weekend_m_withNA_symmetric_fulldates.rds"))


## using beta distribution ====
raw_m = raw_m[order(rownames(raw_m)), order(colnames(raw_m))]

raw_m[is.na(raw_m)] <- 0 
num = c(raw_m)
denom = rep(unname(colSums(raw_m, na.rm=T)), each = nrow(raw_m))

# doesn't matter conf.level or "type" because only using shape parameters
bayes.results = binom.bayes(x = num, n = denom, type = "central", 
                            conf.level = 0.8) 

alpha.mat = matrix(c(bayes.results$shape1), nrow = nrow(raw_m))
beta.mat = matrix(c(bayes.results$shape2), nrow = nrow(raw_m))

colnames(alpha.mat) = colnames(raw_m)
rownames(alpha.mat) = rownames(raw_m)
colnames(beta.mat) = colnames(raw_m)
rownames(beta.mat) = rownames(raw_m)

shape.mat = list(alpha = alpha.mat, beta = beta.mat)

saveRDS(shape.mat, here("out", "FB-analysis","fb_weekend_beta_param_rawm.rds"))
