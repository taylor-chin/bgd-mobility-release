# dependencies
library(here)
library(tidyverse)
library(binom)
library(sf) #st_transform
library(Hmisc) #cut2
library(sp) #SpatialPointsDataFrame
library(rgeos) #gCentroid
library(geosphere) #distVincentyEllipsoid
library(rgdal) #readOGR

source(here("code", "mat_utils.R")) 

## Consistent plot theme ----
export_theme <- theme_classic() + 
  theme(
    ## Axis text and titles
    axis.text.x = element_text(size=8),
    axis.text.y=element_text(size=8),
    axis.title.x=element_text(size=10),
    axis.title.y=element_text(size=10),
    ## Axis lines
    axis.line = element_line(),
    axis.ticks = element_line(),
    ## Legends
    legend.title=element_text(size=10),
    legend.text=element_text(size=8),
    legend.key.size= unit(0.5, "cm"),
    legend.margin = margin(0,0,0,0, "cm"),
    ## Strips for facet_wrap
    strip.text=element_text(size=8),
    strip.background=element_rect(fill="#f0f0f0"),
    ## Tags
    plot.tag=element_text(size=10,face="bold",hjust=0,vjust=-3),
    plot.title=element_text(size=8,hjust=0.5,face="bold",vjust=-3)
  )

## Demographic data processing ----
#### upazila code - upazila name mapping ####
upa <- read.csv(file = here("data","popcount_by_union.csv"), header=T) %>%
  mutate(upazila_code = as.character(Upa_Code)) %>%
  dplyr::select(upazila_code, Upazila) %>% distinct()

#### district-upazila mapping ###
district_upa <- read.csv(file = here("data","popcount_by_union.csv"), header=T) %>%
  mutate(upazila_code = as.character(Upa_Code), Dis_Code = as.character(Dis_Code), District = as.character(District)) %>%
  dplyr::select(upazila_code, Upazila, Dis_Code, District) %>% distinct()

district_upa$District[district_upa$District == "Cox'S Bazar"] <- "Cox's Bazar"

# updated this to remove the double entry for upa code "303334" in the unique Upazila name
district_upa <- district_upa %>% filter(Upazila != "Kaliganj Paurashava") # two entries for this upa_code "303334"

#### get population data ####
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
  left_join(., upa, by = "upazila_code") %>%
  # duplicate entry
  filter(Upazila != "Kaliganj Paurashava")


## Get maps ----
# shapefile source: https://data.humdata.org/dataset/administrative-boundaries-of-bangladesh-as-of-2015
map <- list.files(path = here("data","bgd_admbnda_adm3_bbs_20180410/"), pattern = "*?\\.shp$", full.names = T) %>%
  lapply(st_read) %>%
  purrr::reduce(st_combine) %>% #for multiple map areas
  st_transform(crs = 4326) %>%   #projects to WGS84
  mutate(upazila_code = case_when(substr(ADM3_PCODE, 1, 4) == "4539" ~ paste0("3039", substr(ADM3_PCODE, 5, 7)),
                                  substr(ADM3_PCODE, 1, 4) == "4561" ~ paste0("3061", substr(ADM3_PCODE, 5, 7)),
                                  substr(ADM3_PCODE, 1, 4) == "4572" ~ paste0("3072", substr(ADM3_PCODE, 5, 7)),
                                  substr(ADM3_PCODE, 1, 4) == "4589" ~ paste0("3089", substr(ADM3_PCODE, 5, 7)),
                                  TRUE ~ as.character(ADM3_PCODE))) %>%
  mutate(area_km2 = units::set_units(st_area(.), km^2)) %>%
  left_join(., pop %>% dplyr::select(-"Upazila"), by = "upazila_code") %>%
  mutate(pop_density = as.numeric(pop/area_km2)) %>%
  mutate(pop_density_topcoded = case_when(pop_density > 25000 ~ 25000,
                                          TRUE ~ pop_density)) %>%
  mutate(pop_density_quantile = cut2(pop_density, g=5)) %>%
  distinct(ADM3_PCODE, ADM2_PCODE, ADM1_PCODE, date, upazila_code, .keep_all = TRUE) # remove duplicate row for upa code  "303334"

map_simple <- st_simplify(map, dTolerance = .0015) %>%
  st_transform(crs = 4326) %>%
  rename(map_date = date)

latlong = map %>% st_centroid() %>% st_coordinates() %>% as.data.frame() %>% 
  bind_cols(upazila_name = map_simple$ADM3_EN, upazila_code = map_simple$ADM3_PCODE) %>%
  rename(lat = X,
         long = Y) %>%
  arrange(lat)

map_district <- list.files(path = here("data","bgd_admbnda_adm2_bbs_20180410/"), pattern = "*?\\.shp$", full.names = T) %>%
  lapply(st_read) %>%
  purrr::reduce(st_combine) %>% #for multiple map areas
  st_transform(crs = 4326) 


# convert to long df for plotting and order by latitude
get_df_plot = function(df){
  df_plot = df %>% as.data.frame() %>%
    rownames_to_column(var = "DistrictTo") %>%
    pivot_longer(!DistrictTo, names_to = "DistrictFrom", values_to = "value") %>% 
    mutate(log_value = log(value)) %>%
    mutate(log_value_NA = case_when(is.infinite(log_value) ~ 100000, # random high value for plottin
                                    TRUE ~ log_value)) %>%
    mutate(DistrictTo = factor(DistrictTo, levels = (cent_district_df %>% pull(District))), 
           DistrictFrom = factor(DistrictFrom, levels = (cent_district_df %>% pull(District))))
  
  return(df_plot)
}

## Distance between districts and upazilas
bang_adm2 <- rgdal::readOGR(dsn = here("data","bgd_admbnda_adm2_bbs_20180410"), layer = "bgd_admbnda_adm2_bbs_20180410", stringsAsFactors = F)

bang_adm2@data <- bang_adm2@data %>% mutate(district_code = case_when(ADM2_PCODE == "4539" ~ "3039",
                                                                      ADM2_PCODE == "4561" ~ "3061",
                                                                      ADM2_PCODE == "4572" ~ "3072",
                                                                      ADM2_PCODE == "4589" ~ "3089",
                                                                      TRUE ~ as.character(ADM2_PCODE)))

cent_district <- SpatialPointsDataFrame(gCentroid(bang_adm2, byid=TRUE), bang_adm2@data, match.ID=FALSE)
cent_district_df <- as.data.frame(cbind(long = cent_district$x, lat = cent_district$y, district_code=bang_adm2@data$district_cod,
                                        District=bang_adm2@data$ADM2_EN))

distanceBetweenDistricts = expand.grid(ADM2_PCODE_home = bang_adm2@data$district_code, 
                                       ADM2_PCODE_dest = bang_adm2@data$district_code) %>%
  left_join(., cent_district_df, by = c("ADM2_PCODE_home" = "district_code")) %>%
  rename(long_home = long, 
         lat_home = lat) %>%
  left_join(., cent_district_df, by = c("ADM2_PCODE_dest" = "district_code")) %>%
  rename(long_dest = long, 
         lat_dest = lat) %>% 
  mutate(long_home = as.numeric(as.character(long_home)), 
         lat_home = as.numeric(as.character(lat_home)),
         long_dest = as.numeric(as.character(long_dest)),
         lat_dest = as.numeric(as.character(lat_dest)))

distance = distVincentyEllipsoid(distanceBetweenDistricts[,3:4],
                                 distanceBetweenDistricts[,6:7])

distanceBetweenDistricts = distanceBetweenDistricts %>% mutate(distance_mi = distance*0.000621371) %>%
  mutate(ADM2_PCODE_home = as.numeric(as.character(ADM2_PCODE_home)),
         ADM2_PCODE_dest = as.numeric(as.character(ADM2_PCODE_dest))) 

#saveRDS(distanceBetweenDistricts, here("data", "misc", "distanceBetweenDistricts.rds"))
#saveRDS(cent_district_df, here("data", "misc", "cent_district_df.rds"))

# upazila level
bang <- readOGR(dsn = here("data","bgd_admbnda_adm3_bbs_20180410"), layer = "bgd_admbnda_adm3_bbs_20180410", stringsAsFactors = F)
bang@data <- bang@data %>% mutate(upazila_code = case_when(substr(ADM3_PCODE, 1, 4) == "4539" ~ paste0("3039", substr(ADM3_PCODE, 5, 7)),
                                                           substr(ADM3_PCODE, 1, 4) == "4561" ~ paste0("3061", substr(ADM3_PCODE, 5, 7)),
                                                           substr(ADM3_PCODE, 1, 4) == "4572" ~ paste0("3072", substr(ADM3_PCODE, 5, 7)),
                                                           substr(ADM3_PCODE, 1, 4) == "4589" ~ paste0("3089", substr(ADM3_PCODE, 5, 7)),
                                                           TRUE ~ as.character(ADM3_PCODE))) 
cent <- SpatialPointsDataFrame(gCentroid(bang, byid=TRUE), 
                               bang@data, match.ID=FALSE)

cent_df <- as.data.frame(cbind(long = cent$x, 
                               lat = cent$y,
                               upazila_code=bang@data$upazila_code))

distanceBetweenUpazilas = expand.grid(ADM3_PCODE_home = bang@data$upazila_code, 
                                      ADM3_PCODE_dest = bang@data$upazila_code) %>%
  left_join(., cent_df, by = c("ADM3_PCODE_home" = "upazila_code"))  %>%
  left_join(., cent_df, by = c("ADM3_PCODE_dest" = "upazila_code")) %>%
  mutate(long.x = as.numeric(as.character(long.x)), 
         lat.x = as.numeric(as.character(lat.x)),
         long.y = as.numeric(as.character(long.y)),
         lat.y = as.numeric(as.character(lat.y)))

distance_upa = distVincentyEllipsoid(distanceBetweenUpazilas[,3:4],
                                     distanceBetweenUpazilas[,5:6])
distanceBetweenUpazilas = distanceBetweenUpazilas %>% mutate(distance_mi = distance_upa*0.000621371)

#saveRDS(distanceBetweenUpazilas, here("data", "misc", "distanceBetweenUpazilas.rds"))

## Functions needed for matrices ----
get_district_df = function(df, impute = TRUE){
  df2 = df %>% 
    filter( UpaFrom != "557384" & UpaTo != "557384") %>% # remove one missing upazila from crossmap
    left_join(., district_upa %>% select(upazila_code, District), by = c("UpaFrom" = "upazila_code")) %>%
    rename(DistrictFrom = District) %>%
    left_join(., district_upa %>% select(upazila_code, District), by = c("UpaTo" = "upazila_code")) %>%
    rename(DistrictTo = District) %>%
    group_by(DistrictFrom, DistrictTo, weekend) %>%
    summarise(value = sum(value, na.rm=T))
  
  return(df2)
}

# used in 01_process_cdr.R
get_impute_dhaka_m = function(operator_df, weekday){

  # read in files created from `02_impute_cdr_2017.R` 
  # 538 row x 21 col
  counts.out.weekday.2017 = readRDS(here("data", "cleaned_2017_cdr", "weekday.2017.m.rds"))
  counts.out.weekend.2017 = readRDS(here("data", "cleaned_2017_cdr", "weekend.2017.m.rds"))
  # 21 row x 538 col
  trips.in.weekday.2017 = readRDS(here("data", "cleaned_2017_cdr", "trips.in.weekday.2017.withmult.rds"))
  trips.in.weekend.2017 = readRDS(here("data", "cleaned_2017_cdr", "trips.in.weekend.2017.withmult.rds"))

  if(weekday){
    trips.in.2017 = trips.in.weekday.2017
    counts.out.2017 = counts.out.weekday.2017

  } else {
    trips.in.2017 = trips.in.weekend.2017
    counts.out.2017 = counts.out.weekend.2017
  }
  
  m = get_m(operator_df, weekday)
    
  # merge in 2017 trip counts (adjusted by multiplier) for upa to before calculating proportions
  m_sub = m[!(rownames(m) %in% rownames(trips.in.2017)),] # 523 x 544 - 21 missing dhaka upa
  m_sub_mergein =rbind(m_sub, trips.in.2017[, match(colnames(m_sub), colnames(trips.in.2017)),drop=FALSE])
  m_sub_mergein = m_sub_mergein[order(rownames(m_sub_mergein)), order(colnames(m_sub_mergein))]

  # merge in 2017 counts for upa from - missing column values
  count_m_sub = m_sub_mergein[,!(colnames(m_sub_mergein) %in% colnames(counts.out.2017))] # 544 x 523
  count_m_merge =cbind(count_m_sub, counts.out.2017[match(rownames(count_m_sub), rownames(counts.out.2017)), , drop=FALSE])
  count_m_merge = count_m_merge[order(rownames(count_m_merge)), order(colnames(count_m_merge))]

  return(count_m_merge)
}

# used in 01_process_cdr.R
get_weekend_m_2017 = function(operator_df, weekday){
  m = get_m(operator_df, weekday)
  prop_m = get_prop_m(m)
  
  prop_m = prop_m[order(rownames(prop_m)), order(colnames(prop_m))]
  prop_m[is.nan(prop_m)] <- 0 # if 0 data, treat as missing data
  
  pop_upa_2017 = read.csv(file = here("data","popcount_by_union.csv"), header=T) %>%
    rename("upazila_code" = "Upa_Code") %>%
    mutate(upazila_code = case_when(substr(upazila_code, 1, 4) == "3039" ~ paste0("4539", substr(upazila_code, 5, 7)),
                                    substr(upazila_code, 1, 4) == "3061" ~ paste0("4561", substr(upazila_code, 5, 7)),
                                    substr(upazila_code, 1, 4) == "3072" ~ paste0("4572", substr(upazila_code, 5, 7)),
                                    substr(upazila_code, 1, 4) == "3089" ~ paste0("4589", substr(upazila_code, 5, 7)),
                                    TRUE ~ as.character(upazila_code))) %>%
    group_by(upazila_code) %>%
    summarise(pop = sum(pop.tot)) %>%
    ungroup()
  
  pop_order_2017 = pop_upa_2017 %>% distinct(upazila_code, .keep_all = TRUE) %>% 
    filter(upazila_code %in% colnames(prop_m)) %>%
    arrange(match(upazila_code, colnames(prop_m)))
  
  m_norm = prop_m%*%diag(pop_order_2017$pop)
  colnames(m_norm) = rownames(m_norm)
  print(colSums(m_norm)) # check pop totals
  
  diag(m_norm) <- 0
  
  return(m_norm)
}

# used in 01_process_cdr.R
get_beta_distr_param = function(matrix){
  
  matrix[is.na(matrix)] <- 0
  num = c(matrix)
  denom = rep(unname(colSums(matrix, na.rm=T)), each = nrow(matrix))

  # doesn't matter conf.level or "type" because only using shape parameters
  bayes.results = binom.bayes(x = num, n = denom, type = "central", 
                              conf.level = 0.8) 
  
  alpha.mat = matrix(c(bayes.results$shape1), nrow = nrow(matrix))
  beta.mat = matrix(c(bayes.results$shape2), nrow = nrow(matrix))
  
  colnames(alpha.mat) = colnames(matrix)
  rownames(alpha.mat) = rownames(matrix)
  colnames(beta.mat) = colnames(matrix)
  rownames(beta.mat) = rownames(matrix)

  shape.mat = list(alpha = alpha.mat, beta = beta.mat)
  
  return(shape.mat)
}

## Differentiate zeroes and true NAs ====
# used in 01_process_cdr.R
differentiate_zero_na_df = function(df, 
                                    true_na_upa){ # list of true NA upazilas based on operator

  diff = df %>% 
    filter(UpaFrom != "557384" & UpaTo != "557384") %>%
    rename(raw_value = value) %>% 
    mutate(value =  case_when(UpaFrom != UpaTo & (UpaFrom %in% true_na_upa | UpaTo %in% true_na_upa) ~ NA_real_, # true NA
                              UpaFrom == UpaTo & (UpaFrom %in% true_na_upa | UpaTo %in% true_na_upa) ~ NA_real_, # true NA
                              UpaFrom != UpaTo & !(UpaFrom %in% true_na_upa) & !(UpaTo %in% true_na_upa) & is.na(raw_value) ~ 0, # zeroes
                              UpaFrom == UpaTo & !(UpaFrom %in% true_na_upa) & !(UpaTo %in% true_na_upa) & is.na(raw_value) ~ 0, # zeroes - missing dhaka ones we impute laeter                             # there are only NaNs, no 0
                              UpaFrom != UpaTo & !(UpaFrom %in% true_na_upa) & !(UpaTo %in% true_na_upa) & raw_value == 0 ~ 0, # zeroes # only op1 has mixed 0s and NAs for this
                              UpaFrom == UpaTo & !(UpaFrom %in% true_na_upa) & !(UpaTo %in% true_na_upa) & raw_value == 0 ~ 0, # there should be none of these
                              TRUE ~ raw_value)) %>%
    select(weekend, UpaFrom, UpaTo, value)

  return(diff)
}
