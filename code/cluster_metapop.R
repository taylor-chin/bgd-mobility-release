
require(tidyverse)

setwd('/n/home11/tchin/bgd_mobility')
source("./mat_utils.R")
source("./05_spatial_seir.R")

# load pop data
load("full_pop_age_gender.Rdata")
upa <- read.csv("popcount_by_union.csv", header=T) %>%
  mutate(upazila_code = as.character(Upa_Code)) %>%
  dplyr::select(upazila_code, Upazila) %>% distinct()

upa_pop <- full_pop %>% 
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

district_upa <- read.csv("popcount_by_union.csv", header=T) %>%
  mutate(upazila_code = as.character(Upa_Code), Dis_Code = as.character(Dis_Code), District = as.character(District)) %>%
  dplyr::select(upazila_code, Upazila, Dis_Code, District) %>% distinct()

district_upa$District[district_upa$District == "Cox'S Bazar"] <- "Cox's Bazar"

# updated this to remove the double entry for upa code "303334" in the unique Upazila name
district_upa <- district_upa %>% filter(Upazila != "Kaliganj Paurashava") # two entries for this upa_code "303334"


# parameters
seed_cities = c("Dhaka", "Chittagong", "Panchagarh")
r0 = c(1.2, 1.3, 1.5, 2)
matrix = c("op1_weekday_beta_param,op1_weekend_beta_param",
           "op2_weekday_beta_param,op2_weekend_beta_param",
           "op3_weekday_beta_param,op3_weekend_beta_param",
           "fb_weekday_beta_param_rawm,fb_weekend_beta_param_rawm", 
          "grav_mat.w.lit.params,grav_mat.w.lit.params") # use upa level gravity model if upazila level sims:
           #grav_mat.upa.w.lit.params,grav_mat.upa.w.lit.params")


params = expand.grid(matrix = matrix, seed_city = seed_cities, r0 = r0) %>%
  separate(matrix,c("m1","m2"),sep=",") %>%
  mutate(mobility_source = case_when(grepl( "op1", m1) ~ "Operator 1", 
                                     grepl( "op2", m1) ~ "Operator 2", 
                                     grepl( "op3", m1) ~ "Operator 3",
                                     grepl( "fb", m1)  ~ "Meta",
                                     TRUE ~ "Gravity model")) %>%
  mutate(filename = paste0(sub("_.*", "", m1), ".",seed_city,".r0",r0)) %>%
  arrange(match(mobility_source, c("Gravity model", "Operator 1", "Operator 2", "Operator 3", "Meta")))  #

j <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

set.seed(100)

m1_file = readRDS(paste0(params$m1[j],".rds"))
m2_file = readRDS(paste0(params$m2[j],".rds"))

if(j <= 12){
  # no uncertainty for gravity model
  # change this function to _upa level versions if doing upazila level sims
  results = run_metapop_nomatuncert(m1 = m1_file,
                        m2 = m2_file,
                        mobility_source = params$mobility_source[j],
                        seed_city = params$seed_city[j],
                        pop_data = upa_pop,
                        r0 = params$r0[j])
  saveRDS(results[[1]], file = paste0("./test/",params$filename[j],".district.beta.all.vF",".rds"))
  saveRDS(results[[2]], file = paste0("./test/",params$filename[j],".district.beta.avg.vF",".rds"))
  saveRDS(results[[3]], file = paste0("./test/",params$filename[j],".district.beta.sumbgd.vF",".rds"))
  saveRDS(results[[4]], file = paste0("./test/",params$filename[j],".district.beta.finalsize.vF",".rds"))
  saveRDS(results[[5]], file = paste0("./test/",params$filename[j],".district.beta.numdist30.vF",".rds"))

} else if(j > 12) {
  # change this function to _upa level versions if doing upazila level sims
  results = run_metapop(m1.param = m1_file,
                        m2.param = m2_file,
                        mobility_source = params$mobility_source[j],
                        seed_city = params$seed_city[j],
                        pop_data = upa_pop,
                        r0 = params$r0[j])

  saveRDS(results[[1]], file = paste0("./test/",params$filename[j],".district.beta.all.vF",".rds"))
  saveRDS(results[[2]], file = paste0("./test/",params$filename[j],".district.beta.avg.vF",".rds"))
  saveRDS(results[[3]], file = paste0("./test/",params$filename[j],".district.beta.sumbgd.vF",".rds"))
  saveRDS(results[[4]], file = paste0("./test/",params$filename[j],".district.beta.finalsize.vF",".rds"))
  saveRDS(results[[5]], file = paste0("./test/",params$filename[j],".district.beta.numdist30.vF",".rds"))
  
}


