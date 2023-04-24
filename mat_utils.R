# Basic matrix utilities
# dependencies: tidyverse

get_df_from_m = function(m){
  df = as.data.frame(m) %>%
    rownames_to_column(var = "UpaTo") %>%
    pivot_longer(!UpaTo, names_to = "UpaFrom", values_to = "value")
  
  return(df)
}

half_mat = function(m){
  m[upper.tri(m)] <- -10000
  
  return(m)
}

# convert df to matrix
get_m = function(operator_df, weekday){
  
  if(any(grepl("UpaFrom", colnames(operator_df)))) {
    weekday_m = operator_df %>% 
      filter(weekend == FALSE & UpaFrom != "557384" & UpaTo != "557384") %>% # one unmapped upazila
      select(UpaFrom, UpaTo, value) %>%
      pivot_wider(names_from = UpaFrom, values_from = value) %>%
      remove_rownames %>% column_to_rownames(var="UpaTo") %>%
      data.matrix()
    
    weekend_m = operator_df %>% 
      filter(weekend == TRUE & UpaFrom != "557384" & UpaTo != "557384") %>% 
      select(UpaFrom, UpaTo, value) %>%
      pivot_wider(names_from = UpaFrom, values_from = value) %>%
      remove_rownames %>% column_to_rownames(var="UpaTo") %>%
      data.matrix()
  } else {
    weekday_m = operator_df %>% 
      filter(weekend == FALSE) %>% 
      select(DistrictFrom, DistrictTo, value) %>%
      pivot_wider(names_from = DistrictFrom, values_from = value) %>%
      remove_rownames %>% column_to_rownames(var="DistrictTo") %>%
      data.matrix()
    
    weekend_m = operator_df %>% 
      filter(weekend == TRUE) %>% 
      select(DistrictFrom, DistrictTo, value) %>%
      pivot_wider(names_from = DistrictFrom, values_from = value) %>%
      remove_rownames %>% column_to_rownames(var="DistrictTo") %>%
      data.matrix()
  }
  
  if(weekday){
    return(weekday_m)
  } else {
    return(weekend_m)
  }
  
}

get_prop_m = function(matrix){
  matrix[is.na(matrix)] <- 0 # need this for prop.table
  
  prop_m = prop.table(matrix, 2)
  #colSums(prop_m) # check 1s
  
  return(prop_m)
}

# normalize by pop 2020 - load pop 2020 as pop_input
get_norm_m = function(matrix, pop_input){
  
  # do this to get this function to run; then differentiate between 0 and NA later
  matrix[is.na(matrix)] <- 0
  
  # affix population data
  if(ncol(matrix) > 500) { # upazila
    
    pop_order = pop_input %>% distinct(upazila_code, .keep_all = TRUE) %>% arrange(match(upazila_code, colnames(matrix)))
    
  } else {
    pop_order = pop_input %>%
      left_join(., district_upa %>% select(District, upazila_code), by = "upazila_code") %>%
      group_by(District) %>%
      summarise(pop = sum(pop)) %>%
      arrange(match(District, colnames(matrix)))
    
  }
  m_norm = matrix%*%diag(pop_order$pop)
  colnames(m_norm) = rownames(m_norm)
  print(colSums(m_norm, na.rm=T)) # check pop totals
  
  return(m_norm)
}

get_na_upa_list = function(df){
  missing_dhakaupa = c("302602", "302605", "302606", "302609", "302610", "302611", "302624",
                       "302629", "302632", "302630", "302633", "302637", "302663", "302665",
                       "302667", "302674", "302675", "302680", "302692", "302693", "302696")
  na_list = df %>% 
    filter(weekend == FALSE) %>% # same if use weekend == TRUE
    filter(UpaFrom != "557384" & UpaTo != "557384") %>% # one unmapped upazila
    filter(!(UpaFrom %in% missing_dhakaupa)) %>%
    group_by(UpaFrom) %>%
    tally(is.na(value)) %>%
    filter(n== max(n)) %>% 
    pull(UpaFrom)
  
  return(na_list)
}

get_mat_with_na = function(matrix, pop_input, mobility_source){
  
  # hard-code lists from get_na_upa_list function
  if(mobility_source == "Operator 1"){
    true_na_upa = c("104240", "104243", "104273", "104284", "200395", "201207", "201967", "208421", 
                    "208429", "208447", "302990", "308247", "404794", "501085", "609027")
  } else if (mobility_source == "Operator 2"){
    true_na_upa = c("202245", "208478", "302695", "304802", "304833", "304859", "307218",
                    "307247", "307263", "404147", "553221", "603602", "609086", "609092")
  } else if (mobility_source == "Operator 3"){
    true_na_upa = c( "200304", "200373", "200389", "200391", "200395", "201207", "201933", "204661", "204665", "204677",
                     "208407", "208421", "208429", "208447", "208458", "208475", "302640", "302648", "302990", "307218",
                     "308247", "400138", "404433", "404721", "405039", "405595", "501094", "508185", "552769", "555239",
                     "557364", "558576", "609027", "609127", "609131")
  } else if (mobility_source == "Meta" | mobility_source == "Gravity Model"){
    true_na_upa = NULL
  }
  
  matrix[is.na(matrix)] <- 0 
  matrix[, colnames(matrix) %in% true_na_upa] <- NaN
  matrix[rownames(matrix) %in% true_na_upa,] <- NaN
  
  return(matrix)
  
}

get_norm_sym_m = function(matrix, pop_input, mobility_source, sym = TRUE){

  # proportions
  prop_m = get_prop_m(matrix)
  
  # normalize matrix to 2020 pop
  norm_m = get_norm_m(matrix = prop_m, pop_input = pop_input)

  # make matrix symmetric
  if(sym){
    sym_norm_m = (norm_m + t(norm_m))/2
  } else {
    sym_norm_m = norm_m
  }
  
  if(mobility_source == "Operator 1"){
    true_na_upa = c("104240", "104243", "104273", "104284", "200395", "201207", "201967", "208421", 
                    "208429", "208447", "302990", "308247", "404794", "501085", "609027")
  } else if (mobility_source == "Operator 2"){
    true_na_upa = c("202245", "208478", "302695", "304802", "304833", "304859", "307218",
                    "307247", "307263", "404147", "553221", "603602", "609086", "609092")
  } else if (mobility_source == "Operator 3"){
    true_na_upa = c( "200304", "200373", "200389", "200391", "200395", "201207", "201933", "204661", "204665", "204677",
                     "208407", "208421", "208429", "208447", "208458", "208475", "302640", "302648", "302990", "307218",
                     "308247", "400138", "404433", "404721", "405039", "405595", "501094", "508185", "552769", "555239",
                     "557364", "558576", "609027", "609127", "609131")
  } else if (mobility_source == "Meta"){
    true_na_upa = NULL
  } 
  
  sym_norm_m[is.na(sym_norm_m)] <- 0 
  sym_norm_m[, colnames(sym_norm_m) %in% true_na_upa] <- NaN
  sym_norm_m[rownames(sym_norm_m) %in% true_na_upa,] <- NaN
  
  print(colSums(sym_norm_m, na.rm=T)) # check pop totals
  
  return(sym_norm_m)
  
}

get_district_m = function(matrix){
  
  aggregated_df = get_df_from_m(matrix) %>%
    left_join(., district_upa %>% select(upazila_code, District), by = c("UpaFrom" = "upazila_code")) %>%
    rename(DistrictFrom = District) %>%
    left_join(., district_upa %>% select(upazila_code, District), by = c("UpaTo" = "upazila_code")) %>%
    rename(DistrictTo = District) %>%
    group_by(DistrictFrom, DistrictTo) %>%
    summarise(value = sum(value, na.rm=T)) %>%
    ungroup()
  
  matrix = aggregated_df %>% 
    select(DistrictFrom, DistrictTo, value) %>%
    pivot_wider(names_from = DistrictFrom, values_from = value) %>%
    remove_rownames %>% column_to_rownames(var="DistrictTo") %>%
    data.matrix()
  
  return(matrix)
  
}
