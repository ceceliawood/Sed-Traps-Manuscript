#Stratification Strength/Duration Code

rm(list=ls(all=TRUE))
library(dplyr)
library(tidyr)
library(readxl)
library(lubridate)
library(readr)
library(ggplot2)
library(reshape2)
library(scales)
library(stringr)

#Using Freya's functions


#' get daily targets
#' 
#' @param infiles vector of EDI data files
#' @param interpolate should interpolation be carried out
#' @return a targets dataframe in VERA format

get_targets <- function(infiles, interpolate = T, maxgap = 12) {
  
  
  standard_names <- data.frame(variable_new = c('Temp_C','SpCond_uScm', 'Chla_ugL', 'fDOM_QSU'),
                               variable = c('EXOTemp_C_1', 'EXOSpCond_uScm_1', 'EXOChla_ugL_1', 'EXOfDOM_QSU_1'))
  targets <- NULL
  
  # Load data
  message('Reading data from EDI...')
  for (i in 1:length(infiles)) {
    df <- read_csv(infiles[i], show_col_types = F, progress = F) |> 
      filter(Site == 50) |> 
      # mutate(site_id = ifelse(Reservoir == 'BVR', 'bvre', ifelse(Reservoir == 'FCR', 'fcre', Reservoir))) |> 
      rename(datetime = DateTime,
             site_id = Reservoir) |> 
      filter(site_id %in% c('FCR', 'BVR')) |> 
      select(-Site)
    
    df_flags <- df |> 
      select(any_of(c('datetime', 'site_id')) | contains('Flag') & contains(standard_names$variable)) |> 
      pivot_longer(cols = contains('Flag'),
                   names_to = 'variable', 
                   values_to = 'flag_value', 
                   names_prefix = 'Flag_')
    
    df_observations <- df |> 
      select(any_of(c('datetime', 'site_id')) | contains(standard_names$variable) & !contains('Flag')) |> 
      pivot_longer(cols = -any_of(c('site_id', 'datetime')),
                   names_to = 'variable', 
                   values_to = 'observation', 
                   names_prefix = 'Flag_')
    
    
    # Filter any flagged data
    message('Filtering flags...')
    df_long <- inner_join(df_observations, df_flags, 
                          # by = c('site_id', 'datetime', 'depth_m', 'variable'),
                          relationship = 'many-to-many') |> 
      na.omit() |> 
      filter(!flag_value %in% c(9, 7, 2, 1, 5, 3)) |> 
      select(-contains('flag')) |> 
      mutate(depth_m = as.numeric(str_split_i(variable, "_", 3)), 
             depth_m = ifelse(str_detect(variable, 'EXO') & site_id == 'FCR', 1.6, 
                              ifelse(str_detect(variable, 'EXO') & site_id == 'BVR', 1.5, depth_m)),
             variable = str_split_i(variable, "\\.", 1)) |> 
      full_join(standard_names) |> 
      mutate(variable = variable_new) |> 
      select(-variable_new)
    
    # get DO
    df_DO <- df |>
      select(datetime, site_id, any_of(c('RDO_mgL_9_adjusted', 'RDO_mgL_13', 'EXODO_mgL_1', 'EXODO_mgL_1.5'))) |>
      pivot_longer(cols = contains('DO'), names_to = 'variable', values_to = 'observation') |>
      mutate(depth_m = as.numeric(str_split_i(variable, "_", 3)), 
             depth_m = ifelse(str_detect(variable, 'EXO') & site_id == 'FCR', 1.6, 
                              ifelse(str_detect(variable, 'EXO') & site_id == 'BVR', 1.5, depth_m)),
             variable = 'DO_mgL') 
    
    #get bottom temp
    df_bottomT <- df |>
      select(datetime, site_id, starts_with(c('ThermistorTemp', 'Flag_ThermistorTemp'))) |>
      pivot_longer(cols = contains('Thermistor'), names_to = 'variable', values_to = 'observation') |>
      mutate(depth_m = str_split_i(variable, "_", -1),
             depth_m = as.numeric(ifelse(depth_m == 'surface', 0, depth_m))) |> 
      filter(depth_m == max(depth_m)) |> 
      pivot_wider(names_from = variable, values_from = observation) |> 
      filter(if_any(starts_with("Flag"),  ~.x == 0)) |> 
      select(-starts_with('Flag')) |> 
      rename(observation = starts_with('Thermistor')) |> 
      mutate(variable = 'Temp_C')
    
    # combine
    df_all <- df_long |>
      dplyr::bind_rows(df_DO) |> 
      dplyr::bind_rows(df_bottomT) |>
      na.omit() |> 
      dplyr::select(datetime, site_id, depth_m, observation, variable) |>
      dplyr::mutate(observation = ifelse(!is.finite(observation),NA,observation)) |> 
      tsibble::as_tsibble(key = any_of(c("site_id", "depth_m", "variable")),
                          index = "datetime") |>
      tsibble::fill_gaps() |> 
      as_tibble()
    
    
    targets <- bind_rows(targets, df_all) 
  }
  
  if (interpolate == T) {
    test <- targets |> 
      group_by(site_id, depth_m, variable) |> 
      arrange(datetime) |> 
      mutate(observation = imputeTS::na_interpolation(observation, 'linear', maxgap = maxgap, rule = 1))
  }
  
  return(targets)
}

# ================================================================#

# =================== Temperature profiles ==================

get_temp_profiles <- function(current_file = 'none', historic_file){
  source('R/find_depths.R')
  
  if (current_file != 'none') {
    message('reading ', current_file)
    current_df <- readr::read_csv(current_file, show_col_types = F) |>
      dplyr::filter(Site == 50) |>
      dplyr::select(Reservoir, DateTime,
                    dplyr::starts_with('ThermistorTemp'))
    
    if (current_df$Reservoir[1] == 'BVR') {
      bvr_depths <- find_depths(data_file = current_file,
                                depth_offset = "https://raw.githubusercontent.com/FLARE-forecast/BVRE-data/bvre-platform-data-qaqc/BVR_Depth_offsets.csv",
                                output <- NULL,
                                date_offset <- "2021-04-05",
                                offset_column1<- "Offset_before_05APR21",
                                offset_column2 <- "Offset_after_05APR21") |>
        dplyr::filter(variable == 'ThermistorTemp') |>
        dplyr::select(Reservoir, DateTime, variable, depth_bin, Position)
      
      current_df_1 <- current_df  |>
        tidyr::pivot_longer(cols = starts_with('ThermistorTemp'),
                            names_to = c('variable','Position'),
                            names_sep = '_C_',
                            values_to = 'observation') |>
        dplyr::mutate(date = lubridate::as_date(DateTime),
                      Position = as.numeric(Position)) |>
        na.omit() |>
        dplyr::left_join(bvr_depths,
                         by = c('Position', 'DateTime', 'Reservoir', 'variable')) |>
        dplyr::group_by(date, Reservoir, depth_bin) |>
        dplyr::summarise(observation = mean(observation, na.rm = T),
                         n = dplyr::n(),
                         .groups = 'drop') |>
        dplyr::mutate(observation = ifelse(n < 144/3, NA, observation), # 144 = 24(hrs) * 6(10 minute intervals/hr)
                      Reservoir = 'bvre') |>
        
        dplyr::rename(site_id = Reservoir,
                      datetime = date,
                      depth = depth_bin) |>
        dplyr::select(-n) |> 
        dplyr::mutate(depth = as.character(depth))
    }
    
    # read in differently for FCR
    if (current_df$Reservoir[1] == 'FCR') {
      current_df_1 <- current_df |>
        tidyr::pivot_longer(cols = starts_with('ThermistorTemp'),
                            names_to = 'depth',
                            names_prefix = 'ThermistorTemp_C_',
                            values_to = 'observation') |>
        dplyr::mutate(#Reservoir = ifelse(Reservoir == 'FCR',
          #                 'fcre',
          #                ifelse(Reservoir == 'BVR',
          #                      'bvre', NA)),
          date = lubridate::as_date(DateTime)) |>
        na.omit() |>
        dplyr::group_by(date, Reservoir, depth) |>
        dplyr::summarise(observation = mean(observation, na.rm = T),
                         n = dplyr::n(),
                         .groups = 'drop') |>
        dplyr::mutate(observation = ifelse(n < 144/2, NA, observation),
                      depth = as.character(depth)) |> # 144 = 24(hrs) * 6(10 minute intervals/hr)
        dplyr::rename(site_id = Reservoir,
                      datetime = date) |>
        dplyr::select(-n)
    }
    message('Current file ready')
  } else {
    current_df_1 <- NULL
    message('No current file')
  }
  
  # read in historical data file
  # EDI
  # infile <- tempfile()
  # try(download.file(historic_file, infile, method="curl"))
  # if (is.na(file.size(infile))) download.file(historic_file,infile,method="auto")
  
  historic_df <- readr::read_csv(historic_file, show_col_types = FALSE) |>
    dplyr::filter(Site == 50) |>
    dplyr::select(Reservoir, DateTime,
                  dplyr::starts_with('ThermistorTemp'))
  
  # Extract depths for BVR
  if (historic_df$Reservoir[1] == 'BVR') {
    bvr_depths <- find_depths(data_file = historic_file,
                              depth_offset = "https://raw.githubusercontent.com/FLARE-forecast/BVRE-data/bvre-platform-data-qaqc/BVR_Depth_offsets.csv",
                              output <- NULL,
                              date_offset <- "2021-04-05",
                              offset_column1<- "Offset_before_05APR21",
                              offset_column2 <- "Offset_after_05APR21") |>
      dplyr::filter(variable == 'ThermistorTemp') |>
      dplyr::select(Reservoir, DateTime, variable, depth_bin, Position)
    
    historic_df_1 <- historic_df |>
      tidyr::pivot_longer(cols = starts_with('ThermistorTemp'),
                          names_to = c('variable','Position'),
                          names_sep = '_C_',
                          values_to = 'observation') |>
      dplyr::mutate(date = lubridate::as_date(DateTime),
                    Position = as.numeric(Position)) |>
      na.omit() |>
      dplyr::left_join(bvr_depths,
                       by = c('Position', 'DateTime', 'Reservoir', 'variable')) |>
      dplyr::group_by(date, Reservoir, depth_bin) |>
      dplyr::summarise(observation = mean(observation, na.rm = T),
                       n = dplyr::n(),
                       .groups = 'drop') |>
      dplyr::mutate(observation = ifelse(n < 144/3, NA, observation), # 144 = 24(hrs) * 6(10 minute intervals/hr)
                    Reservoir = 'bvre') |>
      dplyr::rename(site_id = Reservoir,
                    datetime = date,
                    depth = depth_bin) |>
      dplyr::select(-n) |> 
      dplyr::mutate(depth = as.character(depth))
  }
  
  if (historic_df$Reservoir[1] == 'FCR') {
    historic_df_1 <- historic_df |>
      tidyr::pivot_longer(cols = starts_with('ThermistorTemp'),
                          names_to = 'depth',
                          names_prefix = 'ThermistorTemp_C_',
                          values_to = 'observation') |>
      dplyr::mutate(#Reservoir = ifelse(Reservoir == 'FCR',
        #                  'fcre',
        #                 ifelse(Reservoir == 'BVR',
        #                       'bvre', NA)),
        date = lubridate::as_date(DateTime)) |>
      dplyr::group_by(date, Reservoir, depth)  |>
      dplyr::summarise(observation = mean(observation, na.rm = T),
                       n = dplyr::n(),
                       .groups = 'drop') |>
      dplyr::mutate(observation = ifelse(n < 6/2, NA, observation)) |> # 6 = 6(10 minute intervals/hr)
      dplyr::rename(site_id = Reservoir,
                    datetime = date)|>
      dplyr::select(-n) |> 
      dplyr::mutate(depth = as.character(depth))
  }
  
  message('EDI file ready')
  
  ## manipulate the data files to match each other
  
  
  ## bind the two files using row.bind()
  final_df <- dplyr::bind_rows(historic_df_1, current_df_1) |>
    dplyr::mutate(variable = 'Temp_C_mean',
                  depth = as.numeric(ifelse(depth == "surface", 0.1, depth))) |>
    rename(depth_m = depth)
  
  final_df <- final_df |>
    mutate(observation = ifelse(is.nan(observation), NA, observation)) |>
    drop_na(depth_m)
  ## Match data to flare targets file
  # Use pivot_longer to create a long-format table
  # for time specific - use midnight UTC values for daily
  # for hourly
  
  ## return dataframe formatted to match FLARE targets
  return(final_df)
}


calc_strat_dates <- function(density_diff = 0.1,
                             temp_profiles) {
  
  ## extract the depths that will be used to calculate the density difference (surface, bottom)
  depths_use <- temp_profiles |>
    na.omit() |> 
    dplyr::group_by(datetime, site_id) |>
    dplyr::summarise(top = min(as.numeric(depth_m, na.rm = T)),
                     bottom = max(as.numeric(depth_m, na.rm = T)),.groups = 'drop') |>
    tidyr::pivot_longer(cols = top:bottom, 
                        names_to = 'location',
                        values_to = 'depth_m')
  
  sites <- distinct(depths_use, site_id) |> pull()
  
  strat_dates <- NULL
  
  for (site in sites) {
    temp_profile_site <- filter(temp_profiles, site_id == site)
    # need a full timeseries
    all_dates <- data.frame(datetime = seq.Date(min(temp_profile_site$datetime), 
                                                max(temp_profile_site$datetime),
                                                'day'))
    density_obs <-
      filter(depths_use, site_id == site) |> 
      inner_join(na.omit(temp_profile_site), by = join_by(datetime, site_id, depth_m)) |> 
      mutate(density = rLakeAnalyzer::water.density(observation)) |> 
      select(datetime, site_id, density, observation, location) |> 
      pivot_wider(values_from = c(density, observation), names_from = location, id_cols = c(datetime, site_id)) |> 
      full_join(all_dates, by = 'datetime') |> 
      mutate(dens_diff = density_bottom - density_top,
             strat = ifelse(abs(dens_diff > 0.1) & observation_top > observation_bottom, 1, 0),
             strat = imputeTS::na_interpolation(strat, option = 'linear'))
    
    
    # extract the dates of the stratified periods
    #using a loop function to go through each year and do the rle function
    
    strat <- data.frame(year = unique(year(density_obs$datetime)), 
                        length = NA,
                        start = NA,
                        end = NA)
    
    for (i in 1:nrow(strat)) {
      year_use <- strat$year[i]
      
      temp.dens <- density_obs %>%
        filter(year(datetime) == year_use)
      
      if (nrow(temp.dens) >= 300) {
        #run length encoding according to the strat var
        temp.rle <- rle(temp.dens$strat)
        
        #what is the max length for which the value is "norm"
        strat$length[i] <- max(temp.rle$lengths[temp.rle$values==1], 
                               na.rm = T)
        
        #stratification dates
        rle.strat <- data.frame(strat = temp.rle$values, 
                                lengths = temp.rle$lengths)
        
        # Get the end of ech run
        rle.strat$end <- cumsum(rle.strat$lengths)
        # Get the start of each run
        rle.strat$start <- rle.strat$end - rle.strat$lengths + 1
        
        # Sort rows by whehter it is stratified or not
        rle.strat <- rle.strat[order(rle.strat$strat), ]
        
        start.row <- rle.strat$start[which(rle.strat$length == max(rle.strat$lengths)
                                           & rle.strat$strat == 1)] 
        #gets the row with the start date
        #of the run which has the max length and is 1
        
        end.row <- rle.strat$end[which(rle.strat$length == max(rle.strat$lengths)
                                       & rle.strat$strat == 1)] 
        #gets the row with the end date
        #of the run which has the max length and is TRuE
        
        strat$start[which(strat$year == year_use)] <- as.character(temp.dens$datetime[start.row])
        strat$end[which(strat$year == year_use)] <- as.character(temp.dens$datetime[end.row])
        
        strat$site_id <- site
      }
      
    } 
    strat_dates <- bind_rows(strat, strat_dates)
    message(site)
  }
  
  return(na.omit(strat_dates))
}


get_targets_sample  <- function(infiles, start_date, end_date) {
  
  # list of standardised column names
  standard_names <- c(site_id = "Reservoir", 
                      depth_m = "Depth_m",
                      datetime = "DateTime")
  
  final_df <- NULL
  
  for (i in 1:length(infiles)) {
    
    df <- read_csv(infiles[i], show_col_types = F, progress = F) |> 
      filter(Site == 50) |> 
      rename(any_of(standard_names))  |> 
      # mutate(site_id = ifelse(site_id == 'BVR', 'bvre', ifelse(site_id == 'FCR', 'fcre', site_id))) |> 
      filter(site_id %in% c('FCR', 'BVR')) |> 
      select(-Site)
    
    df_flags <- df |> 
      select(any_of(c('datetime', 'site_id', 'depth_m')) | contains('Flag')) |> 
      pivot_longer(cols = contains('Flag'),
                   names_to = 'variable', 
                   values_to = 'flag_value', 
                   names_prefix = 'Flag_')
    
    df_observations <- df |> 
      select(-contains('Flag')) |> 
      pivot_longer(cols = -any_of(c('site_id', 'datetime', 'depth_m', 'Rep')),
                   names_to = 'variable', 
                   values_to = 'observation', 
                   names_prefix = 'Flag_')
    
    df_long <- inner_join(df_observations, df_flags, 
                          # by = c('site_id', 'datetime', 'depth_m', 'variable'),
                          relationship = 'many-to-many') |> 
      na.omit() |> 
      filter(!flag_value %in% c(9, 5, 2)) |> 
      select(-contains('flag')) |> 
      group_by(pick(any_of(c('site_id', 'datetime', 'depth_m', 'variable')))) |> 
      summarise(observation = mean(observation), .groups = 'drop')
    
    
    # Combine with other dataframe
    final_df <- bind_rows(final_df, df_long) |> mutate(method = 'sample')
    
  }
  
  return(final_df) 
}

max_na <- function(ts) {
  rle_na <- rle(is.na(ts))
  
  rle_na <- data.frame(value = rle_na$values,
                       length = rle_na$lengths) 
  
  if (sum(rle_na$value == T) > 0) {
    longest_na <- rle_na |> 
      filter(value == T) |> 
      summarise(max(length)) |> 
      pull()
  } else {
    longest_na <- 0
  }
  return(longest_na)
}


n_cont <- function(ts) {
  
  rle_na <- rle(is.na(ts))
  
  rle_na <- data.frame(value = rle_na$values,
                       length = rle_na$lengths) 
  
  n_cont <- rle_na |> 
    filter(value == F) |> 
    summarise(n()) |> pull()
  
  return(n_cont)
}

fcre_EDI <- "https://pasta.lternet.edu/package/data/eml/edi/271/8/fbb8c7a0230f4587f1c6e11417fe9dce"
bvre_EDI <- "https://pasta.lternet.edu/package/data/eml/edi/725/4/9adadd2a7c2319e54227ab31a161ea12"

fcre_temp_profile <- get_temp_profiles(current_file = 'none', historic_file = fcre_EDI)
bvre_temp_profile <- get_temp_profiles(current_file = 'none', historic_file = bvre_EDI) %>%
  mutate(site_id = ifelse(site_id == "bvre",
                          "BVR",
                          site_id))

fcre_strat_dates <- calc_strat_dates(density_diff = 0.1, temp_profiles = fcre_temp_profile) |> 
  mutate(start = as_date(start),
         end = as_date(end))
bvre_strat_dates <- calc_strat_dates(density_diff = 0.1, temp_profiles = bvre_temp_profile) |> 
  mutate(start = as_date(start),
         end = as_date(end))
strat_dates <- bind_rows(fcre_strat_dates, bvre_strat_dates)

library(rLakeAnalyzer)
jpeg("figures/stratification by site.jpg", res = 300, height = 4, width = 6, units = "in")
fcre_temp_profile %>%
  bind_rows(bvre_temp_profile) %>%
  filter(year(datetime) >= 2021) %>%
  mutate(observation = ifelse(site_id == "BVR" & datetime == "2022-06-21" & depth_m == 0,
                              NA,
                              observation)) %>%
  group_by(datetime, site_id) %>%
  summarize(bf = max(buoyancy.freq(observation, depth_m))) %>%
  ggplot(aes(x = datetime, y = bf, color = site_id))+
  geom_point()+
  geom_vline(aes(xintercept = start, color = site_id, lty = "Onset"), data = strat_dates)+
  geom_vline(aes(xintercept = end, color = site_id, lty = "Turnover"), data = strat_dates)+
  ylab("Maximum buoyancy frequency (1/s)")+
  scale_linetype_discrete(name = "Stratification")+
  facet_wrap(~site_id, ncol = 1)+
  theme_bw()+
  theme(axis.title.x = element_blank())
dev.off()

jpeg("figures/bf_one2one.jpg", res = 300, height = 4, width = 4, units = "in")
fcre_temp_profile %>%
  bind_rows(bvre_temp_profile) %>%
  filter(year(datetime) >= 2021) %>%
  mutate(observation = ifelse(site_id == "BVR" & datetime == "2022-06-21" & depth_m == 0,
                              NA,
                              observation)) %>%
  group_by(datetime, site_id) %>%
  summarize(bf = max(buoyancy.freq(observation, depth_m))) %>%
  pivot_wider(names_from = site_id, values_from = bf) %>%
  ggplot(aes(x = FCR, y = BVR))+
  geom_point()+
  geom_abline()+
  geom_smooth(method = "lm")+
  theme_bw()+
  xlab("FCR buoyancy frequency (1/s)")+
  ylab("BVR buoyancy frequency (1/s)")
dev.off()

strat_dates %>%
  select(-length) %>%
  pivot_longer(start:end) %>%
  mutate(value = yday(value)) %>%
  pivot_wider(names_from = site_id, values_from = value) %>%
  ggplot(aes(x = FCR, y = BVR, color = name))+
  geom_point()+
  geom_abline()+
  facet_wrap(~name, scales = "free")
