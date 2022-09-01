##########################################################################
# Updated version of the R code for the analysis in:
#
#  "Heatwave exposure and preterm birth in China: attributable disease burden and
#   human capital consequences"
# 
#
# Update: 30 Aug. 2022
# * an updated version of this code is available at:
#   https://github.com/gasparrini/2019_vicedo-cabrera_Epidem_Rcodedata
#
# Requirements:
#  - R version 4.1.0
#  - 'readr' package version 1.4.0
#  - 'tidyverse' package version 1.3.1
##########################################################################
rm(list=ls())
gc()
library(ncdf4)
library(readr)
library(tidyverse)

##########################################################################
# 03. CALIBRATE THE SIMULATED TEMPERATURES, IDENTIFY HEATWAVES FOR EACH MODEL
#     AND CALCULATE THE HEATWAVE DAYS FOR FACTUAL AND COUNTERFACTUAL SCENARIOS
#
# 031.Before calculation, we first reorganize the observed and simulated  
#     temperature series into the form we need for calibration
#        
#
# LOAD RELEVANT DATA
#  -ERA5_1979-2019_DailyX4_Tmax_daymax.nc, providing gridded daily maximum temperature
#     series during 1979-2019.
#  -Daily_Tmax_2020.nc, providing gridded daily maximum temperature series in 2020.
#  -shp file, providing the longitudes and latitudes at a 0.5
#     degree resolution for all provinces in China
#  -Simulated tempertures in "hist" and "nat" files 

##########################################################################
##########################################################################
## REORGANIZE THE OBSERVED TEMPERATURE SERIES.
# Load the observed temperatures
ncintmax <- nc_open("/Volumes/Grace_PSSD/analysis/data/2.ERA5/era5_1979-2019_dailyx4_tmax_daymax.nc")

lon <- ncvar_get(ncintmax,"lon")
nlon <- dim(lon)
lat <- ncvar_get(ncintmax,"lat")
nlat <- dim(lat)
time <- ncvar_get(ncintmax,"time")
tmax_array1 <- ncvar_get(ncintmax,"mx2t_NON_CDM")
dim(tmax_array1) 


ncintmax3 <- nc_open("/Volumes/Grace_PSSD/analysis/data/2.ERA5/daily_tmax_2020.nc")
tmax_array3 <- ncvar_get(ncintmax3,"Tmax")

tmax_array2 <- tmax_array3
for (i in 1:141) {
  for (k in 1:366) {
    for (j in 1:111) {
      tmax_array2[i,(112-j),k] <- tmax_array3[i,j,k]
    }
  }
}

tmax_array_obs <- abind(tmax_array1,tmax_array2,along = 3)

# Due to simulated temperature only have 365 days for each year,so we extract the 365*42 days from
#  observed data correspondingly.
text <- tmax_array_obs[,,1:15330] #365*42 YEARS = 15330 DAYS
dim(tmax_array_obs)

# Reorganize the observe data 
tmax_obs_matrix <- matrix(nrow = nlon*nlat, ncol = (dim(text)[3])) #
for (yr in c(1:(dim(text)[3]))) {
  tmax_obs_matrix[ , yr] <- c(text[ , , yr])
}

colnames(tmax_obs_matrix) <- 1:15330
tmax_obs_tibble <- as.data.frame.matrix(tmax_obs_matrix)

country_shp <- sf::st_read("/Users/grace/data/Analysed/zyl/shp/全国数据.shp")
province_shp <- sf::st_read("/Users/grace/data/Analysed/zyl/shp/省级数据-带标签.shp")


all_lon_lat <- expand.grid("lon" = lon, "lat" = lat)
lon_lat_with_locaion <- data.frame(matrix(ncol = length(province_shp$ENAME)+1))
colnames(lon_lat_with_locaion) <- c(province_shp$SNAME, "china_boundary")
lon_lat_with_locaion <- all_lon_lat %>% bind_cols(lon_lat_with_locaion)
lon_lat_to_st_point <- st_sfc(lapply(1:nrow(all_lon_lat), 
                                     FUN = function(i){return(st_point(c(all_lon_lat$lon[i], all_lon_lat$lat[i])))}), 
                              crs = 4326)


for (temp_province in c(province_shp$SNAME, "china_boundary")) {
  if (temp_province == "china_boundary") {
    temp_geometry <- country_shp$geometry[1]
  } else {
    temp_geometry <- province_shp %>% filter(SNAME == temp_province) %>% pull(geometry) 
  }
  
  temp_result <- st_intersects(lon_lat_to_st_point, st_transform(temp_geometry, crs = 4326), sparse=FALSE)
  lon_lat_with_locaion[, temp_province] <- temp_result[, 1]
  
}

final_data <- lon_lat_with_locaion %>% bind_cols(tmax_obs_tibble)

grid_tmax <- final_data %>% filter(china_boundary == TRUE) %>% 
  select(-c( "china_boundary")) %>% 
  pivot_longer(cols = 3:36, names_to = "province", values_to = "contain") %>% 
  filter(contain == TRUE) %>% select(-contain)

grid_tmax_obs <- grid_tmax %>% pivot_longer(cols = 3:15332, names_to = "date", values_to = "obs")
grid_tmax_obs$obs <- grid_tmax_obs$obs-273.15

write_excel_csv(grid_tmax_obs,'/Volumes/Grace_PSSD/data analysis/hist-cal/grid_tmax_obs.csv')


##########################################################################################
## REORGANIZE THE SIMULATED TEMPERATURE SERIES IN FACTUAL AND CONTERFACTUAL SCENARIOS
# "hist" FILE FOLDER INCLUDES ALL THE SIMULATED TEMPERATURE DATA OF FACTUAL SCENARIO 
# "nat" FILE FOLDER INCLUDES ALL THE SIMULATED TEMPERATURE DATA OF COUNTERFACTUAL SCENARIO 
#  THE PROCESS IS SIMILAR TO THAT ABOVE
# 
## Simulated temperature data from models we used for scenarios
rm(list=ls())
gc()
ncintmax1 <- "ACCESS-ESM1-5_19790101-20201231.nc"
ncintmax2 <- "BCC-CSM2-MR_19790101-20201231.nc"
ncintmax3 <-"CanESM5_19790101-20201231.nc"
ncintmax4 <- "CNRM-CM6-1_19790101-20201231.nc"
ncintmax5 <- "FGOALS-g3_19790101-20201231.nc"

ncintmax7 <-"GFDL-ESM4_19790101-20201231.nc"
ncintmax9 <- "IPSL-CM6A-LR_19790101-20201231.nc"
ncintmax10 <- "MIROC6_19790101-20201231.nc"
ncintmax11 <-"MRI-ESM2-0_19790101-20201231.nc"
ncintmax12 <- "NorESM2-LM_19790101-20201231.nc"

# Take one of the simulated temperature data for factual scenario for example
# Load the data
ncintmax <-  ncintmax1
ncintmax_hist <- nc_open(paste0('hist/',ncintmax))

lon <- ncvar_get(ncintmax_hist,"lon") #74.0-134.5
nlon <- dim(lon) 
lat <- ncvar_get(ncintmax_hist,"lat") #18.5-53.0
nlat <- dim(lat) 
time <- ncvar_get(ncintmax_hist,"time")
tmax_array_hist <- ncvar_get(ncintmax_hist,"tasmax")
dim(tmax_array_hist) # 122 70  15330  ##1979/1/1-2020/12/31=15330(each as 365)

# Reorganize the simulated data 
tmax_hist <- tmax_array_hist-273.15 
dim(tmax_hist)

tmax_hist_matrix <- matrix(nrow = nlon*nlat, ncol = (dim(tmax_hist)[3]))
for (yr in c(1:(dim(tmax_hist)[3]))) {
  tmax_hist_matrix[ , yr] <- c(tmax_hist[ , , yr])
}

colnames(tmax_hist_matrix) <- 1:15330
tmax_hist_tibble <- as.data.frame.matrix(tmax_hist_matrix)

country_shp <- sf::st_read("shp/全国数据.shp")
province_shp <- sf::st_read("shp/省级数据-带标签.shp")


all_lon_lat <- expand.grid("lon" = lon, "lat" = lat)
lon_lat_with_locaion <- data.frame(matrix(ncol = length(province_shp$ENAME)+1))
colnames(lon_lat_with_locaion) <- c(province_shp$SNAME, "china_boundary")
lon_lat_with_locaion <- all_lon_lat %>% bind_cols(lon_lat_with_locaion)
lon_lat_to_st_point <- st_sfc(lapply(1:nrow(all_lon_lat), 
                                     FUN = function(i){return(st_point(c(all_lon_lat$lon[i], all_lon_lat$lat[i])))}), 
                              crs = 4326)


for (temp_province in c(province_shp$SNAME, "china_boundary")) {
  if (temp_province == "china_boundary") {
    temp_geometry <- country_shp$geometry[1]
  } else {
    temp_geometry <- province_shp %>% filter(SNAME == temp_province) %>% pull(geometry) 
  }
  
  temp_result <- st_intersects(lon_lat_to_st_point, st_transform(temp_geometry, crs = 4326), sparse=FALSE)
  lon_lat_with_locaion[, temp_province] <- temp_result[, 1]
  
}

final_data <- lon_lat_with_locaion %>% bind_cols(tmax_hist_tibble)

grid_tmax <- final_data %>% filter(china_boundary == TRUE) %>% 
  select(-c( "china_boundary")) %>% 
  pivot_longer(cols = 3:36, names_to = "province", values_to = "contain") %>% 
  filter(contain == TRUE) %>% select(-contain)

grid_tmax_mod <- grid_tmax %>% pivot_longer(cols = 3:15332, names_to = "date", values_to = "model")

write_excel_csv(grid_tmax_mod,'hist-cal/1_hist_tmax_mod.csv')




