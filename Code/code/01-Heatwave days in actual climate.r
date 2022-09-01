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
#  - 'ncdf4' package version 1.17
#  - 'dplyr' package version  1.0.6
#  - 'lattice' package version 0.20.44
#  - 'readr' package version 1.4.0
#  - 'doParallel' package version 1.0.16
#  - 'foreach' package version 1.5.1
#  - 'sf' package version 0.9.8
#  - 'tidyverse' package version 1.3.1
#  - 'hrbrthemes' package version 0.8.0
#  - 'abind' package version 1.4.5
##########################################################################
rm(list=ls())
gc()
library(ncdf4)
library(dplyr)
library(lattice)
library(readr)
library(doParallel)

library(foreach)
library(sf)
library(tidyverse)
library(hrbrthemes)
library(abind)

setwd("/Volumes/Grace_PSSD/analysis/data/")
getwd()
##########################################################################
# 01 IDENTIFICATION OF HEATWAVE DAYS IN THE ACTUAL CLIMATE
################################################################################
# Use the observed daily maximum temperature series in warm seasons during 
#     1979-2020 from ERA5 datasets.
# Identify heatwaves based on the definition that the consecutive 2 or more
#     days with daily maximum temperature exceed the 90th distribution of daily 
#     maximum temperature series.
# We identify Heatwaves for each grid in China and calculate the heatwave days.
# 
#
# LOAD RELEVANT DATA
#  -ERA5_1979-2019_DailyX4_Tmax_daymax.nc, providing gridded daily maximum temperature
#     series during 1979-2019.
#  -Daily_Tmax_2020.nc, providing gridded daily maximum temperature series in 2020.
#  -shp file, providing the longitudes and latitudes at a 0.5
#     degree resolution for all provinces in China.
#     

##########################################################################
## LOAD THE TEMPERATURE SERIES AND STORE IN A ARRAY
# Load temperatures in 1979-2019
ncintmax <- nc_open("era5/ERA5_1979-2019_DailyX4_Tmax_daymax.nc")

lon <- ncvar_get(ncintmax,"lon")
nlon <- dim(lon)
lat <- ncvar_get(ncintmax,"lat")
nlat <- dim(lat)
time <- ncvar_get(ncintmax,"time")
tmax_array1 <- ncvar_get(ncintmax,"mx2t_NON_CDM")
dim(tmax_array1) # 141 111  14975  ##1979/1/1-2019/12/31=14975


ncintmax3 <- nc_open("era5/Daily_Tmax_2020.nc")
tmax_array3 <- ncvar_get(ncintmax3,"Tmax")

# Combine temperatures in 2020
tmax_array2 <- tmax_array3
for (i in 1:nlon) {
  for (k in 1:366) { # 366 days in 2020
    for (j in 1:nlat) {
      tmax_array2[i,(112-j),k] <- tmax_array3[i,j,k]
    }
  }
}

tmax_array <- abind(tmax_array1,tmax_array2,along = 3)


##########################################################################
##EXTRACT TEMPERATURES IN WARM SEASONS (MONTH 5-10)
Date <- as.Date(NA)
day <- dim(tmax_array)[3]
a <- as.Date('1978-12-31"')
Date <- 1:day + a

year <- as.numeric(format(Date,format="%Y"))
month <- as.numeric(format(Date,format="%m"))
warm_month <- month %in% 5:10

tmax_array_warm <- tmax_array[,,warm_month] 
dim(tmax_array_warm)

## CALCULATE THE GRIDDED THRESHODS BASED ON THE TEMPERATURE DISTRIBUTION IN WARM SEASONS OF 1985-2005
threshold_year <- year %in% c(1985:2005) 
# extract the daily temperatures in warm seasons between 1985 and 2005
tmax_array_threswarm <- tmax_array[,,threshold_year&warm_month]
dim(tmax_array_threswarm)
# calculate the threshold for each grid (the 90th of the cell temperature distribution)
a <- matrix(nrow=nlon,ncol=nlat,dimnames = list(lon,lat))

for (i in 1:nlon) {
  for (j in 1:nlat) {
    a[i,j] <- quantile(tmax_array_threswarm[i,j, ],0.90)
  }
}
a


####################################################################
## IDENTIFY HEATWAVES FOR EACH GRID 
# Define a function that if there are consecutive two or more days with daily maximum temperatures
#    exceed the cell threshold, then these days are labeled with 1 and deemed as heatwave days.
fun1 <- function(dat, threshold, i) {
  f = dat[i,,]
  for (j in 1:dim(dat)[2]) {
    d1 = (f[j, ] > (threshold[i,j]))
    d2 = c(FALSE, d1[-dim(dat)[3]])
    f[j, ] = d1 & d2 
    for (k in 2:dim(dat)[3]) {
      f[j, (k - 1)] = ifelse(f[j, k] == 1, 1, f[j, (k - 1)])
    }
  }
  return(f)
}

# Apply the function with ERA5 daily temperatures and thresholds we calculated
# We use parallel processing to increase efficiency.
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

ll <- foreach(m = 1:dim(tmax_array_warm)[1],.inorder = T) %dopar%  fun1(dat = tmax_array_warm,threshold = a,i = m)
stopCluster(cl)

# Store the results from parallel processing
result = tmax_array_warm
for (i in 1:dim(tmax_array_warm)[1]) {
  result[i,,]=ll[[i]]
}

dim(result)

##########################################################################
## CALCULATE HEATWAVE DAYS FOR EACH YEAR AT THE PROVINCIAL AND NATIONAL LEVEL.
# Sum the heatwave days for each year at the grid level.
year <- (1979:2020)
nyear <- 42 # total 42 years in the period of 1979-2020
heatwavedays_array <- array(NA,dim = c(nlon,nlat,nyear),
                            dimnames = list(lon,lat,year))
for (p in 1:42) {
  p_days <- (184*(p-1)+1):(184*p) 
  for (i in 1:nlon) {
    for (j in 1:nlat) {
      heatwavedays_array[i,j,p] <- sum(result[i,j,p_days]) 
    }
  }
}

# Store the result in the form of data frame
HWdays_yr_matrix <- matrix(nrow = nlon*nlat, ncol = (dim(heatwavedays_array)[3])) 
for (yr in c(1:(dim(heatwavedays_array)[3]))) {
  HWdays_yr_matrix[ ,yr] <- c(heatwavedays_array[ , , yr]) 
}
colnames(HWdays_yr_matrix) <- paste0("y",c(1979:2020))
HWdays_yr_tibble <- as.data.frame.matrix(HWdays_yr_matrix)   

# Fill in the longitude and latitude information from shp files
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


final_data <- lon_lat_with_locaion %>% bind_cols(HWdays_yr_tibble)

grid_HWdays_data <- final_data %>% filter(china_boundary == TRUE) %>% 
  select(-c( "china_boundary")) %>% 
  pivot_longer(cols = 3:36, names_to = "province", values_to = "contain") %>% 
  filter(contain == TRUE) %>% select(-contain)


## Calculate the yearly heatwave days for provinces in China 
province_HWdays_data_temp <- final_data %>% filter(china_boundary == TRUE) %>% 
  select(-c("lon", "lat", "china_boundary")) %>% 
  pivot_longer(cols = 1:34, names_to = "province", values_to = "contain") %>% 
  filter(contain == TRUE) %>% select(-contain) %>%
  pivot_longer(cols = -c("province"))
names(province_HWdays_data_temp) <- c("province","year","days")


province_HWdays_data <- province_HWdays_data_temp %>% group_by(province,year) %>% 
  dplyr::summarise(n = n(), mean = mean(days,na.rm = TRUE))%>% 
  pivot_wider(names_from = year, values_from = mean)


## Calculate the yearly heatwave days for the whole nation
nation_HWdays_data <- province_HWdays_data_temp %>% group_by(province,year) %>% 
  dplyr::summarise(n = n(), mean = mean(days,na.rm = TRUE)) %>% 
  group_by(year) %>% dplyr::summarise(mean = mean(mean,na.rm = TRUE))


write_excel_csv(grid_HWdays_data,"grid_HWdays_data.csv")
write_excel_csv(province_HWdays_data,"province_HWdays_data.csv",row.names = FALSE)
write_excel_csv(nation_HWdays_data,"nation_HWdays_data.csv",row.names = FALSE)

