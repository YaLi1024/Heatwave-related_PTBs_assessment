####################################################################################################
# Updated version of the R code for the analysis in:
#
#  "Heatwave exposure and preterm birth in China: attributable disease burden and
#   human capital consequences"
#
# Update: 30 Aug. 2022
# * an updated version of this code is available at:
#   https://github.com/gasparrini/2019_vicedo-cabrera_Epidem_Rcodedata
#
# Requirements:
#  - R version 4.1.0
#  - 'readr' package version 1.4.0
#  - 'dplyr' package version  1.0.6
#  - 'tidyverse' package version 1.3.1
#  - 'abind' package version 1.4.5
##########################################################################
library(readr)
library(dplyr)
library(tidyverse)
library(abind)

###############################################################################
# 03. CALIBRATE THE SIMULATED TEMPERATURES, IDENTIFY HEATWAVES FOR EACH MODEL
#     AND CALCULATE THE HEATWAVE DAYS FOR FACTUAL AND COUNTERFACTUAL SCENARIOS
#
# 032.With the data obtained in "031-data format conversion.R", then we calibrate  
#     the simulated temperatures with observed temperatures
# 033.Identify the heatwaves based on the calibrated data of each model for the scenarios
#    -The identification include 1) the heatwaves during 1979-2020; 2) heatwaves during 2010-2020
#     to 3) determine if there was at least one heatwave in the week before each day (t), which would
#     be used for the the calculation of attributable PTB in simulated scenarios.
#
# LOAD RELEVANT DATA AND FILE
#  -Data obtained from "031-data format conversion.R"
#     1. grid_tmax_obs.csv; 
#     2. ten files named "x_hist_tmax_mod.csv" stored in "hist" file folder;X is the num of 1 to 10  
#     3. ten files named "x_nat_tmax_mod.csv"stored in "nat" file folder;X is the num of 1 to 10 
#  -Method file from the "Hands-on Tutorial on a Modeling Framework for Projections 
#     of Climate Change Impacts on Health" : fhempel_main.R;fhempel_corr.R
#

##########################################################################
## CALIBRATE SIMULATED TEMPERATURES SEREIES BASED ON THE OBSERVED TEMPERATURES

##Generate the series of date during 1979-2020
Date <- as.Date(NA)
for (i in 1:42){ 
  a <- as.Date(paste0(as.character(1977+i),'-12-31"'))
  Date[(365*i-364):(365*i)] <- (1:365) + a
}

# Use the data obtained in "031-data format conversion.R"
tmax_obs <- as_tibble(read_csv('grid_tmax_obs.csv'))
tmax_obs <- as.data.frame(tmax_obs %>% group_by(lon,lat) %>% mutate(Date)%>%arrange(lon, lat, Date))

files1=list.files('/data/hist/',pattern="*.csv")
files2=list.files('/data/nat/',pattern="*.csv")

# There are 1O models for each scenario;n refers to the num of the model which ranges from 1 to 10;
for (n in 1:10){
  s1=files1[n] 
  s2=files2[n]
  grid_mod_hist <- read_csv(s1)
  grid_mod_hist <- as.data.frame(grid_mod_hist %>% group_by(lon,lat) %>% mutate(Date)%>%arrange(lon, lat, Date))
  grid_mod_nat <- read_csv(s2)
  grid_mod_nat <- as.data.frame(grid_mod_nat %>% group_by(lon,lat) %>% mutate(Date)%>%arrange(lon, lat, Date))
  
  grid_mod_hist$bias <- NA
  grid_mod_nat$bias<- NA

  source("/data/hist/fhempel_main.R")
  source("/data/hist/fhempel_corr.R")
 
  # calibrate the simulated daily temperature series during 1979-2020 for each grid
  # 3854 refers to the number of grids in China
  # 15330 refers to the total days which equal to 365 (days) *42 (years)
  # "grid_mod_hist" and "grid_mod_nat" are paired results for factual and counterfactual scenarios respectively.
  for (i in 1:3854){ 
    obs <- tmax_obs[(i*15330-15329): (i*15330),c(6,5)]
    mod_hist <- mod_hist[(i*15330-15329): (i*15330),c(6,5)]
    mod_nat <- mod_nat[(i*15330-15329): (i*15330),c(6,5)]
    corr <- fhempel_main(obs,mod_hist,output = "correction")
    
    grid_mod_hist$bias[(i*15330-15329): (i*15330)] <- fhempel_corr(mod_hist,corr)[,2]
    grid_mod_nat$bias[(i*15330-15329): (i*15330)] <- fhempel_corr(mod_nat,corr)[,2]
  }
  
#################################################################################################
### USE THE CALIBRATED TEMPERATURES TO IDENTIFY  HEATWAVES IN THE FACTUAL SCENARIO
  # The process are similar to that in actual climate;
  # Some detailed explanation of the code can be seen in "01-Heatwave days in actual climate.R"
  
## 1) Identify heatwaves for the facutal scenario during 1979-2020  
  grid_tmax_mod <- grid_mod_hist
  
  month <- as.numeric(format(grid_tmax_mod$Date,format="%m"))
  molcal_warm <- grid_tmax_mod[month %in% c(5,6,7,8,9,10),]
  molcal_warm$hwday <- NA
  molcal_warm$model <- NULL
  molcal_warm$date <- NULL
  
  rm(grid_tmax_mod)
  gc()
  
  # Threshold for heatwave identification during 1979-2020,based on the temperature distribution in warm 
  #   season during 1986-2005
  threshold1 <- NA
  for (i in 1:3854){
    grid_cal_warm <- molcal_warm[(7728*i-7727):(7728*i),]
    threshold1[i]  <- quantile(grid_cal_warm$bias[(184*7+1):(184*27)],0.90) 
    gridday0 <- ifelse (grid_cal_warm$bias > threshold1[i], 1, 0)
    hwday <- NA
    for (j in 1:42){ # total 42 year during 1979-2020
      gridday1 <- gridday0[(184*j-183):(184*j)] # 184 days in the warm season per year 
      gridday2 <- c(0,gridday1[-184])
      gridday <- NA 
      gridday[1] <- 0
      for (k in 2:184){
        gridday[k] <- ifelse(gridday1[k]+gridday2[k]==2,1,0)
        if (gridday[k]==1) gridday[(k-1)]=1
      }
      hwday[(184*j-183):(184*j)] <-gridday
    }
    molcal_warm$hwday[(7728*i-7727):(7728*i)] <- hwday 
  }
 
  # Sum for the yearly heatwave days for each grid
  year <- substr(as.character(molcal_warm$Date),1,4)
  molcal_warm$year <-year
  
  hist_hwday_modcal<- as_tibble(molcal_warm) %>% group_by(lon,lat,province,year) %>% summarise(yr_hwday=sum(hwday))%>%
    select(lon,lat,province,yr_hwday,year)%>% pivot_wider(names_from = year, values_from = yr_hwday)
  
  write_excel_csv(hist_hwday_modcal,paste0('hist-cal/42y_hwdays/',substr(s,1,2),'cal hist_42yrgrid_HWdays.csv'))
  
  rm(hist_hwday_modcal)
  gc()
  
  
## 2) Identify heatwaves for the factual scenario during 2010-2020
  molcal_warm_11yr <- molcal_warm[year %in% as.character(2010:2020),]
  molcal_warm_11yr$hwday <- NA
  
  rm(molcal_warm)
  gc()
  # threshold for heatwave identification during 2010-2020
  threshold2 <- NA
  for (i in 1:3854){ # 3854 grids in China
    grid_cal_warm <- molcal_warm_11yr[(2024*i-2023):(2024*i),] # total 2024 (184*11) days in warm season during 2010-2020 
    threshold2[i]  <- quantile(grid_cal_warm$bias,0.90) 
    gridday0 <- ifelse (grid_cal_warm$bias > threshold2[i], 1, 0)
    hwday <- NA
    for (j in 1:11){
      gridday1 <- gridday0[(184*j-183):(184*j)] 
      gridday2 <- c(0,gridday1[-184]) 
      gridday <- NA
      gridday[1] <- 0
      for (k in 2:184){
        gridday[k] <- ifelse(gridday1[k]+gridday2[k]==2,1,0)
        if (gridday[k]==1) gridday[(k-1)]=1
      }
      hwday[(184*j-183):(184*j)] <-gridday
    }
    molcal_warm_11yr$hwday[(2024*i-2023):(2024*i)] <- hwday 
  }
 
  
## 3) determine if there was at least one heatwave in the week before each day (t)
  molcal_warm_11yr$hwday_pre1w <- NA
  
  for (i in 1:3854){
    grid_cal_warm <- molcal_warm_11yr[(2024*i-2023):(2024*i),]
    hwday_pre1w <- NA
    for(j in 1:11){
      a <- grid_cal_warm$hwday[(184*j-183):(184*j)] #extract the warm season
      b <- c(0,a[-184])
      c <- (a & b) 
      pre1w <- NA
      for (k in 1:184){
        ident_temp <- ifelse(k < 7,sum(c[1:k]), sum(c[(k-6):k])) 
        pre1w[k] <- ifelse(ident_temp==0,0,1)
      }
      hwday_pre1w[(184*j-183):(184*j)] <- pre1w
    }
    molcal_warm_11yr$hwday_pre1w[(2024*i-2023):(2024*i)] <- hwday_pre1w
  }
  
  # total number of days for which there was at least a heatwave in the week before in each month, which would
  #       be used for the the calculation of attributable PTB in simulated scenarios.
  molcal_warm_11yr$month <- as.numeric(substr(molcal_warm_11yr$Date,6,7))
  
  grid_hwpre1w_sumcal<- as_tibble(molcal_warm_11yr) %>% group_by(lon,lat,province,year,month) %>%
    summarise(hwpre1w_sum=sum(hwday_pre1w))%>%select(lon,lat,province,year,month,hwpre1w_sum)%>%
    pivot_wider(names_from = year, values_from = hwpre1w_sum)
  
  write_excel_csv(grid_hwpre1w_sumcal,paste0('hist-cal/pre1w/',substr(s,1,2),'cal hist_gridhwpre1w_month.csv'))

  a<-ls()
  rm(list=a[which(a!='threshold1'& a!='threshold2'& a!='n'& a!='grid_mod_nat')])
  gc()

  
#################################################################################################
### USE THE CALIBRATED TEMPERATURES TO IDENTIFY  HEATWAVES IN THE COUNTERFACTUAL SCENARIO
## 1) Identify heatwaves for the facutal scenario during 1979-2020 
  grid_tmax_mod <- grid_mod_nat
    
  month <- as.numeric(format(grid_tmax_mod$Date,format="%m"))
  molcal_warm <- grid_tmax_mod[month %in% c(5,6,7,8,9,10),]
  molcal_warm$hwday <- NA
  molcal_warm$model <- NULL
  molcal_warm$date <- NULL
  
  rm(grid_tmax_mod)
  gc()
  
  molcal_warm <- as_tibble(molcal_warm) %>% arrange(lon, lat, Date)
  molcal_warm <- as.data.frame(molcal_warm)
  
  for (i in 1:3854){
    grid_cal_warm <- molcal_warm[(7728*i-7727):(7728*i),] 
    # The gridded thresholds for heatwave identification in the counterfactual scenario were 
    #     derived from that in the factual scenario.
    gridday0 <- ifelse (grid_cal_warm$bias > threshold1[i], 1, 0)
    hwday <- NA
    for (j in 1:42){
      gridday1 <- gridday0[(184*j-183):(184*j)] 
      gridday2 <- c(0,gridday1[-184]) 
      gridday <- NA 
      gridday[1] <- 0
      for (k in 2:184){
        gridday2[k] <- gridday1[(k-1)]
        gridday[k] <- ifelse(gridday1[k]+gridday2[k]==2,1,0)
        if (gridday[k]==1) gridday[(k-1)]=1
      }
      hwday[(184*j-183):(184*j)] <-gridday
    }
    molcal_warm$hwday[(7728*i-7727):(7728*i)] <- hwday 
  }
  
  # Sum for the yearly heatwave days for each grid
  year <- substr(as.character(molcal_warm$Date),1,4)
  molcal_warm$year <-year
  
  hist_hwday_modcal<- as_tibble(molcal_warm) %>% group_by(lon,lat,province,year) %>% summarise(yr_hwday=sum(hwday))%>%
    select(lon,lat,province,yr_hwday,year)%>% pivot_wider(names_from = year, values_from = yr_hwday)
  
  write_excel_csv(hist_hwday_modcal,paste0('nat-cal/42y_hwdays/',substr(s,1,2),'cal nat_42yrgrid_HWdays.csv'))
  
  rm(hist_hwday_modcal)
  gc()
  
## 2) Identify heatwaves for the factual scenario during 2010-2020
  molcal_warm_11yr <- molcal_warm[year %in% as.character(2010:2020),]
  molcal_warm_11yr$hwday <- NA
  
  rm(molcal_warm)
  gc()
  
  molcal_warm_11yr <- as_tibble(molcal_warm_11yr) %>% arrange(lon, lat, Date)
  molcal_warm_11yr <- as.data.frame(molcal_warm_11yr)
  
  for (i in 1:3854){
    grid_cal_warm <- molcal_warm_11yr[(2024*i-2023):(2024*i),] 
    # The gridded thresholds for heatwave identification in the counterfactual scenario were 
    #     derived from that in the factual scenario.
    gridday0 <- ifelse (grid_cal_warm$bias > threshold2[i], 1, 0)
    hwday <- NA
    for (j in 1:11){
      gridday1 <- gridday0[(184*j-183):(184*j)]
      gridday2 <- c(0,gridday1[-184])
      gridday <- NA 
      gridday[1] <- 0
      for (k in 2:184){
        gridday[k] <- ifelse(gridday1[k]+gridday2[k]==2,1,0)
        if (gridday[k]==1) gridday[(k-1)]=1
      }
      hwday[(184*j-183):(184*j)] <-gridday
    }
    molcal_warm_11yr$hwday[(2024*i-2023):(2024*i)] <- hwday 
  }
  
  
## 3) determine if there was at least one heatwave in the week before each day (t)
  molcal_warm_11yr$hwday_pre1w <- NA
  
  for (i in 1:3854){
    grid_cal_warm <- molcal_warm_11yr[(2024*i-2023):(2024*i),]
    hwday_pre1w <- NA
    for(j in 1:11){
      a <- grid_cal_warm$hwday[(184*j-183):(184*j)] 
      b <- c(0,a[-184])
      c <- (a & b) 
      pre1w <- NA
      for (k in 1:184){
        ident_temp <- ifelse(k < 7,sum(c[1:k]), sum(c[(k-6):k])) 
        pre1w[k] <- ifelse(ident_temp==0,0,1)
      }
      hwday_pre1w[(184*j-183):(184*j)] <- pre1w
    }
    molcal_warm_11yr$hwday_pre1w[(2024*i-2023):(2024*i)] <- hwday_pre1w
  }
  
  # total number of days for which there was at least a heatwave in the week before in each month, which would
  #       be used for the the calculation of attributable PTB in simulated scenarios.
  molcal_warm_11yr$month <- as.numeric(substr(molcal_warm_11yr$Date,6,7))
  
  grid_hwpre1w_sumcal<- as_tibble(molcal_warm_11yr) %>% group_by(lon,lat,province,year,month) %>%
    summarise(hwpre1w_sum=sum(hwday_pre1w))%>%select(lon,lat,province,year,month,hwpre1w_sum)%>%
    pivot_wider(names_from = year, values_from = hwpre1w_sum)
  
  write_excel_csv(grid_hwpre1w_sumcal,paste0('nat-cal/pre1w/',substr(s,1,2),'cal nat_gridhwpre1w_month.csv'))
  
  rm(list=ls())
  gc()
 
} 
  
  
  

