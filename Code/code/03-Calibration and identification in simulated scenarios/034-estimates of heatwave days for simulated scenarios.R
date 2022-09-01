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
##########################################################################################################
library(readr)
library(dplyr)
library(tidyverse)
library(abind)

###########################################################################################################
# 03. CALIBRATE THE SIMULATED TEMPERATURES, IDENTIFY HEATWAVES FOR EACH MODEL AND CALCULATE THE HEATWAVE 
#     DAYS FOR FACTUAL AND COUNTERFACTUAL SCENARIOS AND THOSE INDUCED BY ANTHROPOGENIC CLIMATE CHANGE.
#
# 034.With the results obtained in "032-033-calibration and identification for each model.R", we calculated  
#     heatwave days for 1) the factual scenario, 2) counterfactual scenario and 3) those induced by
#     anthropogenic climate change.
#
# LOAD RELEVANT DATA AND FILE
#  -Results obtained from "032-033-calibration and identification for each model.R"
#     1. ten files named "X cal hist_gridhwpre1w_month.csv" stored in folder "hist-cal/42y_hwdays/"; X is num of 1 to 10
#     2. ten files named ""X cal nat_gridhwpre1w_month.csv" stored in folder "nat-cal/42y_hwdays/"; X is num of 1 to 10 


############################################################################################################
## 1) HEATWAVE DAYS FOR THE FACTUAL SCENARIO
# Mean heatwave days of ten models for factual scenario were deemed as the final results for this scenario
# Put results of ten models into a data frame
path = 'hist-cal/42yr_hwdays/'
files=list.files(path,pattern="*.csv")
n=length(files)
grid_HW_data_temp <- read_csv("1_cal_hist_42yrgrid_HWdays")

grid_HW_data_temp <- as_tibble(grid_HW_data_temp) %>% pivot_longer(cols = 4:45, names_to = "year", values_to = "days")
grid_HW_data_temp <- grid_HW_data_temp[,-5]

for (i in 1:n){
  s=files[i]
  grid_hw_data <- as_tibble(read_csv(paste0(path,s))) %>% 
    pivot_longer(cols = 4:45, names_to = "year", values_to = "days")
  
  grid_HW_data_temp <- merge(grid_HW_data_temp,grid_hw_data, by=c("lon","lat","province","year"))  
  rm(grid_hw_data)
}
names(grid_HW_data_temp)[5:14] <- paste0("model_",1:10)

write_excel_csv(grid_HW_data_temp,"42y_hist_grid_hwdays_allmodelcal.csv")


# Mean heatwave days of ten models for each grid of China 
grid_mean_hwdata_temp <-as_tibble(grid_HW_data_temp)%>%
  pivot_longer(cols=5:14, names_to="model",values_to="days")%>%
  group_by(lon,lat,province,year) %>% 
  dplyr::summarise(mean1=mean(days))

# Mean heatwave days from 10 models for each grids, further used to calculate mean heatwave days  
#   of grids in a province to obtain the results at provincial level
province_mean_hwdata <- grid_mean_hwdata_temp  %>%
  group_by(province,year) %>% dplyr::summarise(mean2=mean(mean1),na.rm = TRUE)%>%
  select(-c('na.rm'))%>% pivot_wider(names_from = year, values_from = mean2)


# Mean heatwave days of all grids in a nation for each model, which can further be used to calculate  
#   the mean and 95%CI for heatwave days at the national level (Figure 1.b)
hist_grid_hwdays <-as_tibble(grid_HW_data_temp)%>%
  pivot_longer(cols=5:14, names_to="model",values_to="days")%>%
  group_by(province,year,model) %>% dplyr::summarise(meantimes=mean(days),na.rm = TRUE)%>%
  group_by(year,model) %>% dplyr::summarise(days=mean(meantimes),na.rm = TRUE)%>%
  select(-c('na.rm'))

write.csv(hist_grid_hwdays,"42y_hist_nation_hwdays_models.csv",row.names = FALSE)


############################################################################################################
## 2) HEATWAVE DAYS FOR THE COUNTERFACTUAL SCENARIO
# The process is same to that for the facutal scenarios
# We just need to change the path = 'nat-cal/42yr_hwdays/'
# Two results are saved: "42y_nat_grid_hwdays_allmodelcal.csv";"42y_nat_nation_hwdays_models.csv"

############################################################################################################
## 3) HEATWAVE DAYS INDUCED BY ANTHROPOGENIC CLIMATE CHANGE.

# Load results saved in process 1) and 2)
X42y_hist_grid_hwdays_allmodelcal <-read_csv('42y_hist_nation_hwdays_models.csv')
X42y_nat_grid_hwdays_allmodelcal <- read_csv('42y_nat_nation_hwdays_models.csv')

# transform the format and put them into a data frame.  
hist_grid_hwday_model <- X42y_hist_grid_hwdays_allmodelcal%>%
  pivot_longer(cols=5:14, names_to="model",values_to="hist_days")

nat_grid_hwday_model <- X42y_nat_grid_hwdays_allmodelcal%>%
  pivot_longer(cols=5:14, names_to="model",values_to="nat_days")

grid_HWdays_diff <- hist_grid_hwday_model %>% merge(nat_grid_hwday_model, by=c('lon','lat','province','year','model'))

# Calculate the heatwave days induced by anthropogenic climate change.
# first obtain the difference between the results of paired models for factual and counterfactual scenarios
# then obtain the mean of the difference by models and the yearly mean heatwave days
# finally obtain the mean heatwave days in a decade at the grid level to draw Figure 1.d
grid_HWdays_diff$decade <- NA
grid_HWdays_diff$decade[which(grid_HWdays_diff$year %in% 1979:1989)] <- '1979-1989'
grid_HWdays_diff$decade[which(grid_HWdays_diff$year %in% 1990:1999)] <- '1990-1999'
grid_HWdays_diff$decade[which(grid_HWdays_diff$year %in% 2000:2009)] <- '2000-2009'
grid_HWdays_diff$decade[which(grid_HWdays_diff$year %in% 2010:2020)] <- '2010-2020'

grid_HWdays_diff <- grid_HWdays_diff %>%
  group_by(lon,lat,province,year,model,decade) %>%
  dplyr::summarise(HumanInduced_hwdays1=(hist_days-nat_days)) %>%
  group_by(lon,lat,province,year,decade)%>%
  dplyr::summarise(HumanInduced_hwdays2=mean(HumanInduced_hwdays1)) %>%  
  group_by(lon, lat,province,decade)%>% 
  dplyr::summarise(HInduced_hwdays_decade=mean(HumanInduced_hwdays2))

write.csv(grid_HWdays_diff,"grid_HWdays_diff.csv",row.names = FALSE)



