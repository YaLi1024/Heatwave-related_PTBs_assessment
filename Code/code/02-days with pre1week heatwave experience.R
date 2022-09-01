######################################################################################
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

##############################################################################################
# 02 DETERMINE IF THERE WAS AT LEAST ONE HEATWAVE IN THE WEEK BEFORE EACH DAY (t)
##############################################################################################
# For further calculation of attributable PTBs in the actual climate, we need to determine if 
#    there was at least one heatwave event in the week before each day t during 2010-2020, 
#    namely during the time between day (t-6) and day (t).
#
#
# LOAD RELEVANT DATA
#   - file "result" obtained in "01-Heatwave days in the actual climate.R" (line 152)
#   - "01grid_HWdays_May.csv" used for extract the longitude and latitude of China
#   
#################################################################################################
# 1) Determine if there was at least one heatwave in the week before each day during 2010-2020

# extract the longitude and latitude of China from X01grid_HWdays_May and create an array to store
X01grid_HWdays_May <- read_csv("grid_HWdays_May.csv")
lon_china <- as.data.frame(X01grid_HWdays_May)[,1]
lat_china <- as.data.frame(X01grid_HWdays_May)[,2]
lon_china <- sort(unique(lon_china))
lat_china <- sort(unique(lat_china))


year <- paste0("y",2010:2020)
day <- c(paste0("May",1:31),paste0("Jun",1:30),paste0("Jul",1:31),paste0("Aug",1:31),
         paste0("Sep",1:30),paste0("Oct",1:31))

gchina_hw_warms <- array(NA,dim = c(122,70,11,184),
                                dimnames = list(lon_china,lat_china,year,day))


# determine if there was at least one heatwave in the week before each day based on the "result" in 
#   "01-Heatwave days in the actual climate.R" (line 152)
m2o <- as.data.frame(read_csv("m2o.csv"))

for (p in 1:11) {
  gchina_hw_warms[ , ,p, ] <- result[which(lon %in% lon_china),which(lat %in% lat_china),(m2o[1,p]:m2o[2,p])]  
} # to extract result in the warm season 

gchina_dayhw_ident <- gchina_hw_warms
gchina_dayhw_identYN <- gchina_dayhw_ident

 for(i in 1:122){
   for(j in 1:70){
     for(k in 1:11){
       for(l in 1:184){
         
         ifelse(l<7,gchina_dayhw_ident[i,j,k,l] <- sum(gchina_hw_warms[i,j,k,(1:l)]),
         gchina_dayhw_ident[i,j,k,l] <- sum(gchina_hw_warms[i,j,k,(l-6):l]) )
         
        if (l>=7) {   
          if(gchina_dayhw_ident[i,j,k,l]==0|gchina_dayhw_ident[i,j,k,l]==1){
          gchina_dayhw_identYN[i,j,k,l]<- 0}
          else if (gchina_dayhw_ident[i,j,k,l]==2){ if (gchina_hw_warms[i,j,k,l]==1 & gchina_hw_warms[i,j,k,(l-6)]==1){gchina_dayhw_identYN[i,j,k,l]<- 0}
                                                    else gchina_dayhw_identYN[i,j,k,l]<- 1}

          else gchina_dayhw_identYN[i,j,k,l] <-1
         }
         
         else if (l<=2){
           if(gchina_dayhw_ident[i,j,k,l]==2) {gchina_dayhw_identYN[i,j,k,l]<- 1}
           else {gchina_dayhw_identYN[i,j,k,l]<- 0}
         }
         
         else {
           if(gchina_dayhw_ident[i,j,k,l]==0|gchina_dayhw_ident[i,j,k,l]==1){
             gchina_dayhw_identYN[i,j,k,l]<- 0}
           
           else if (gchina_dayhw_ident[i,j,k,l]==2){ if (gchina_hw_warms[i,j,k,l]==1 & gchina_hw_warms[i,j,k,1]==1) gchina_dayhw_identYN[i,j,k,l]<- 0
                                                     else gchina_dayhw_identYN[i,j,k,l]<- 1}
           
           else gchina_dayhw_identYN[i,j,k,l] <-1
         }
         
       }
     }
   }
 }
#CHECK IT: RIGHT!
gchina_dayhw_identYN[which(lon_china==119.5),which(lat_china==23.0),1,62:92]
gchina_dayhw_ident[which(lon_china==100),which(lat_china==22.0),1,32:61]
gchina_hw_warms[which(lon_china==100),which(lat_china==22.0),1,32:61]
result[ which(lon==100.5),which(lat==22.0),(m2o[7,1]:m2o[8,1])] 

######################################################################################
# 2) Calculate the total number of days for which there was at least one heatwave in 
#       the week before in each month, which would be used for the the calculation of 
#       attributable PTB in the acutal climate.

month <-c("May","Jun","Jul","Aug","Sep","Oct")
gchina_with_hwpre1w_sum<- array(NA,dim = c(122,70,11,6),
                                dimnames = list(lon_china,lat_china,year,month))
day_location <- paste0(month[1],1:31)
for(i in 1:122){
  for(j in 1:70){
    for(k in 1:11){
      
      for(l in 1:6){
        if(l %in% c(2,5)) day_location <-paste0(month[l],1:30) else day_location <- paste0(month[l],1:31)
        gchina_with_hwpre1w_sum[i,j,k,l] <-sum(gchina_dayhw_identYN[i,j,k,day_location])
      }
      
    }
  }
}
#CHECK IT: RIGHT!
gchina_with_hwpre1w_sum[which(lon_china==101.5),which(lat_china==22.0),4,1:2]
l=2
if(l %in% c(2,5)) day_location <-paste0(month[l],1:30) else day_location <- paste0(month[l],1:31)
sum(gchina_dayhw_identYN[which(lon_china==101.5),which(lat_china==22.0),4,day_location])
gchina_dayhw_identYN[which(lon_china==101.5),which(lat_china==22.0),4,32:61]
result[ which(lon==101.5),which(lat==22.0),(m2o[5,4]:m2o[6,4])] 


######################################################################################
# 2) Store the result in form of data frame.
grid_HWdays_yr_month <- read_csv("grid_HWdays_yr_month.csv")
grid_with_hwpre1w_sum <- as.data.frame(grid_HWdays_yr_month)


for(i in 1:23052){
  for(k in 1:11){
    lon_location <- grid_with_hwpre1w_sum[i,1]
    lat_location <- grid_with_hwpre1w_sum[i,2]
    # the months in the gchina_with_hwpre1w_sum were represented as 1-6, which were 5-10 in grid_with_hwpre1w_sum
    month <-grid_with_hwpre1w_sum[i,15]-4  
    grid_with_hwpre1w_sum[i,(k+2)] <- gchina_with_hwpre1w_sum[which(lon_china==lon_location),which(lat_china==lat_location),k,month]
  }
}

write_excel_csv(grid_with_hwpre1w_sum,"grid_with_hwpre1w_sum.2.csv")


