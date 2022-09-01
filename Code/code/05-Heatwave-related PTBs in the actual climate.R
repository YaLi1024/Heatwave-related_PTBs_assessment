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
#  - 'readr' package version 1.4.0
#  - 'readxl' package version 1.3.1
########################################################################################
rm(list=ls())
library(ncdf4)
library(readr)
library(readxl)

########################################################################################
# 05.CALCULATION OF HEATWAVE-ATTRIBUTABLE PTBS IN WARM SEASONS DURING 2010-2020
#      OF CHINA IN THE ACTUAL CLIAMTE.
########################################################################################
# We first estimate the baseline preterm birth (PTB) number for each day during 
#     2010-2020 for each grid in China.
# By identifying if there was at least a heatwave between day(t-6) to day(t) 
#     for each day(t), we can know if premature infants in day (t) whose mother 
#     have experienced heatwaves in their last gestational week.
# If mothers have experienced heatwave in their last gestational week (day (t-6) - day t),
#     then a part of baseline PTBs on the day(t) can be attributed to heatwaves.
#     This proportion is the attributable fraction and computed as ((RR-1))/RR. 
# 
# LOAD REVELANT DATA
#  -global_demographics_1950_2020.nc, providing global gridded population number. 
#  -grid_HWdays_May.csv,providing the longitude and latitude for grids of China.
#  -1979-2020nation_birth_rate.xlsx,providing provincial birth rates in China.
#  -nation_preterm_rate.csv, estimated by findings of previous studies and providing
#     the provincial PTB rates in China during 2010-2020.
#  -birth propotion by month.xlsx, providing the proportion of living birth in each 
#     month to that in a year.
#  -point estimate of RR, relative risk of PTB due to the exposure of heatwaves in the
#     last gestational week, obtained from the work (under review) of our research team.
#  -

############################################################################################
##  ESTIMATE THE BASELINE GRIDDED PTB NUMBER FOR EACH DAY IN WARM SEASONS DURING 2010-2020   

## 1) Load global population number during 1950-2020 ######################################
  ncpop<-nc_open("global_demographics_1950_2020.nc")
  
  lon <- ncvar_get(ncpop,"longitude")
   dim(lon)
  lat <- ncvar_get(ncpop,"latitude")
   dim(lat)
   # the data include the population stratified by 14 age groups
  age <- ncvar_get(ncpop,"age_band_lower_bound")
   dim(age)
  time <- ncvar_get(ncpop, "year")
  which(time==2010) # the index of 2010 
  pop_totals_array <- ncvar_get(ncpop,"demographic_totals")
  dim(pop_totals_array) 
  # calculate and store the population from all age groups
  pop_sumage_array <-  array(NA, dim = c(71,720,290), dimnames = list(time,lon,lat))
  for ( i in 1:dim(pop_totals_array)[1] ) {
    for ( j in 1:dim(pop_totals_array)[3]){
      for ( k in 1:dim(pop_totals_array)[4]){
        pop_sumage_array[i,j,k] <-sum(pop_totals_array[i,,j,k], na.rm=TRUE)
      }
    }
  }

## 2) Obtain gridded population number during 2010-2020 for China #########################
# obtain the longitude and latitude,then extract the data for China
  grid_hwdays_data = read_csv("grid_HWdays_May.csv")
  grid_pop_china=grid_hwdays_data
  
   for (i in 1: dim(grid_hwdays_data)[1]){
     for(j in 61:71){ # extract years from 2010-2020 
       longitude <- which(lon==grid_hwdays_data[i,1])
       latitude <- which(lat==grid_hwdays_data[i,2])
        grid_pop_china[i,(j-60+2)]<-pop_sumage_array[j,longitude,latitude]
     }
   }
  head(grid_pop_china) 

## 3) Calculate the annual living birth number during 2010-2020 for China ###################
# 1979-2020nation_birth_rate.xlsx provides the provincial birth rates (‰,per thousand).
# We assumed that the birth rate of specific grid is equal to birth rate of the province 
#     that the grid belongs to.
  nation_birth_rate0 <- read_excel("1979-2020nation_birth_rate.xlsx")
  
  nation_birth_rate <- as.data.frame(lapply(nation_birth_rate0,as.numeric))
  nation_birth_rate <- nation_birth_rate[1:11,-2]

  grid_birth_rate <- grid_pop_china
  # match the provincial birth rates to grids
  for (i in 1:dim(grid_birth_rate)[1]){
    for (j in 1:11){
      grid_birth_rate[i,(j+2)] <- (nation_birth_rate[(12-j), grid_birth_rate[[i,14]]])/1000 
    }
  }
  
  # calculate the grid birth number per year during 2010-2020
  grid_birth_number <- grid_pop_china
  
  for (i in 1:dim(grid_birth_number)[1]){
    for (j in 1:11){
      grid_birth_number[i,(j+2)] <- grid_birth_rate[i,(j+2)] * (grid_pop_china[i,(j+2)]) 
      }
  }
  
  # sum the grid birth number to provinces and the whole nation
  province_birth_number <- as.data.frame(matrix(NA, nrow=32, ncol=12))
  colnames(province_birth_number) <- c("province",paste0("y",2010:2020))
  
  province_birth_number$province <- c("nation", levels(prov))
  
  for (i in 1:11){
    province_birth_number[2:32,i+1] <-tapply(grid_birth_number[,i+2], grid_birth_number$province,sum,na.rm=TRUE)
    province_birth_number[1,i+1] <- sum(province_birth_number[-1,i+1])
  }
  
   write_excel_csv(province_birth_number,"province_birth_number.csv")
  

## 4) Calculate the baseline PTB birth number in warm seasons during 2010-2020 for China ##############
  # estimate the annual PTB birth number at the grid level
  nation_preterm_rate <- read_csv("nation_preterm_rate.csv")
  nation_preterm_rate <- as.data.frame(nation_preterm_rate)
  rownames(nation_preterm_rate) <- nation_preterm_rate$location
  
  grid_PTB_number <- grid_birth_number

  for (i in 1:dim(grid_PTB_number)[1]) {
    province <- grid_PTB_number[i,14]
    grid_PTB_number[i,3:13] <- grid_birth_number[i,3:13] * nation_preterm_rate[province, 2]
  }
  # the annual PTB birth number at the provincial and national level
  province_PTB_number = province_birth_number

  for (i in 1:11){
    province_PTB_number[2:32,i+1] <-tapply(grid_PTB_number[,i+2], grid_PTB_number$province,sum,na.rm=TRUE)
    province_PTB_number[1,i+1] <- sum(province_PTB_number[-1,i+1])
  }
  
  write_excel_csv(province_PTB_number,"province_PTB_number.csv")
  write_excel_csv(grid_PTB_number,"grid_PTB_number.csv")

 # estimate the PTB birth number in warm seasons per year at the grid level
  birth_proportion_by_month <- read_excel("birth propotion by month.xlsx")
  birth_proportion_by_month <- as.data.frame(birth_proportion_by_month[,1:7])
  rownames(birth_proportion_by_month) <-birth_proportion_by_month$location
  
  grid_PTB_month <- matrix(NA, nrow=23052, ncol=15) # 23052 = grid number(3842) * warm months(6)
  grid_PTB_month <- as.data.frame(grid_PTB_month)
  colnames(grid_PTB_month) <- c(colnames(grid_PTB_number), "month")
  grid_PTB_month[,15] <- 5:10
  
    
  for (i in 1:23052){
    grid_PTB_month[i,1:14] <- grid_PTB_number[trunc((i+5)/6),1:14]
    
    province <- grid_PTB_month[i,14]
    month <- as.character(grid_PTB_month[i,15])
    
    grid_PTB_month[i,3:13] <- grid_PTB_month[i,3:13] * birth_proportion_by_month[province, month]
  }
 
   write_excel_csv(grid_PTB_month,"grid_PTB_month.csv")
  
  
  # estimate the PTB birth number in warm seasons per year at the provincial and national level
  province_PTB_warm = province_PTB_number
  
  for (i in 1:11){
    province_PTB_warm[2:32,i+1] <-tapply(grid_PTB_month[,i+2], grid_PTB_month$province,sum,na.rm=TRUE)
    province_PTB_warm[1,i+1] <- sum(province_PTB_warm[-1,i+1])
  }
  
  write_excel_csv(province_PTB_warm,"province_PTB_warm.csv")


#################################################################################################
# ESTIMATE THE HEATWAVE-RELATED PTBS BASED ON THE BASELINE PTB NUMBER ON EACH DAY FOR EACH GRID 
# We have dertermined whether premature infants in each day (t) whose mother have experienced 
#     heatwaves in their last gestational week in "02-days with pre1week heatwave experience.R"
# grid_with_hwpre1w_sum.2.csv provides the information of heatwave exposures of mothers who 
#     delivered prematurely in each day (corresponding to baseline daily PTBs we estimated above).
# We used the RR,including coef and se from the model, from the work (under review) of our reseach 
#     team, and based on the Monte Carlo simulation to obtain the samples of coefficients
  
  ## 1) Create data frame to store results based on the data frame used above.  
  grid_aPTB_month = grid_PTB_month
  grid_aPTB_month_ci.l = grid_PTB_month
  grid_aPTB_month_ci.u = grid_PTB_month
  grid_with_hwpre1w_sum.2 <- read_csv("grid_with_hwpre1w_sum.2.csv")

  ## 2) Monte Carlo simulation to obtain the samples of coefficients
  set.seed(1234)
  nsim = 1000
  # coef = 0.172, se = 0.04345
  coefsim_m <- rnorm(nsim, mean = 0.172, sd = 0.04345)
  RR_estimates = exp(coefsim_m)
  quantile(RR_estimates,0.025)
  
  ## 3) based in the Monte Carlo simulations to calculate the gridded attributable PTBs and 95%CI 
  #     for each warm-season month during 2010-2020
  for (i in 1:nrow(grid_PTB_month)){
    for (j in 3:13){
      grid_aPTB_estimates = NA
      for (k in 1:nsim){
        RR = RR_estimates[k]
        AF = (RR-1)/RR
        #assume that there are 31 days in each month to simplify the calculation
        est <- (grid_PTB_month[i,j]/31)*grid_with_hwpre1w_sum.2[i,j]*AF
        grid_aPTB_estimates[k] <- as.numeric(est)
      }
      grid_aPTB_month[i,j] <- quantile(grid_aPTB_estimates,0.5)
      grid_aPTB_month_ci.l[i,j] <- quantile(grid_aPTB_estimates,0.025)
      grid_aPTB_month_ci.u[i,j] <- quantile(grid_aPTB_estimates,0.975)
    }
  }
    
  # 4) calculate the total of attributable PTBs in warm seasons per year for each grid
  grid_aPTB_number <- as.data.frame(grid_PTB_number)
  grid_aPTB_number_ci.l <- grid_aPTB_number
  grid_aPTB_number_ci.u <- grid_aPTB_number
  
  level = rep(1:3842,each=6) # 3842 grid number, 6 warm-season months
  
  for (i in 1:11){
    grid_aPTB_number[,i+2] <- tapply(grid_aPTB_month[,i+2],level,sum)
  }
  for (i in 1:11){
    grid_aPTB_number_ci.l[,i+2] <- tapply(grid_aPTB_month_ci.l[,i+2],level,sum)
  }
  for (i in 1:11){
    grid_aPTB_number_ci.u[,i+2] <- tapply(grid_aPTB_month_ci.u[,i+2],level,sum)
  }
  
  write_excel_csv(grid_aPTB_number,"grid_aPTB_number.csv")
  write_excel_csv(grid_aPTB_number_ci.l,"grid_aPTB_number_cil.csv")
  write_excel_csv(grid_aPTB_number_ci.u,"grid_aPTB_number_ciu.csv")
  
  
  # 5) calculate the number of attributable PTBs in warm seasons per year for provinces and nation
  province_aPTB = province_PTB_number
  province_aPTB_ci.l = province_aPTB
  province_aPTB_ci.u = province_aPTB
  
  for (i in 1:11){
    province_aPTB[2:32,i+1] <-tapply(grid_aPTB_month[,i+2], grid_aPTB_month$province,sum,na.rm=TRUE)
    province_aPTB[1,i+1] <- sum(province_aPTB[-1,i+1])
  }
  
  for (i in 1:11){
    province_aPTB_ci.l[2:32,i+1] <-tapply(grid_aPTB_number_ci.l[,i+2], grid_aPTB_number_ci.l$province,sum,na.rm=TRUE)
    province_aPTB_ci.l[1,i+1] <- sum(province_aPTB_ci.l[-1,i+1])
  }
  
  for (i in 1:11){
    province_aPTB_ci.u[2:32,i+1] <-tapply(grid_aPTB_number_ci.u[,i+2], grid_aPTB_number_ci.u$province,sum,na.rm=TRUE)
    province_aPTB_ci.u[1,i+1] <- sum(province_aPTB_ci.u[-1,i+1])
  }
  
  write_excel_csv(province_aPTB,"province_aPTB.csv")
  write_excel_csv(province_aPTB_ci.l,"province_aPTB_cil.csv")
  write_excel_csv(province_aPTB_ci.u,"province_aPTB_ciu.csv")
  

  
  