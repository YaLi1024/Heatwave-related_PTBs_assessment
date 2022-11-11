###################################################################################################
# Updated version of the R code for the analysis in:
#
#  "The burden of heatwave-related preterm births and associated human capital losses in China"
# 
#
# Update: 11 Nov. 2022
# * an updated version of this code is available at:
#   https://github.com/YaLi1024/Heatwave-related_PTBs_assessment
#
# Requirements:
#  - R version 4.1.0
#  - 'tidyverse' package version 1.3.1
#  - 'readr' package version 1.4.0
#  - 'doParallel' package version 1.0.16
#  - 'foreach' package version 1.5.1
##################################################################################################
rm(list=ls())
gc()
library(tidyverse)
library(readr)
library(doParallel)
library(foreach)

##########################################################################################
## 03 CALCULATE THE HEATWAVE-RELATED PTBS IN THE FACTUAL AND COUNTERFACTUAL SCENARIOS AND
##    THOSE CAUSED BY ANTHROPOGENIC CLIMATE CHANGE
##########################################################################################
# LOAD RELEVANT DATA
#   - grid_PTB_month.csv calculated and saved in "02-Heatwave-related PTBs in the actual climate.R"
#   - Data generated from "012-013-calibration and identification for each model.R"
#      1. ten 'X cal hist_gridhwpre1w_month.csv' stored in 'hist-cal/pre1w/', X is the num of models 
#      ranging from 1 to 10
#      2. ten 'X cal nat_gridhwpre1w_month.csv' stored in 'nat-cal/pre1w/', X is the num of models 
#      ranging from 1 to 10

#####################################################################################################
# CALCULATE THE ATTRIBUTABLE PTBS

files=list.files('hist-cal/',pattern="*.csv")

grid_PTB_month <- read_csv('grid_PTB_month.csv')
# sort the data 
grid_PTB_month <- as_tibble(grid_PTB_month[,c(1,2,14,15,3:13)]) %>% arrange(lon,lat,province,month)
grid_pre1w_aPTB <- grid_PTB_month
grid_pre1w_aPTB_ci.l <- grid_PTB_month
grid_pre1w_aPTB_ci.u <- grid_PTB_month


# Monte carlo simulations
set.seed(1234)
nsim = 1000
# coef = 0.172, se = 0.04345
coefsim_m <- rnorm(nsim, mean = 0.172, sd = 0.04345)
RR_estimates = exp(coefsim_m)
quantile(RR_estimates,0.025)

# 1) Create a function to calculate the attributable PTBs and obtain the 95%CIs from the 2.5th and
#    97.5th percentiles of the empirical distribution across coefficients samples and scenario models.
fun <- function(dat,pathway,province,i){
  library(tidyverse)
  library(readr)
  
  f = as.data.frame(matrix(NA,nrow=3,ncol =16),dimnames= c(colname(dat),'est'))
  f[,1:4] <- dat[i,1:4]
  f[,16] <- c('point','ci.l','ci.u')
  
  files=list.files(pathway,pattern="*.csv")
  for (j in 5:15){
    
    estimates <- NA
    for (k in 1:length(files)){
      s=files[k]
      grid_pre1w_mod <- read_csv(paste0(pathway,s))
      grid_pre1w_mod <- as_tibble(grid_pre1w_mod[grid_pre1w_mod$province ==province,]) %>% arrange(lon,lat,province,month)
      est0 <- (dat[i,j]/31)*grid_pre1w_mod[i,j]
      
      for (l in 1:nsim){
        RR = RR_estimates[l]
        AF = (RR-1)/RR
        est <- est0*AF
        estimates[(k-1)*nsim+l] <- as.numeric(est)
      }
    }
    f[1,j] = quantile(estimates,0.5) 
    f[2,j] <- quantile(estimates,0.025)
    f[3,j] <- quantile(estimates,0.975)
  } 
  
  return(f)
}

## 2) Calculate the results for the factual scenario
cl <- makeCluster(no_cores)
registerDoParallel(cl)

system.time(
  f <- foreach(m=1:dim(grid_PTB_month)[1], .inorder = T) %dopar% fun(dat = grid_PTB_month,pathway='hist-cal/',province=prov,i = m)
)
stopCluster(cl)

# transform into the format of data frame
grid_attriPTB_month <-f[[1]]
for (i in 2:dim(grid_PTB_test)[1]){
  grid_PTB_month <- rbind(f[[i]],grid_attriPTB_month)
}
grid_attriPTB_month
colnames(grid_attriPTB_month) <- c('lon','lat','province','month',paste0('y',2010:2020),"est")

write_excel_csv(grid_attriPTB_month,"hist_aPTB/grid_aPTB_hist.csv")

## 3) Calculate the results for the counterfactual scenario
cl <- makeCluster(no_cores)
registerDoParallel(cl)

system.time(
  f <- foreach(m=1:dim(grid_PTB_month)[1], .inorder = T) %dopar% fun(dat = grid_PTB_month,pathway='nat-cal/',province=prov,i = m)
)
stopCluster(cl)

# transform into the format of data frame
grid_attriPTB_month <-f[[1]]
for (i in 2:dim(grid_PTB_test)[1]){
  grid_PTB_month <- rbind(f[[i]],grid_attriPTB_month)
}
grid_attriPTB_month
colnames(grid_attriPTB_month) <- c('lon','lat','province','month',paste0('y',2010:2020),"est")

write_excel_csv(grid_attriPTB_month,"hist_aPTB/grid_aPTB_nat.csv")


## 4) Calculate the attributable PTBs caused by anthropogenic climate change
#  The differences of heatwave-related PTBs between the factual and counterfactual scenarios
#     are the attributable PTBs caused by anthropogenic climate change.






