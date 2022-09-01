##########################################################################
# Updated version of the R code for the analysis in:
#
#  "Heatwave exposure and preterm birth in China: attributable disease burden and
#   human capital consequences"
# 
#
# Update: 30 Aug. 2022
# * an updated version of this code is available at:
#   https://github.com/YaLi1024/Heatwave-related_PTBs_assessment
#
########################################################################################
##07 CALCULATE THE HUMAN CAPITAL CONSEQUENCES OF ATTRIBUTABLE PTBS IN THE ACTUAL CLIMATE OR 
##  THOSE CAUSED BY ANTHROPOGENIC CLIMATE CHANGE
########################################################################################
# The data we used include:
#   - the prevalence of specific human capital (P_conseq) obtained from existing research
#   - national PTB rate (P_ptb) from existing research
#   - the relative risk of specific human capital outcome due to PTB (RR) which also obtained
#     from existing research  
#   - the number of attributable PTBs (nPTB) in the actual cliamte or those from anthropogenic 
#     climate change, which have been calculated in "05-Heatwave-related PTBs in the actual 
#     climate.R" and "06-heatwave-related PTBs in simualted scenarios.R".

# create a function to calculate the impacts on human capital indicators (without reduced IQ)
additional_HCs = function(P_conseq,P_ptb,RR,nPTB){
  I0 <- P_conseq/(P_ptb*RR + (1-P_ptb))
  I0 <- round(I0,digits = 3)
  Ie <- I0*RR 
  add_hcconseq <- nPTB*(Ie-P_conseq)
  return(add_hcconseq)
}


neonates_mor <- c(0.0083,0.0078,0.0069,0.0063,0.0059,0.0054,0.0049,0.0045,0.0039,0.0035,0.0035)
neonates_mor_mean <- mean(neonates_mor)

death_mean <-trunc(nPTB * (0.019 - neonates_mor_mean))
asthma_mean <- additional_HCs(0.0212,0.067,1.46,nPTB)
Diabetes_T1D_mean <- additional_HCs(0.00545,0.067,1.18,nPTB)
Diabetes_T2D_mean <- additional_HCs(0.1036,0.067,1.51,nPTB)
reduced_IQ <- nPTB*8.4
ASD_mean <-nPTB*(0.07-0.0039)
ADHD_mean <- additional_HCs(0.065,0.067,1.6,nPTB)

