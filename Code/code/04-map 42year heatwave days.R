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
#  - 'sf' package version 1.3.1
#  - 'tidyverse' package version 1.3.1
#  - 'readr' package version 1.4.0
#  - 'raster' package version  3.4.13
#  - 'RColorBrewer' package version  1.1.2
#  - 'ggtext' package version  0.1.1
#  - 'cowplot' package version  1.1.1
#  - 'showtext' package version  0.9.4
#  - 'hrbrthemes' package version  0.8.0
#  - 'ncdf4' package version  1.17
#  - 'rgdal' package version 1.5.23
#  - 'maptools' package version 1.1.1
#  - 'ggplot2' package version 3.3.3
##########################################################################
rm(list=ls())
gc()

library(sf)
library(tidyverse)
library(readr)
library(raster)
library(RColorBrewer)
library(ggtext)
library(cowplot)
library(showtext)
library(hrbrthemes)
library(ncdf4)
library(rgdal)
library(maptools)
library(ggplot2)
font_add("Tms", "/Library/Fonts/Times New Roman.ttf")
showtext_auto()

setwd("/Users/grace/DATA/")
getwd()

#######################################################################################
# 04. CODE FOR FIGURE 1
#
#
# LOAD REVELANT DATA
#  -FIGURE 1.a,c: nation_HWdays_data.csv & grid_HWdays_data.csv derived from 
#      "01-Heatwave days in actual climate.R"
#  -FIGURE 1.b: 42y_hist_nation_hwdays_models.csv & 42y_hist_nation_hwdays_models.csv
#      from "034-estimates of hwdays for simulated scenarios.R"
#  -FIGURE 1.d: grid_HWdays_diff.csv from "034-estimates of hwdays for simulated scenarios.R"
#  -file folder 'gridmap' for drawing the map at grid level
######################################################################################
## DRAW THE LINE GRAPH AND MAP FOR HEATWAVE DAYS IN THE ACTUAL CLIMATE

## Draw the line graph for heatwaves in the actual climate (Figure 1.a)
nation_hwdays <- read_csv("nation_HWdays_data.csv")
nation_hwdays$year <- as.numeric(substr(nation_hwdays$year, 2,5))

summary(nation_hwdays$mean)

breakx2 <- c(1980,1990,2000,2010,2020)
breaky2 <- c(0,5,10,15,20,25,30)

hw_line <- ggplot(nation_hwdays,aes(x=year,y=mean,group=1))+
  stat_summary(fun = 'mean',geom = 'line',size=0.5, color='red')+
  geom_smooth(method="lm",fill=NA, size=0.5, linetype=2,color='red')+
  labs(x="Year",y="National heatwave days per year")+
  scale_x_continuous(breaks=breakx2, labels =c('1980','1990','2000','2010','2020'))+
  scale_y_continuous(breaks=breaky2, labels =c('0','5',"10",'15',"20","25","30"))+
  theme(text = element_text(colour = "black",family = "Tms"), 
        panel.grid = element_blank(),
        panel.background = element_rect(fill = '#ffffff', color = NA),
        axis.line = element_line(size=0.5,color = 'grey10')
  ) +
  theme(legend.text = element_text(size = rel(1),family = 'Tms'),legend.title =element_blank()) +
  theme(aspect.ratio = 0.6,
        axis.text.x=element_text(color="black",size = 16),
        axis.text.y= element_text(color="black",size = 16),
        axis.title =element_text(color="black",size = 16))+
  coord_cartesian(ylim=c(0,30))

hw_line

ggsave("/Users/grace/data/biye/line_hwdays_90th.pdf", width =9, height = 6, device = NULL)


## Draw the map for heatwaves in the actual climate (Figure 1.c)
#    Preprocess the data
yr42_grid_HWdays_data <- read_csv('grid_HWdays_data.csv')

decade_grid_HWdays_data <- yr42_grid_HWdays_data %>% pivot_longer(cols = 3:44, names_to = "year", values_to = "days")
summary(decade_grid_HWdays_data$days)

year <- decade_grid_HWdays_data$year
decade_grid_HWdays_data$decade <- NA
decade_grid_HWdays_data$decade[which(year %in% paste0('y',1979:1989))] <- '1979-1989'
decade_grid_HWdays_data$decade[which(year %in% paste0('y',1990:1999))] <- '1990-1999'
decade_grid_HWdays_data$decade[which(year %in% paste0('y',2000:2009))] <- '2000-2009'
decade_grid_HWdays_data$decade[which(year%in% paste0('y',2010:2020))] <- '2010-2020'
decade_grid_HWdays_data <- decade_grid_HWdays_data %>% group_by(lon,lat,province,decade)%>% dplyr::summarise(meandays=mean(days))


# Draw the map 
CHN_sp = readOGR("gridmap/grid map/省/31个省面.shp",use_iconv=TRUE, encoding = "UTF-8")
CHNshp = fortify(CHN_sp)

china_border <- readOGR('gridmap/grid map/国界线及岛屿/bou2_4l.shp') %>%
  spTransform('+proj=longlat +datum=WGS84 +no_defs') %>%
  fortify()

china_border2 <- china_border[china_border$lat>15, ]
breakx <- c(80,90,100,110,120,130)
breaky <- c(10,20,30,40,50)

my_colormap <- colorRampPalette(brewer.pal(11,'YlOrRd'))(32)
my_colormap <-c("#FFF5B5", "#FEE288", "#FEC25D","#FD913E","#FC542B","#AF0026")


decade <- decade_grid_HWdays_data$decade
grid_HWdays_data <- decade_grid_HWdays_data[which(decade == '2010-2020'),-4]
summary(grid_HWdays_data)

grid_HWdays_data$value <- cut(grid_HWdays_data$meandays, breaks=c(0,18,21,24,27,30,80),
                              labels=c('0~18','18~21','21~24','24~27','27~30','30~80'),right=FALSE,order=TRUE)

plot1 <- ggplot() + 
  theme(plot.margin=unit(rep(1,4),'cm'))+
  geom_raster(data=grid_HWdays_data, aes(x=lon, y=lat, fill=value))+ #, na.color = "White"
  scale_fill_manual(values = my_colormap,name="Human-induced \nHeatwave days", na.value="grey")+
  #scale_fill_gradientn(colors=my_colormap,name='Human-induced\nheatwave days',breaks=c(0,2,4,6,8,10,12))+
  #scale_fill_brewer(palette = 'YlOrRd')+
  coord_quickmap()+

  geom_path(aes(x = long, y = lat, group = group), color = 'grey20',size=0.15,
            data = china_border2) +theme_classic()+

  #facet_wrap(~decade,ncol = 2)+
  theme(text = element_text(colour = "black",family = "Tms"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = '#ffffff', color = NA),
        legend.background = element_rect(fill = '#ffffff', color = NA),
        axis.line = element_blank())+
  theme(legend.title = element_text(size = rel(1.2),vjust = 0.75,hjust = 0.5,family = 'Tms')) +
  theme(legend.text = element_text(size = rel(1.0),family = 'Tms',)) +
  theme(legend.key.size = unit(0.5,"cm"),
        legend.spacing = unit(1.0,"cm"),
        panel.border = element_rect(fill=NA,color="grey20",linetype=1,size=0.8),
        #legend.position = "bottom",legend.direction = "horizontal",
  ) +
  theme(strip.background.x = element_blank(),strip.background.y = element_blank())+
  theme(strip.text = element_text(size = 11,family = 'Tms',colour = "black")) +
  theme(strip.text = element_blank()) +
  theme(plot.margin = unit(c(5,5,5,5),"mm"))+
  theme(axis.text.x=element_text(color="black",size = 12),
        axis.text.y= element_text(color="black",size =12),
        axis.title =element_blank() )+
  scale_x_continuous(breaks=breakx, labels = paste(breakx,"°E",sep = ""))+
  scale_y_continuous(breaks=breaky, labels = paste(breaky,"°N",sep = ""))

plot1
ggsave("map_hwdays.eps", width = 16, height =16, device = cairo_ps)

######################################################################################
## DRAW THE LINE GRAPH AND MAP FOR HEATWAVE DAYS IN THE SIMULATED SCENARIOS 
## The process to draw Figure 1.b,d is similar to that in factual climate (above)
##
## Relevant data for FIGURE 1.b (42y_hist_nation_hwdays_models.csv & 
##    42y_hist_nation_hwdays_models.csv) and FIGURE 1.d (grid_HWdays_diff.csv) are 
##    all the results stored in "034-estimates of hwdays for simulated scenarios.R"

