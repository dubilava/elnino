# Relationship between ENSO and economic growth
## Original code: Christopher Callahan 
## Christopher.W.Callahan.GR@dartmouth.edu
## Augmentation: David ubilava
## david.ubilava@sydney.edu.au

rm(list=ls())
library(ggplot2)
library(margins)
library(sandwich)
library(stargazer)
library(lfe)
library(dplyr)
library(lemon)
library(lspline)
library(fixest)
library(cowplot)
library(texreg)
library(data.table) #DU
library(fixest) #DU
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

# location
# setwd("/Volumes/rc/lab/C/CMIG/ccallahan/Variability_Economics/Replication/Scripts/") #Data/")
loc_panel <- "CallahanMankin_ENSOEconomics-main/Data/Panel/"
loc_save <- "CallahanMankin_ENSOEconomics-main/Data/RegressionResults/"

# helper functions
source("CallahanMankin_ENSOEconomics-main/Scripts/HelperFunctions.R")

# read in data
y1 <- 1960
y2 <- 2019
panel <- read.csv(paste(loc_panel,"ENSO_Growth_Panel_",y1,"-",y2,".csv",sep=""))

panel_dt <- data.table(panel) #DU

## tau index, which is abs(cor(e,p))+abs(cor(e,t)) where cor() is the
## the max of the three-month running mean of correlations between the
## DJF enso index and monthly precipitation and temperature (respectively)
## during the Jun[t-1]--Aug[t] 15-month period

tau_dt <- panel_dt[year==2010,.(iso,t_p_corr_e)]
tau_dt <- tau_dt[complete.cases(tau_dt)]
tau_dt

ggplot(tau_dt,aes(x=t_p_corr_e))+
  geom_histogram(color="lightgray",fill="dimgray")+
  geom_vline(xintercept=c(mean(tau_dt$t_p_corr_e),quantile(tau_dt$t_p_corr_e,.5)),color=c("indianred","steelblue"),linewidth=1)


panel_dt[,`:=`(e_0=e,e_1=shift(e,1),e_2=shift(e,2),e_3=shift(e,3),e_4=shift(e,4),e_5=shift(e,5)),by=.(iso)] #DU

panel_dt[,`:=`(de_1=e_1-e_0,de_2=e_2-e_0,de_3=e_3-e_0,de_4=e_4-e_0,de_5=e_5-e_0),by=.(iso)] #DU

re <- feols(growth~e_0+de_1+de_2+de_3+de_4+de_5+(e_0+de_1+de_2+de_3+de_4+de_5):t_p_corr_e | iso, data=panel_dt) #DU
summary(re,vcov="cluster") #DU

t0 <- re$coeftable["e_0","Estimate"]
t1 <- re$coeftable["e_0:t_p_corr_e","Estimate"]

tau_dt$effect <- t0+tau_dt$t_p_corr_e*t1

ggplot(tau_dt,aes(x=effect))+
  geom_histogram(color="lightgray",fill="dimgray")


world <- merge(world,tau_dt,by.x="adm0_a3",by.y="iso",all=T)

ggplot(data = world) +
  geom_sf(aes(fill=effect))+
  scale_fill_viridis_c(option="B")


re <- feols(growth~e_0+e_1+e_2+e_3+e_4+e_5+(e_0+e_1+e_2+e_3+e_4+e_5):t_p_corr_e | iso,data=panel_dt[year%in%c(1960:1997,1998:2019)]) #DU
summary(re,vcov="cluster") #DU

re <- feols(growth~e_0+(e_0):t_p_corr_e | iso,data=panel_dt) #DU
summary(re,vcov="cluster") #DU

re <- feols(growth~e_0+(e_0):t_p_corr_e | iso,data=panel_dt[year%in%c(1960:1997,1999:2019)]) #DU
summary(re,vcov="cluster") #DU

# i <- i+1
# ggplot(panel_dt[iso==unique(iso)[i]],aes(x=year,y=growth))+
#   geom_line(na.rm=T)+
#   labs(title=unique(panel_dt$iso)[i])

sum(re$coeftable[7:12,1]) #DU

oni <- fread("https://psl.noaa.gov/data/correlation/oni.data",skip=1) #DU

oni_djf <- oni[,.(year=V1,oni=V2)] #DU

ggplot(oni_djf,aes(x=year,y=oni))+
  geom_line()


e_dt <- panel_dt[iso==unique(iso)[1],.(year,e,c,nino3,nino34,oni)] #DU
e_lg <- melt(e_dt,id.vars = "year") #DU

ggplot(e_lg,aes(x=year,y=value,color=variable))+ #DU
  geom_line() #DU



panel_dt <- merge(panel_dt,oni_djf,by="year",all.x=T) #DU

panel_dt[,`:=`(e_0=oni,e_1=shift(oni,1),e_2=shift(oni,2),e_3=shift(oni,3),e_4=shift(oni,4),e_5=shift(oni,5)),by=.(iso)] #DU

re <- feols(growth~e_0+e_1+e_2+e_3+e_4+e_5+(e_0+e_1+e_2+e_3+e_4+e_5):t_p_corr_e | iso,data=panel_dt) #DU
summary(re,vcov=conley(2000,distance="spherical")) #DU

sum(re$coeftable[7:12,1]) #DU

# panel_dt[,`:=`(growth_1=shift(growth,1))] #DU
# 
# summary(feols(growth~e_0+e_1+e_2+e_3+e_4+e_5+(e_0+e_1+e_2+e_3+e_4+e_5):t_p_corr_e+growth_1 | iso,data=panel_dt)) #DU
