library(data.table)
library(fixest)
library(zoo)
library(ggplot2)
library(ggforce)
library(ggnewscale)
library(cowplot)
library(viridis)
library(stringr)
library(sf)
library(haven)
library(splines)
library(rnaturalearth)
library(rnaturalearthdata)
library(modelsummary)

rm(list=ls())
gc()

"%!in%" <- Negate("%in%")


# these bits are for tables
pstars <- function(ps){
  p_stars <- ifelse(ps<.01,"***",ifelse(ps<.05,"**",ifelse(ps<.1,"*","")))
  return(p_stars)
}

f1 <- function(x) format(round(x,3),big.mark=",")
f2 <- function(x) format(round(x,0),big.mark=",")
gm <- list(list("raw"="nobs","clean"="Obs.","fmt"=f2),
           list("raw"="r.squared","clean"="R2","fmt"=f1))


# load the data
load("data/climatecropsconflict_agri.RData")


# some minor data wranglings ----

## make sure there are no 'extra' years
climatecropsconflict_dt <- climatecropsconflict_dt[elnino_year%!in%c("1996","2024")]

## area in 10,000 ha and postharvest period (3 months)
climatecropsconflict_dt[,`:=`(area10=area/10000,ph3=ifelse(harv%in%1:3,1,0))]

## trend (for no particular reason)
climatecropsconflict_dt[,trend:=as.numeric(as.factor(yearmo))]


## total number of incidents by conflict category
climatecropsconflict_dt[,.(sum(incidents)),by=.(type,civilian_targeting)]



# Eq.2 ----

impact2 <- function(x){
  r <- feols(incidents~area10:oni_crop_prec:gs_tc:ph3+area10:oni_crop_prec:gs_tc:I(1-ph3)+prec+tmax | xy + country^yearmo, data=x,vcov=conley(500)~x+y)
  m <- x[area10>0 & gs_tc>0,.(area10=mean(area10),gs_tc=mean(gs_tc),incidents=mean(incidents))]
  s <- 100*m$area10*m$gs_tc/m$incidents
  h1_coef <- round(r$coeftable["area10:oni_crop_prec:gs_tc:ph3","Estimate"]*s,3)
  h1_se <- round(r$coeftable["area10:oni_crop_prec:gs_tc:ph3","Std. Error"]*s,3)
  h1_stars <- pstars(r$coeftable["area10:oni_crop_prec:gs_tc:ph3","Pr(>|t|)"])
  h1_est <- paste0(format(round(h1_coef,1),nsmall=1),h1_stars)
  h1_std <- paste0("(",format(round(h1_se,1),nsmall=1),")")
  return(list(descriptive=c(area_tc=round(m$area10,3),gs_tc=round(m$gs_tc,3),incidents_tc=round(m$incidents,3)),effect=c(h1_est,h1_std),output=c(h1_coef,h1_se)))
}

reg2_pc <- feols(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg2_pn <- feols(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)

reg2p_pc <- fepois(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg2p_pn <- fepois(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)


modelsummary(list(reg2_pc,reg2_pn,reg2p_pc,reg2p_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))


impact2(climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"])

impact2(climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"])




