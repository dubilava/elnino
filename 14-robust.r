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
load("data/climatecropsconflict.RData")


# some minor data wranglings ----

## make sure there are no 'extra' years
climatecropsconflict_dt <- climatecropsconflict_dt[elnino_year%!in%c("1996","2024")]

## area in 10,000 ha and postharvest period (3 months)
climatecropsconflict_dt[,`:=`(incidence=ifelse(incidents>0,1,0),area10=area/10000,ph3=ifelse(harv%in%1:3,1,0),pp3=ifelse(plant%in%1:3,1,0))]

## trend (for no particular reason)
climatecropsconflict_dt[,trend:=as.numeric(as.factor(yearmo))]


# R1 - fixed effects ----

## B1: year-month ----
impact3b1 <- function(x){
  r <- feols(incidents~area10:oni_crop:gs_tc:ph3+area10:oni_crop:gs_tc:I(1-ph3)+prec+tmax | xy + yearmo, data=x,vcov=conley(500)~x+y)
  m <- x[area10>0 & gs_tc>0,.(area10=mean(area10),gs_tc=mean(gs_tc),incidents=mean(incidents))]
  s <- 100*m$area10*m$gs_tc/m$incidents
  h1_coef <- round(r$coeftable["area10:oni_crop:gs_tc:ph3","Estimate"]*s,3)
  h1_se <- round(r$coeftable["area10:oni_crop:gs_tc:ph3","Std. Error"]*s,3)
  h1_stars <- pstars(r$coeftable["area10:oni_crop:gs_tc:ph3","Pr(>|t|)"])
  h1_est <- paste0(format(round(h1_coef,1),nsmall=1),h1_stars)
  h1_std <- paste0("(",format(round(h1_se,1),nsmall=1),")")
  return(list(descriptive=c(area_tc=round(m$area10,3),gs_tc=round(m$gs_tc,3),incidents_tc=round(m$incidents,3)),effect=c(h1_est,h1_std),output=c(h1_coef,h1_se)))
}

reg3_pc <- feols(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg3_pn <- feols(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)


modelsummary(list(reg3_pc,reg3_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))

impact3b1(climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"])

impact3b1(climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"])


## B2: country-year ----
impact3b2 <- function(x){
  r <- feols(incidents~area10:oni_crop:gs_tc:ph3+area10:oni_crop:gs_tc:I(1-ph3)+prec+tmax | xy + country^year, data=x,vcov=conley(500)~x+y)
  m <- x[area10>0 & gs_tc>0,.(area10=mean(area10),gs_tc=mean(gs_tc),incidents=mean(incidents))]
  s <- 100*m$area10*m$gs_tc/m$incidents
  h1_coef <- round(r$coeftable["area10:oni_crop:gs_tc:ph3","Estimate"]*s,3)
  h1_se <- round(r$coeftable["area10:oni_crop:gs_tc:ph3","Std. Error"]*s,3)
  h1_stars <- pstars(r$coeftable["area10:oni_crop:gs_tc:ph3","Pr(>|t|)"])
  h1_est <- paste0(format(round(h1_coef,1),nsmall=1),h1_stars)
  h1_std <- paste0("(",format(round(h1_se,1),nsmall=1),")")
  return(list(descriptive=c(area_tc=round(m$area10,3),gs_tc=round(m$gs_tc,3),incidents_tc=round(m$incidents,3)),effect=c(h1_est,h1_std),output=c(h1_coef,h1_se)))
}

reg3_pc <- feols(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^year, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg3_pn <- feols(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^year, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)

modelsummary(list(reg3_pc,reg3_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))

impact3b2(climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"])

impact3b2(climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"])


## B3: no weather ----
impact3b3 <- function(x){
  r <- feols(incidents~area10:oni_crop_prec:gs_tc:ph3+area10:oni_crop_prec:gs_tc:I(1-ph3) | xy + country^yearmo, data=x,vcov=conley(500)~x+y)
  m <- x[area10>0 & gs_tc>0,.(area10=mean(area10),gs_tc=mean(gs_tc),incidents=mean(incidents))]
  s <- 100*m$area10*m$gs_tc/m$incidents
  h1_coef <- round(r$coeftable["area10:oni_crop_prec:gs_tc:ph3","Estimate"]*s,3)
  h1_se <- round(r$coeftable["area10:oni_crop_prec:gs_tc:ph3","Std. Error"]*s,3)
  h1_stars <- pstars(r$coeftable["area10:oni_crop_prec:gs_tc:ph3","Pr(>|t|)"])
  h1_est <- paste0(format(round(h1_coef,1),nsmall=1),h1_stars)
  h1_std <- paste0("(",format(round(h1_se,1),nsmall=1),")")
  return(list(descriptive=c(area_tc=round(m$area10,3),gs_tc=round(m$gs_tc,3),incidents_tc=round(m$incidents,3)),effect=c(h1_est,h1_std),output=c(h1_coef,h1_se)))
}

reg3_pc <- feols(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3 | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg3_pn <- feols(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3 | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)

modelsummary(list(reg3_pc,reg3_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))

impact3b3(climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"])

impact3b3(climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"])


# R2 - exposure measures ----

climatecropsconflict_dt[,`:=`(area_dum=ifelse(area10>=.5,1,0),gs_tc_dum=ifelse(gs_tc>.33,1,0))]

## B4: area ----
impact3b4 <- function(x){
  r <- feols(incidents~area_dum:oni_crop:gs_tc:ph3+area_dum:oni_crop:gs_tc:I(1-ph3)+prec+tmax | xy + country^yearmo, data=x,vcov=conley(500)~x+y)
  m <- x[area_dum>0 & gs_tc>0,.(area_dum=mean(area_dum),gs_tc=mean(gs_tc),incidents=mean(incidents))]
  s <- 100*m$area_dum*m$gs_tc/m$incidents
  h1_coef <- round(r$coeftable["area_dum:oni_crop:gs_tc:ph3","Estimate"]*s,3)
  h1_se <- round(r$coeftable["area_dum:oni_crop:gs_tc:ph3","Std. Error"]*s,3)
  h1_stars <- pstars(r$coeftable["area_dum:oni_crop:gs_tc:ph3","Pr(>|t|)"])
  h1_est <- paste0(format(round(h1_coef,1),nsmall=1),h1_stars)
  h1_std <- paste0("(",format(round(h1_se,1),nsmall=1),")")
  return(list(descriptive=c(area_tc=round(m$area_dum,3),gs_tc=round(m$gs_tc,3),incidents_tc=round(m$incidents,3)),effect=c(h1_est,h1_std),output=c(h1_coef,h1_se)))
}

reg3_pc <- feols(incidents~area_dum:oni_crop:gs_tc+area_dum:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg3_pn <- feols(incidents~area_dum:oni_crop:gs_tc+area_dum:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)

modelsummary(list(reg3_pc,reg3_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))

impact3b4(climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"])

impact3b4(climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"])


## B5: teleconnection ----
impact3b5 <- function(x){
  r <- feols(incidents~area10:oni_crop:gs_tc_dum:ph3+area10:oni_crop:gs_tc_dum:I(1-ph3)+prec+tmax | xy + country^yearmo, data=x,vcov=conley(500)~x+y)
  m <- x[area10>0 & gs_tc_dum>0,.(area10=mean(area10),gs_tc_dum=mean(gs_tc_dum),incidents=mean(incidents))]
  s <- 100*m$area10*m$gs_tc_dum/m$incidents
  h1_coef <- round(r$coeftable["area10:oni_crop:gs_tc_dum:ph3","Estimate"]*s,3)
  h1_se <- round(r$coeftable["area10:oni_crop:gs_tc_dum:ph3","Std. Error"]*s,3)
  h1_stars <- pstars(r$coeftable["area10:oni_crop:gs_tc_dum:ph3","Pr(>|t|)"])
  h1_est <- paste0(format(round(h1_coef,1),nsmall=1),h1_stars)
  h1_std <- paste0("(",format(round(h1_se,1),nsmall=1),")")
  return(list(descriptive=c(area_tc=round(m$area10,3),gs_tc_dum=round(m$gs_tc_dum,3),incidents_tc=round(m$incidents,3)),effect=c(h1_est,h1_std),output=c(h1_coef,h1_se)))
}

reg3_pc <- feols(incidents~area10:oni_crop:gs_tc_dum+area10:oni_crop:gs_tc_dum:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg3_pn <- feols(incidents~area10:oni_crop:gs_tc_dum+area10:oni_crop:gs_tc_dum:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)

modelsummary(list(reg3_pc,reg3_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))

impact3b5(climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"])

impact3b5(climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"])


## B6: area & teleconnection ----
impact3b6 <- function(x){
  r <- feols(incidents~area_dum:oni_crop:gs_tc_dum:ph3+area_dum:oni_crop:gs_tc_dum:I(1-ph3)+prec+tmax | xy + country^yearmo, data=x,vcov=conley(500)~x+y)
  m <- x[area_dum>0 & gs_tc_dum>0,.(area_dum=mean(area_dum),gs_tc_dum=mean(gs_tc_dum),incidents=mean(incidents))]
  s <- 100*m$area_dum*m$gs_tc_dum/m$incidents
  h1_coef <- round(r$coeftable["area_dum:oni_crop:gs_tc_dum:ph3","Estimate"]*s,3)
  h1_se <- round(r$coeftable["area_dum:oni_crop:gs_tc_dum:ph3","Std. Error"]*s,3)
  h1_stars <- pstars(r$coeftable["area_dum:oni_crop:gs_tc_dum:ph3","Pr(>|t|)"])
  h1_est <- paste0(format(round(h1_coef,1),nsmall=1),h1_stars)
  h1_std <- paste0("(",format(round(h1_se,1),nsmall=1),")")
  return(list(descriptive=c(area_tc=round(m$area_dum,3),gs_tc_dum=round(m$gs_tc_dum,3),incidents_tc=round(m$incidents,3)),effect=c(h1_est,h1_std),output=c(h1_coef,h1_se)))
}

reg3_pc <- feols(incidents~area_dum:oni_crop:gs_tc_dum+area_dum:oni_crop:gs_tc_dum:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg3_pn <- feols(incidents~area_dum:oni_crop:gs_tc_dum+area_dum:oni_crop:gs_tc_dum:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)

modelsummary(list(reg3_pc,reg3_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))

impact3b6(climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"])

impact3b6(climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"])


# R3 - teleconnections ----

## B7: precipitation ----
impact3b7 <- function(x){
  r <- feols(incidents~area10:oni_crop:gs_tcp:ph3+area10:oni_crop:gs_tcp:I(1-ph3)+prec+tmax | xy + country^yearmo, data=x,vcov=conley(500)~x+y)
  m <- x[area10>0 & gs_tcp>0,.(area10=mean(area10),gs_tcp=mean(gs_tcp),incidents=mean(incidents))]
  s <- 100*m$area10*m$gs_tcp/m$incidents
  h1_coef <- round(r$coeftable["area10:oni_crop:gs_tcp:ph3","Estimate"]*s,3)
  h1_se <- round(r$coeftable["area10:oni_crop:gs_tcp:ph3","Std. Error"]*s,3)
  h1_stars <- pstars(r$coeftable["area10:oni_crop:gs_tcp:ph3","Pr(>|t|)"])
  h1_est <- paste0(format(round(h1_coef,1),nsmall=1),h1_stars)
  h1_std <- paste0("(",format(round(h1_se,1),nsmall=1),")")
  return(list(descriptive=c(area_tc=round(m$area10,3),gs_tcp=round(m$gs_tcp,3),incidents_tc=round(m$incidents,3)),effect=c(h1_est,h1_std),output=c(h1_coef,h1_se)))
}

reg3_pc <- feols(incidents~area10:oni_crop:gs_tcp+area10:oni_crop:gs_tcp:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg3_pn <- feols(incidents~area10:oni_crop:gs_tcp+area10:oni_crop:gs_tcp:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)

modelsummary(list(reg3_pc,reg3_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))

impact3b7(climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"])

impact3b7(climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"])


## B8: temperature ----
impact3b8 <- function(x){
  r <- feols(incidents~area10:oni_crop:gs_tct:ph3+area10:oni_crop:gs_tct:I(1-ph3)+prec+tmax | xy + country^yearmo, data=x,vcov=conley(500)~x+y)
  m <- x[area10>0 & gs_tct>0,.(area10=mean(area10),gs_tct=mean(gs_tct),incidents=mean(incidents))]
  s <- 100*m$area10*m$gs_tct/m$incidents
  h1_coef <- round(r$coeftable["area10:oni_crop:gs_tct:ph3","Estimate"]*s,3)
  h1_se <- round(r$coeftable["area10:oni_crop:gs_tct:ph3","Std. Error"]*s,3)
  h1_stars <- pstars(r$coeftable["area10:oni_crop:gs_tct:ph3","Pr(>|t|)"])
  h1_est <- paste0(format(round(h1_coef,1),nsmall=1),h1_stars)
  h1_std <- paste0("(",format(round(h1_se,1),nsmall=1),")")
  return(list(descriptive=c(area_tc=round(m$area10,3),gs_tct=round(m$gs_tct,3),incidents_tc=round(m$incidents,3)),effect=c(h1_est,h1_std),output=c(h1_coef,h1_se)))
}

reg3_pc <- feols(incidents~area10:oni_crop:gs_tct+area10:oni_crop:gs_tct:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg3_pn <- feols(incidents~area10:oni_crop:gs_tct+area10:oni_crop:gs_tct:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)

modelsummary(list(reg3_pc,reg3_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))

impact3b8(climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"])

impact3b8(climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"])


# R4 - conflict incidence ----

## B9: conflict incidence ----
impact3b9 <- function(x){
  r <- feols(incidence~area10:oni_crop:gs_tc:ph3+area10:oni_crop:gs_tc:I(1-ph3)+prec+tmax | xy + country^yearmo, data=x,vcov=conley(500)~x+y)
  m <- x[area10>0 & gs_tc>0,.(area10=mean(area10),gs_tc=mean(gs_tc),incidence=mean(incidence))]
  s <- 100*m$area10*m$gs_tc/m$incidence
  h1_coef <- round(r$coeftable["area10:oni_crop:gs_tc:ph3","Estimate"]*s,3)
  h1_se <- round(r$coeftable["area10:oni_crop:gs_tc:ph3","Std. Error"]*s,3)
  h1_stars <- pstars(r$coeftable["area10:oni_crop:gs_tc:ph3","Pr(>|t|)"])
  h1_est <- paste0(format(round(h1_coef,1),nsmall=1),h1_stars)
  h1_std <- paste0("(",format(round(h1_se,1),nsmall=1),")")
  return(list(descriptive=c(area_tc=round(m$area10,3),gs_tc=round(m$gs_tc,3),incidence_tc=round(m$incidence,3)),effect=c(h1_est,h1_std),output=c(h1_coef,h1_se)))
}

reg3_pc <- feols(incidence~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg3_pn <- feols(incidence~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)

modelsummary(list(reg3_pc,reg3_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))

impact3b9(climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"])

impact3b9(climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"])


# R5 - poisson regression ----

## B10: poisson regression ----
### coefficients
reg2p_pc <- fepois(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg2p_pn <- fepois(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)

modelsummary(list(reg2p_pc,reg2p_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))

### effect
reg2p_pc <- fepois(incidents~area10:oni_crop:gs_tc:I(1-ph3)+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg2p_pn <- fepois(incidents~area10:oni_crop:gs_tc:I(1-ph3)+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)

modelsummary(list(reg2p_pc,reg2p_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))


## B11: poisson year-month ----
### coefficients
reg2p_pc <- fepois(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg2p_pn <- fepois(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)

modelsummary(list(reg2p_pc,reg2p_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))

### effect
reg2p_pc <- fepois(incidents~area10:oni_crop:gs_tc:I(1-ph3)+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg2p_pn <- fepois(incidents~area10:oni_crop:gs_tc:I(1-ph3)+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)

modelsummary(list(reg2p_pc,reg2p_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))


# R6 - hotspots ----
impact3 <- function(x){
  r <- feols(incidents~area10:oni_crop:gs_tc:ph3+area10:oni_crop:gs_tc:I(1-ph3)+prec+tmax | xy + country^yearmo, data=x,vcov=conley(500)~x+y)
  m <- x[area10>0 & gs_tc>0,.(area10=mean(area10),gs_tc=mean(gs_tc),incidents=mean(incidents))]
  s <- 100*m$area10*m$gs_tc/m$incidents
  h1_coef <- round(r$coeftable["area10:oni_crop:gs_tc:ph3","Estimate"]*s,3)
  h1_se <- round(r$coeftable["area10:oni_crop:gs_tc:ph3","Std. Error"]*s,3)
  h1_stars <- pstars(r$coeftable["area10:oni_crop:gs_tc:ph3","Pr(>|t|)"])
  h1_est <- paste0(format(round(h1_coef,1),nsmall=1),h1_stars)
  h1_std <- paste0("(",format(round(h1_se,1),nsmall=1),")")
  return(list(descriptive=c(area_tc=round(m$area10,3),gs_tc=round(m$gs_tc,3),incidents_tc=round(m$incidents,3)),effect=c(h1_est,h1_std),output=c(h1_coef,h1_se)))
}

## cells
cells1_dt <- climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting",.(incidents=sum(incidents)),by=.(xy)]

cells1_dt <- cells1_dt[order(-incidents)]

one <- round(.01*nrow(cells1_dt))

excluded_cells <- cells1_dt[1:one]

xy01_dt <- cells1_dt[(one+1):nrow(cells1_dt)]
xy1 <- xy01_dt$xy


cells2_dt <- climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting",.(incidents=sum(incidents)),by=.(xy)]

cells2_dt <- cells2_dt[order(-incidents)]

one <- round(.01*nrow(cells2_dt))

excluded_cells <- cells2_dt[1:one]

xy02_dt <- cells2_dt[(one+1):nrow(cells2_dt)]
xy2 <- xy02_dt$xy


## countries
countries1_dt <- climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting",.(incidents=sum(incidents)),by=.(country)]

countries1_dt <- countries1_dt[order(-incidents)]
countries1_dt[,incidents_csum:=cumsum(incidents)]

num <- 3

excluded_countries1 <- countries1_dt[1:num]

countries1_dt <- countries1_dt[(num+1):nrow(countries1_dt)]

cty1 <- countries1_dt$country


countries2_dt <- climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting",.(incidents=sum(incidents)),by=.(country)]

countries2_dt <- countries2_dt[order(-incidents)]
countries2_dt[,incidents_csum:=cumsum(incidents)]

num <- 3

excluded_countries2 <- countries2_dt[1:num]

countries2_dt <- countries2_dt[(num+1):nrow(countries2_dt)]

cty2 <- countries2_dt$country


## B12: cells ----
reg3_pc <- feols(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting" & xy%in%xy1],vcov=conley(500)~x+y)

reg3_pn <- feols(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting" & xy%in%xy2],vcov=conley(500)~x+y)

modelsummary(list(reg3_pc,reg3_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))

impact3(climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting" & xy%in%xy1])

impact3(climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting" & xy%in%xy2])


## B13: countries ----
reg3_pc <- feols(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting" & country%in%cty1],vcov=conley(500)~x+y)

reg3_pn <- feols(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting" & country%in%cty2],vcov=conley(500)~x+y)

modelsummary(list(reg3_pc,reg3_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))

impact3(climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting" & country%in%cty1])

impact3(climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting" & country%in%cty2])


# R7 - croplands ----
cells_dt <- climatecropsconflict_dt[,.(area10=mean(area10)),by=.(xy)]
cells_dt <- cells_dt[order(-area10)]
one <- round(.01*nrow(cells_dt))
xy01_dt <- cells_dt[(one+1):nrow(cells_dt)]
xy1 <- xy01_dt$xy

countries_dt <- unique(climatecropsconflict_dt[yearmo=="2020-01" & area10>0,.(country,x,y,area10)])
countries_dt <- countries_dt[,.(cells=.N,cropland=sum(area10)),by=.(country)]
countries_dt[,cropland_cell:=cropland/cells]
countries_dt <- countries_dt[order(-cropland_cell)]
countries_dt <- countries_dt[cropland_cell>.5 & cropland_cell<2]
cty <- countries_dt$country

## B14: cells ----
reg3_pc <- feols(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting" & xy%in%xy1],vcov=conley(500)~x+y)

reg3_pn <- feols(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting" & xy%in%xy1],vcov=conley(500)~x+y)

modelsummary(list(reg3_pc,reg3_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))

impact3(climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting" & xy%in%xy1])

impact3(climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting" & xy%in%xy1])


## B15: countries ----
reg3_pc <- feols(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting" & country%in%cty],vcov=conley(500)~x+y)

reg3_pn <- feols(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting" & country%in%cty],vcov=conley(500)~x+y)

modelsummary(list(reg3_pc,reg3_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))

impact3(climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting" & country%in%cty])

impact3(climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting" & country%in%cty])


# R8 - postplanting ----

## B16: postplanting conflict ----
impact3pp <- function(x){
  r <- feols(incidents~area10:oni_croplag1:gs_tc:pp3+area10:oni_croplag1:gs_tc:I(1-pp3)+prec+tmax | xy + country^yearmo, data=x,vcov=conley(500)~x+y)
  m <- x[area10>0 & gs_tc>0,.(area10=mean(area10),gs_tc=mean(gs_tc),incidents=mean(incidents))]
  s <- 100*m$area10*m$gs_tc/m$incidents
  h1_coef <- round(r$coeftable["area10:oni_croplag1:gs_tc:pp3","Estimate"]*s,3)
  h1_se <- round(r$coeftable["area10:oni_croplag1:gs_tc:pp3","Std. Error"]*s,3)
  h1_stars <- pstars(r$coeftable["area10:oni_croplag1:gs_tc:pp3","Pr(>|t|)"])
  h1_est <- paste0(format(round(h1_coef,1),nsmall=1),h1_stars)
  h1_std <- paste0("(",format(round(h1_se,1),nsmall=1),")")
  return(list(descriptive=c(area_tc=round(m$area10,3),gs_tc=round(m$gs_tc,3),incidents_tc=round(m$incidents,3)),effect=c(h1_est,h1_std),output=c(h1_coef,h1_se)))
}

reg3_pc <- feols(incidents~area10:oni_croplag1:gs_tc+area10:oni_croplag1:gs_tc:pp3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg3_pn <- feols(incidents~area10:oni_croplag1:gs_tc+area10:oni_croplag1:gs_tc:pp3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)

modelsummary(list(reg3_pc,reg3_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))

impact3pp(climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"])

impact3pp(climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"])



