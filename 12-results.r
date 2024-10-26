library(data.table)
library(fixest)
library(zoo)
library(ggplot2)
library(ggforce)
library(ggnewscale)
library(cowplot)
library(Cairo)
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
climatecropsconflict_dt[,`:=`(area10=area/10000,ph3=ifelse(harv%in%1:3,1,0))]

## trend (for no particular reason)
climatecropsconflict_dt[,trend:=as.numeric(as.factor(yearmo))]

## total number of incidents by conflict category
climatecropsconflict_dt[,.(sum(incidents)),by=.(type,civilian_targeting)]


# Eq. 1 / Table 1 ----

impact1 <- function(x){
  r <- feols(incidents~area10:oni_crop:gs_tc+prec+tmax | xy + country^yearmo, data=x,vcov=conley(500)~x+y)
  m <- x[area10>0 & gs_tc>0,.(area10=mean(area10),gs_tc=mean(gs_tc),incidents=mean(incidents))]
  s <- 100*m$area10*m$gs_tc/m$incidents
  h1_coef <- round(r$coeftable["area10:oni_crop:gs_tc","Estimate"]*s,3)
  h1_se <- round(r$coeftable["area10:oni_crop:gs_tc","Std. Error"]*s,3)
  h1_stars <- pstars(r$coeftable["area10:oni_crop:gs_tc","Pr(>|t|)"])
  h1_est <- paste0(format(round(h1_coef,1),nsmall=1),h1_stars)
  h1_std <- paste0("(",format(round(h1_se,1),nsmall=1),")")
  return(list(descriptive=c(area_tc=round(m$area10,3),gs_tc=round(m$gs_tc,3),incidents_tc=round(m$incidents,3)),effect=c(h1_est,h1_std),output=c(h1_coef,h1_se)))
}

reg1_pc <- feols(incidents~area10:oni_crop:gs_tc+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg1_pn <- feols(incidents~area10:oni_crop:gs_tc+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)

modelsummary(list(reg1_pc,reg1_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))

impact1(climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"])

impact1(climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"])


# Eq. 2 / Table 2 ----

impact2 <- function(x){
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

reg2_pc <- feols(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg2_pn <- feols(incidents~area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)

modelsummary(list(reg2_pc,reg2_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))

impact2(climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"])

impact2(climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"])


# Eq. 3 / Table 3 ----

climatecropsconflict_dt[,neg:=ifelse(gs_tc_prec<0 & gs_tc_tmax>0,1,0)]
climatecropsconflict_dt[,pos:=1-neg]

climatecropsconflict_dt[is.na(neg)]$neg <- 0
climatecropsconflict_dt[is.na(pos)]$pos <- 0

impact2cells <- function(x){
  r <- feols(incidents~(area10:oni_crop_prec:gs_tc:ph3+area10:oni_crop_prec:gs_tc:I(1-ph3)):neg+(area10:oni_crop_prec:gs_tc:ph3+area10:oni_crop_prec:gs_tc:I(1-ph3)):pos+prec+tmax | xy + country^yearmo, data=x,vcov=conley(500)~x+y)
  m1 <- x[area10>0 & gs_tc>0 & neg==1,.(area10=mean(area10),gs_tc=mean(gs_tc),incidents=mean(incidents))]
  m2 <- x[area10>0 & gs_tc>0 & pos==1,.(area10=mean(area10),gs_tc=mean(gs_tc),incidents=mean(incidents))]
  s1 <- 100*m1$area10*m1$gs_tc/m1$incidents
  s2 <- 100*m2$area10*m2$gs_tc/m2$incidents
  h1_coef <- round(r$coeftable["area10:oni_crop_prec:gs_tc:ph3:neg","Estimate"]*s1,3)
  h1_se <- round(r$coeftable["area10:oni_crop_prec:gs_tc:ph3:neg","Std. Error"]*s1,3)
  h1_stars <- pstars(r$coeftable["area10:oni_crop_prec:gs_tc:ph3:neg","Pr(>|t|)"])
  h1_est <- paste0(format(round(h1_coef,1),nsmall=1),h1_stars)
  h1_std <- paste0("(",format(round(h1_se,1),nsmall=1),")")
  h2_coef <- round(r$coeftable["area10:oni_crop_prec:gs_tc:ph3:pos","Estimate"]*s2,3)
  h2_se <- round(r$coeftable["area10:oni_crop_prec:gs_tc:ph3:pos","Std. Error"]*s2,3)
  h2_stars <- pstars(r$coeftable["area10:oni_crop_prec:gs_tc:ph3:pos","Pr(>|t|)"])
  h2_est <- paste0(format(round(h2_coef,1),nsmall=1),h2_stars)
  h2_std <- paste0("(",format(round(h2_se,1),nsmall=1),")")
  return(list(descriptive=c(area_n=round(m1$area10,3),area_p=round(m2$area10,3),gs_tc_n=round(m1$gs_tc,3),gs_tc_p=round(m2$gs_tc,3),incidents_n=round(m1$incidents,3),incidents_p=round(m2$incidents,3)),effect=c(h1_est,h1_std,h2_est,h2_std),output=c(h1_coef,h1_se,h2_coef,h2_se)))
}

reg2_pc <- feols(incidents~(area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3):neg+(area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3):pos+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg2_pn <- feols(incidents~(area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3):neg+(area10:oni_crop:gs_tc+area10:oni_crop:gs_tc:ph3):pos+prec+tmax | xy + country^yearmo, data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)

modelsummary(list(reg2_pc,reg2_pn),estimate="{estimate}{stars}",stars=c('*'=.1,'**'=.05,'***'=.01),fmt=fmt_sprintf("%.4f"))

impact2cells(climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"])

impact2cells(climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"])


# Eq. 4 / Figure 5 ----

reg5_pc <- feols(incidents~area10:oni_cropevent:gs_tc:i(harv,keep=c(1:12))+prec+tmax | xy+country^yearmo,data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting"],vcov=conley(500)~x+y)

reg5_pn <- feols(incidents~area10:oni_cropevent:gs_tc:i(harv,keep=c(1:12))+prec+tmax | xy+country^yearmo,data=climatecropsconflict_dt[type=="Political violence" & civilian_targeting!="Civilian_targeting"],vcov=conley(500)~x+y)

m_pc <- climatecropsconflict_dt[area10>0 & type=="Political violence" & civilian_targeting=="Civilian_targeting",.(area10=mean(area10),gs_tc=mean(gs_tc),incidents=mean(incidents))]
s_pc <- 100*m_pc$area10*m_pc$gs_tc/m_pc$incidents
r_pc <- reg5_pc$coeftable[-c(1:2),][c(9:12,1:8),1:2]
h_pc <- data.table(r_pc*s_pc)
colnames(h_pc) <- c("Effect","SE")
h_pc[,Conflict:="One-sided"]
h_pc[,Period:=c(-4:7)]

m_pn <- climatecropsconflict_dt[area10>0 & type=="Political violence" & civilian_targeting=="Civilian_targeting",.(area10=mean(area10),gs_tc=mean(gs_tc),incidents=mean(incidents))]
s_pn <- 100*m_pn$area10*m_pn$gs_tc/m_pn$incidents
r_pn <- reg5_pn$coeftable[-c(1:2),][c(9:12,1:8),1:2]
h_pn <- data.table(r_pn*s_pn)
colnames(h_pn) <- c("Effect","SE")
h_pn[,Conflict:="Two-sided"]
h_pn[,Period:=c(-4:7)]

effect_dt <- data.table(Reduce(rbind,list(h_pc,h_pn)))
effect_dt[,StatSig:=ifelse(abs(Effect/SE)>1.96,1,0)]

effect_dt$Period <- factor(effect_dt$Period,levels=c(-4:7))

gg <- ggplot(effect_dt,aes(x=Period,y=Effect))+
  geom_hline(yintercept = 0,linewidth=.4,color="dimgray")+
  geom_pointrange(aes(ymin=Effect-1.96*SE,ymax=Effect+1.96*SE),position=position_dodge(width=.4),color="black",fill="black",shape=21,size=.5,stroke=.2,linewidth=.6)+
  annotate("rect",xmin=4.5,xmax=7.5,ymin=-6.5,ymax=3.5,color="dimgray",fill=NA,linetype=2,linewidth=.4)+
  facet_wrap(.~Conflict,ncol=2,scales="free",labeller = label_wrap_gen(multi_line=F))+
  labs(x="Month from harvest",y="% chg (relative to baseline conflict)")+
  theme_minimal_hgrid()+
  theme(plot.title.position="plot",panel.background = element_rect(color=NA,fill="white"),plot.background = element_rect(color=NA,fill="white"),axis.text=element_text(size=10,color="dimgray"),panel.grid.major.y = element_line(linewidth=.4,linetype=3,color="dimgray"),panel.spacing=unit(1,"lines"),axis.ticks = element_blank(),axis.line.x=element_blank(),axis.title=element_text(size=12),plot.subtitle=element_text(size=12))

ggsave("figures/fig5_event.png",gg,width=6.5,height=6.5*9/16,dpi="retina")
ggsave("figures/fig5_event.eps",gg,width=6.5,height=6.5*9/16,dpi="retina")

