library(data.table)
library(sf)
library(stringr)
library(zoo)
library(lmtest)
library(sandwich)
library(rnaturalearth)
library(rnaturalearthdata)

rm(list=ls())
gc()

world <- ne_countries(returnclass="sf",scale="large")
africa <- ne_countries(returnclass="sf",scale="large",continent="Africa")
europe <- ne_countries(returnclass="sf",scale="large",continent="Europe")
asia <- ne_countries(returnclass="sf",scale="large",continent="Asia")

crs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"


sf_use_s2(FALSE)

world <- st_set_crs(world,crs)
africa <- st_set_crs(africa,crs)
europe <- st_set_crs(europe,crs)
asia <- st_set_crs(asia,crs)

lakes <- ne_download(returnclass="sf",scale="large",type="lakes",category="physical")
lakes <- st_set_crs(lakes,crs)
lakes <- st_make_valid(lakes)
lakes <- st_intersection(africa,lakes)

rivers <- ne_download(returnclass="sf",scale="large",type="rivers_lake_centerlines",category="physical")
rivers <- st_set_crs(rivers,crs)
rivers <- st_make_valid(rivers)
rivers <- st_intersection(africa,rivers)

ocean <- ne_download(returnclass="sf",scale="large",type="ocean",category="physical")
ocean <- st_set_crs(ocean,crs)

sf_use_s2(TRUE)


enso_yrs <- 1979:2023

"%!in%" <- Negate("%in%")

reg_out <- function(y,x){
  tr <- 1:length(y)
  r <- lm(y~x+tr)
  b <- coefficients(summary(r))["x","Estimate"]
  s <- coefficients(summary(r))["x","t value"]
  p <- coefficients(summary(r))["x","Pr(>|t|)"]
  return(list(b,s,p))
}

# load crops
load("data/calendar.RData")
load("data/spam.RData")

# load climate
load("data/precipitation.RData")
load("data/temperature.RData")


# fetch oni (enso) index 
oni_dt <- fread("https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt")

oni_dt <- oni_dt[SEAS=="NDJ",.(elnino_year=as.numeric(YR),oni=ANOM)]


# sort out the temperature dataset
tmax_dt[,xy:=paste(x,y,sep=",")]

dw <- dcast(unique(tmax_dt[,.(xy,year)]),xy~year,value.var="xy")#[,.N,by=.(year)]
dw <- dw[complete.cases(dw)]
dl <- melt(dw,id.vars="xy")

sub_dt <- dl[,.(xy,year=as.numeric(as.character(variable)))]
tmax_dt[,year:=as.numeric(as.character(year))]

tmax_dt <- merge(tmax_dt,sub_dt,by=c("year","xy"))
unique(tmax_dt[,.(xy,year)])[,.N,by=.(year)]

tmax_dt$xy <- NULL


# sort out the precipitation dataset
prec_dt[,xy:=paste(x,y,sep=",")]

dw <- dcast(unique(prec_dt[,.(xy,year)]),xy~year,value.var="xy")#[,.N,by=.(year)]
dw <- dw[complete.cases(dw)]
dl <- melt(dw,id.vars="xy")

sub_dt <- dl[,.(xy,year=as.numeric(as.character(variable)))]
prec_dt[,year:=as.numeric(as.character(year))]

prec_dt <- merge(prec_dt,sub_dt,by=c("year","xy"))
unique(prec_dt[,.(xy,year)])[,.N,by=.(year)]

prec_dt$xy <- NULL

# put the weather variables together
weather_dt <- merge(prec_dt,tmax_dt,by=c("year","mo","x","y"))
weather_dt <- weather_dt[order(x,y,year,mo)]

weather_dt[,`:=`(prec3=frollmean(prec,3,align="c"),tmax3=frollmean(tmax,3,align="c")),by=.(x,y)]

weather_dt[is.na(prec3)]$prec3 <- weather_dt[is.na(prec3)]$prec
weather_dt[is.na(tmax3)]$tmax3 <- weather_dt[is.na(tmax3)]$tmax

weather_dt[,elnino_year:=year]
weather_dt[mo%in%unique(weather_dt$mo)[1:5]]$elnino_year <- NA
weather_dt[,elnino_year:=c(rep(NA,5),na.locf(elnino_year)),by=.(x,y)]
weather_dt[year==min(enso_yrs) & mo%in%unique(weather_dt$mo)[1:5]]$elnino_year <- min(enso_yrs)-1
weather_dt <- weather_dt[order(x,y,year,mo)]


# bring in ENSO
climate_dt <- merge(weather_dt,oni_dt,by="elnino_year")

# combine croplands and calendars
crops_dt <- merge(gs_mid_dt,spam_comb_dt,by=c("x","y"))

crops_dt[,crop:=substring(colnames(crops_dt[,.(area_maize,area_sorghum,area_millet,area_rice,area_wheat)])[max.col(crops_dt[,.(area_maize,area_sorghum,area_millet,area_rice,area_wheat)],ties.method="first")],first=6)]

crops_dt[,area:=apply(crops_dt[,.(area_maize,area_sorghum,area_millet,area_rice,area_wheat)],1,max)]

crops_dt[,crop:=ifelse(area>0,crop,"none")]

crops_dt[,harv:=ifelse(area>0 & area==area_maize,Maize_harvest_mid,ifelse(area>0 & area==area_sorghum,Sorghum_harvest_mid,ifelse(area>0 & area==area_wheat,Wheat.Winter_harvest_mid,ifelse(area>0 & area==area_rice,Rice_harvest_mid,ifelse(area>0,Millet_harvest_mid,0)))))]

crops_dt[,gs:=ifelse(area>0 & area==area_maize,Maize_gs_mid,ifelse(area>0 & area==area_sorghum,Sorghum_gs_mid,ifelse(area>0 & area==area_wheat,Wheat.Winter_gs_mid,ifelse(area>0 & area==area_rice,Rice_gs_mid,ifelse(area>0,Millet_gs_mid,0)))))]

crops_dt[,gs:=ifelse(gs==0,0,ifelse(gs==1 & harv==1,0,1))]

crops_dt[,plant:=ifelse(area>0 & area==area_maize,Maize_plant_mid,ifelse(area>0 & area==area_sorghum,Sorghum_plant_mid,ifelse(area>0 & area==area_wheat,Wheat.Winter_plant_mid,ifelse(area>0 & area==area_rice,Rice_plant_mid,ifelse(area>0,Millet_plant_mid,0)))))]

crops_dt <- crops_dt[,.(x,y,mo,month=Month,crop,area,harv,gs,plant,area_maize,area_sorghum,area_millet,area_rice,area_wheat,area_sum=area_maize+area_sorghum+area_millet+area_rice+area_wheat)]

crops_dt <- crops_dt[complete.cases(crops_dt)]

crops_dt[,mo:=sprintf("%02d",mo)]

combined_dt <- merge(climate_dt,crops_dt,by=c("x","y","mo"),all.x=T)

# add countries
sub_dt <- unique(combined_dt[,.(x,y)])
points_sf <- st_as_sf(sub_dt,coords=c("x","y"),crs=st_crs(crs))
countries_sf <- st_make_valid(africa)
indices <- st_join(points_sf,countries_sf)
sub_dt$country <- indices$name
sub_dt <- sub_dt[complete.cases(sub_dt)]
combined_dt <- merge(combined_dt,sub_dt,by=c("x","y"),all.x=T)
combined_dt <- combined_dt[!is.na(country)]
combined_dt[is.na(combined_dt)] <- 0
combined_dt[crop==0]$crop <- "none"
combined_dt[crop=="0"]$crop <- "none"
combined_dt[,month:=month.abb[as.numeric(mo)]]

combined_dt <- combined_dt[order(x,y,year,mo)]

combined_dt <- combined_dt[elnino_year%in%enso_yrs]


combined_dt[,gs_temp:=ifelse(gs==1,1,NA)]

combined_dt <- combined_dt[,`:=`(yr_prec=mean(prec3),yr_tmax=mean(tmax3)),by=.(x,y)]

combined_dt <- combined_dt[,`:=`(gs_prec=mean(prec3*gs_temp,na.rm=T),gs_tmax=mean(tmax3*gs_temp,na.rm=T)),by=.(x,y)]

combined_dt <- combined_dt[,`:=`(yr_prec_year=mean(prec3),yr_tmax_year=mean(tmax3)),by=.(x,y,elnino_year)]


sub_dt <- unique(combined_dt[,.(x,y,year,mo,prec3,tmax3,gs)])
seasons_dt <- sub_dt[,.(prec_mo=mean(prec3),tmax_mo=mean(tmax3)),by=.(x,y,mo)]

combined_dt <- merge(combined_dt,seasons_dt,by=c("x","y","mo"))

combined_dt <- combined_dt[order(country,x,y,year,mo)]


## monthly teleconnections
combined_dt[prec_mo>=1 | gs==1,c("mo_prec_b","mo_prec_s","mo_prec_p"):=reg_out(y=prec3,x=oni),by=.(x,y,mo)]
combined_dt[,c("mo_tmax_b","mo_tmax_s","mo_tmax_p"):=reg_out(y=tmax3,x=oni),by=.(x,y,mo)]

combined_dt[,`:=`(mo_tc=ifelse(mo_prec_p<.05 | mo_tmax_p<.05,1,0),mo_tcp=ifelse(mo_prec_p<.05,1,0),mo_tct=ifelse(mo_tmax_p<.05,1,0))]

combined_dt[,`:=`(tc=sum(mo_tc,na.rm=T)/12,tcp=sum(mo_tcp,na.rm=T)/12,tct=sum(mo_tct,na.rm=T)/12,gs_tc=ifelse(sum(gs)>0,sum(mo_tc*gs_temp,na.rm=T)/sum(gs_temp,na.rm=T),0),gs_tcp=ifelse(sum(gs)>0,sum(mo_tcp*gs_temp,na.rm=T)/sum(gs_temp,na.rm=T),0),gs_tct=ifelse(sum(gs)>0,sum(mo_tct*gs_temp,na.rm=T)/sum(gs_temp,na.rm=T),0),tc_prec=mean(mo_prec_b,na.rm=T),tc_tmax=mean(mo_tmax_b,na.rm=T),gs_tc_prec=mean(mo_prec_b*gs_temp,na.rm=T),gs_tc_tmax=mean(mo_tmax_b*gs_temp,na.rm=T)),by=.(x,y,elnino_year)]


climatecrops_dt <- combined_dt#merge(combined_dt,subset_dt,by=c("x","y"))

climatecrops_dt <- climatecrops_dt[order(country,x,y,year,mo)]

climatecrops_dt[,crop_year:=year]
climatecrops_dt[harv!=1,crop_year:=NA]
climatecrops_dt[year==min(year) & mo=="06" & harv!=1,crop_year:=min(year)-1]

climatecrops_dt[,crop_year:=na.locf(crop_year),by=.(x,y)]
climatecrops_dt[crop=="none",crop_year:=year]


## the effective ONI (for crops with gs across two enso years)
teleconnected_dt <- climatecrops_dt[crop!="none" & gs_tc>0,.(x,y,gs_temp,mo_tc,oni,crop_year)]
teleconnected_dt[,mo_tc:=ifelse(mo_tc==1,1,NA)]
teleconnected_dt <- teleconnected_dt[complete.cases(teleconnected_dt)]
teleconnected_dt$mo_tc <- NULL

teleconnected_dt <- teleconnected_dt[,.(oni_crop=as.numeric(names(which.max(table(factor(oni,levels=rev(unique(oni)))))))),by=.(x,y,crop_year)]

teleconnected_dt[,oni_crop:=shift(oni_crop,1),by=.(x,y)]


nonteleconnected_dt <- climatecrops_dt[crop!="none" & gs_tc==0,.(x,y,gs_temp,oni,crop_year)]
nonteleconnected_dt <- nonteleconnected_dt[complete.cases(nonteleconnected_dt)]

nonteleconnected_dt <- nonteleconnected_dt[,.(oni_crop=as.numeric(names(which.max(table(factor(oni,levels=rev(unique(oni)))))))),by=.(x,y,crop_year)]

nonteleconnected_dt[,oni_crop:=shift(oni_crop,1),by=.(x,y)]


comb_dt <- rbind(teleconnected_dt,nonteleconnected_dt)
colnames(comb_dt) <- c("x","y","crop_year","oni_crop")

teleconnected_dt <- climatecrops_dt[crop!="none" & gs_tcp>0,.(x,y,gs_temp,mo_tcp,oni,crop_year)]
teleconnected_dt[,mo_tcp:=ifelse(mo_tcp==1,1,NA)]
teleconnected_dt <- teleconnected_dt[complete.cases(teleconnected_dt)]
teleconnected_dt$mo_tcp <- NULL

teleconnected_dt <- teleconnected_dt[,.(oni_crop=as.numeric(names(which.max(table(factor(oni,levels=rev(unique(oni)))))))),by=.(x,y,crop_year)]

teleconnected_dt[,oni_crop:=shift(oni_crop,1),by=.(x,y)]


nonteleconnected_dt <- climatecrops_dt[crop!="none" & gs_tcp==0,.(x,y,gs_temp,oni,crop_year)]
nonteleconnected_dt <- nonteleconnected_dt[complete.cases(nonteleconnected_dt)]

nonteleconnected_dt <- nonteleconnected_dt[,.(oni_crop=as.numeric(names(which.max(table(factor(oni,levels=rev(unique(oni)))))))),by=.(x,y,crop_year)]

nonteleconnected_dt[,oni_crop:=shift(oni_crop,1),by=.(x,y)]


combp_dt <- rbind(teleconnected_dt,nonteleconnected_dt)
colnames(combp_dt) <- c("x","y","crop_year","oni_crop_prec")


teleconnected_dt <- climatecrops_dt[crop!="none" & gs_tct>0,.(x,y,gs_temp,mo_tct,oni,crop_year)]
teleconnected_dt[,mo_tct:=ifelse(mo_tct==1,1,NA)]
teleconnected_dt <- teleconnected_dt[complete.cases(teleconnected_dt)]
teleconnected_dt$mo_tct <- NULL

teleconnected_dt <- teleconnected_dt[,.(oni_crop=as.numeric(names(which.max(table(factor(oni,levels=rev(unique(oni)))))))),by=.(x,y,crop_year)]

teleconnected_dt[,oni_crop:=shift(oni_crop,1),by=.(x,y)]


nonteleconnected_dt <- climatecrops_dt[crop!="none" & gs_tct==0,.(x,y,gs_temp,oni,crop_year)]
nonteleconnected_dt <- nonteleconnected_dt[complete.cases(nonteleconnected_dt)]

nonteleconnected_dt <- nonteleconnected_dt[,.(oni_crop=as.numeric(names(which.max(table(factor(oni,levels=rev(unique(oni)))))))),by=.(x,y,crop_year)]

nonteleconnected_dt[,oni_crop:=shift(oni_crop,1),by=.(x,y)]

combt_dt <- rbind(teleconnected_dt,nonteleconnected_dt)
colnames(combt_dt) <- c("x","y","crop_year","oni_crop_tmax")

climatecrops_dt <- merge(climatecrops_dt,comb_dt,by=c("x","y","crop_year"),all.x=T)
climatecrops_dt <- merge(climatecrops_dt,combp_dt,by=c("x","y","crop_year"),all.x=T)
climatecrops_dt <- merge(climatecrops_dt,combt_dt,by=c("x","y","crop_year"),all.x=T)

climatecrops_dt[crop=="none",`:=`(oni_crop=oni)]
climatecrops_dt[crop=="none",`:=`(oni_crop_prec=oni)]
climatecrops_dt[crop=="none",`:=`(oni_crop_tmax=oni)]

climatecrops_dt[,date:=as.Date(paste0(year,"-",mo,"-01"))]

climatecrops_dt[date>="1981-06-01",`:=`(oni_crop=na.locf(oni_crop,fromLast=F),oni_crop_prec=na.locf(oni_crop_prec,fromLast=F),oni_crop_tmax=na.locf(oni_crop_tmax,fromLast=F)),by=.(x,y)]

climatecrops_dt[,`:=`(oni_crop=na.locf(oni_crop,fromLast=T),oni_crop_prec=na.locf(oni_crop_prec,fromLast=T),oni_crop_tmax=na.locf(oni_crop_tmax,fromLast=T)),by=.(x,y)]

climatecrops_dt[,`:=`(oni_croplag1=shift(oni_crop,12,type="lag"),oni_croplag2=shift(oni_crop,24,type="lag"),oni_croplead1=shift(oni_crop,12,type="lead"),oni_croplead2=shift(oni_crop,24,type="lead")),by=.(x,y)]


## 4-month lag (to enable event study like analysis)
climatecrops_dt[,`:=`(oni_cropevent=shift(oni_crop,4,type="lead")),by=.(x,y)]
climatecrops_dt[,`:=`(oni_cropevent=na.locf(oni_cropevent,fromLast=F)),by=.(x,y)]


## growing season weather
climatecrops_dt <- climatecrops_dt[,`:=`(gs_prec_year1=mean(prec3*gs_temp,na.rm=T),gs_tmax_year1=mean(tmax3*gs_temp,na.rm=T)),by=.(x,y,crop_year)]

climatecrops_dt[is.na(gs_prec_year1)]$gs_prec_year1 <- 0
climatecrops_dt[is.na(gs_tmax_year1)]$gs_tmax_year1 <- 0

climatecrops_dt[,`:=`(gs_prec_year=shift(gs_prec_year1,12),gs_tmax_year=shift(gs_tmax_year1,12)),by=.(x,y)]

climatecrops_dt <- climatecrops_dt[order(country,x,y,year,mo)]

climatecrops_dt[,`:=`(gs_prec_year=na.locf(gs_prec_year,fromLast=T),gs_tmax_year=na.locf(gs_tmax_year,fromLast=T)),by=.(x,y)]
climatecrops_dt[,`:=`(gs_prec_year=na.locf(gs_prec_year,fromLast=F),gs_tmax_year=na.locf(gs_tmax_year,fromLast=F)),by=.(x,y)]


climatecrops_dt$gs_temp <- NULL
climatecrops_dt$gs_last <- NULL


save(climatecrops_dt,file="data/climatecrops.RData")

