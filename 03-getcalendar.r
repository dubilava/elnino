# load the packages
library(data.table)
library(ggplot2)
library(raster)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ncdf4)
library(R.utils)

rm(list=ls())
gc()

# you will need to have downloaded the relevant netCDF files from 
# https://sage.nelson.wisc.edu/data-and-models/datasets/crop-calendar-dataset/netcdf-0-5-degree/
# and stored it in the same folder as the present R file;
# then proceed as follows.

# load the Africa map
africa <- ne_countries(continent="Africa",returnclass="sf",scale="large")
africa <- st_set_crs(africa,"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

crop_list <- c("Maize","Sorghum","Millet","Rice","Wheat.Winter")

plant_srt_list <- list()
plant_mid_list <- list()
plant_end_list <- list()

harvest_srt_list <- list()
harvest_mid_list <- list()
harvest_end_list <- list()

area_list <- list()

for(i in 1:length(crop_list)){
  
  crop <- crop_list[i]
  
  gunzip(paste(crop,".crop.calendar.fill.nc.gz",sep=""),remove=F,skip=T)
  gunzip(paste(crop,".crop.calendar.nc.gz",sep=""),remove=F,skip=T)
  
  plant_srt <- raster(paste(crop,".crop.calendar.nc",sep=""),varname="plant.start")
  plant_mid <- raster(paste(crop,".crop.calendar.nc",sep=""),varname="plant")
  plant_end <- raster(paste(crop,".crop.calendar.nc",sep=""),varname="plant.end")
  
  harvest_srt <- raster(paste(crop,".crop.calendar.nc",sep=""),varname="harvest.start")
  harvest_mid <- raster(paste(crop,".crop.calendar.nc",sep=""),varname="harvest")
  harvest_end <- raster(paste(crop,".crop.calendar.nc",sep=""),varname="harvest.end")
  
  
  raster_mask <- mask(plant_srt,africa)
  rm <- rasterToPoints(raster_mask)
  rm[,1:2] <- round(rm[,1:2],2)
  dm_ps <- data.table(rm)
  colnames(dm_ps)[3] <- "plant_srt"
  dm_ps$plant_srt_date <- as.Date(dm_ps$plant_srt,origin="2000-01-01")
  dm_ps$plant_srt <- as.numeric(substr(dm_ps$plant_srt_date,6,7))
  dm_ps <- dm_ps[order(x,y)]
  
  raster_mask <- mask(plant_mid,africa)
  rm <- rasterToPoints(raster_mask)
  rm[,1:2] <- round(rm[,1:2],2)
  dm_pm <- data.table(rm)
  colnames(dm_pm)[3] <- "plant_mid"
  dm_pm$plant_date <- as.Date(dm_pm$plant_mid,origin="2000-01-01")
  dm_pm$plant_mid <- as.numeric(substr(dm_pm$plant_date,6,7))
  dm_pm <- dm_pm[order(x,y)]
  
  raster_mask <- mask(plant_end,africa)
  rm <- rasterToPoints(raster_mask)
  rm[,1:2] <- round(rm[,1:2],2)
  dm_pe <- data.table(rm)
  colnames(dm_pe)[3] <- "plant_end"
  dm_pe$plant_end_date <- as.Date(dm_pe$plant_end,origin="2000-01-01")
  dm_pe$plant_end <- as.numeric(substr(dm_pe$plant_end_date,6,7))
  dm_pe <- dm_pe[order(x,y)]
  
  
  raster_mask <- mask(harvest_srt,africa)
  rm <- rasterToPoints(raster_mask)
  rm[,1:2] <- round(rm[,1:2],2)
  dm_hs <- data.table(rm)
  colnames(dm_hs)[3] <- "harvest_srt"
  dm_hs$harvest_srt_date <- as.Date(dm_hs$harvest_srt,origin="2000-01-01")
  dm_hs$harvest_srt <- as.numeric(substr(dm_hs$harvest_srt_date,6,7))
  dm_hs <- dm_hs[order(x,y)]
  
  raster_mask <- mask(harvest_mid,africa)
  rm <- rasterToPoints(raster_mask)
  rm[,1:2] <- round(rm[,1:2],2)
  dm_hm <- data.table(rm)
  colnames(dm_hm)[3] <- "harvest_mid"
  dm_hm$harvest_date <- as.Date(dm_hm$harvest_mid,origin="2000-01-01")
  dm_hm$harvest_mid <- as.numeric(substr(dm_hm$harvest_date,6,7))
  dm_hm <- dm_hm[order(x,y)]
  
  
  raster_mask <- mask(harvest_end,africa)
  rm <- rasterToPoints(raster_mask)
  rm[,1:2] <- round(rm[,1:2],2)
  dm_he <- data.table(rm)
  colnames(dm_he)[3] <- "harvest_end"
  dm_he$harvest_end_date <- as.Date(dm_he$harvest_end,origin="2000-01-01")
  dm_he$harvest_end <- as.numeric(substr(dm_he$harvest_end_date,6,7))
  dm_he <- dm_he[order(x,y)]
  
  crop_dt <- Reduce(function(...) merge(...,by=c("x","y"),all=T),list(dm_ps,dm_pm,dm_pe,dm_hs,dm_hm,dm_he))
  crop_dt[,`:=`(plant_length=ifelse(plant_end_date-plant_srt_date+1>0,plant_end_date-plant_srt_date+1,plant_end-plant_srt+1+365),harvest_length=ifelse(harvest_end-harvest_srt+1>0,harvest_end_date-harvest_srt_date+1,harvest_end_date-harvest_srt_date+1+365))]
  
  
  plant_dt <- crop_dt[,.(x,y,plant_srt,plant_mid,plant_end)]
  colnames(plant_dt)[3:5] <- paste0(crop,c("_plant_srt","_plant_mid","_plant_end"),sep="")
  
  ## plant calendars
  plant_srt_dt <- matrix(0,nrow=nrow(plant_dt),ncol=12)
  colnames(plant_srt_dt) <- as.character(1:12)
  for(j in 1:nrow(plant_srt_dt)){
    if(!is.na(plant_dt[j,3])){
      plant_srt_dt[j,((as.numeric(plant_dt[j,3])):12)] <- 1:(12-as.numeric(plant_dt[j,3])+1)
      if(as.numeric(plant_dt[j,3])!=1){
        plant_srt_dt[j,(1:(as.numeric(plant_dt[j,3])-1))] <- (12-as.numeric(plant_dt[j,3])+2):12
      }
    }
  }
  plant_srt_dt <- as.data.table(plant_srt_dt)
  
  plant_mid_dt <- matrix(0,nrow=nrow(plant_dt),ncol=12)
  colnames(plant_mid_dt) <- as.character(1:12)
  for(j in 1:nrow(plant_mid_dt)){
    if(!is.na(plant_dt[j,4])){
      plant_mid_dt[j,((as.numeric(plant_dt[j,4])):12)] <- 1:(12-as.numeric(plant_dt[j,4])+1)
      if(as.numeric(plant_dt[j,4])!=1){
        plant_mid_dt[j,(1:(as.numeric(plant_dt[j,4])-1))] <- (12-as.numeric(plant_dt[j,4])+2):12
      }
    }
  }
  plant_mid_dt <- as.data.table(plant_mid_dt)
  
  plant_end_dt <- matrix(0,nrow=nrow(plant_dt),ncol=12)
  colnames(plant_end_dt) <- as.character(1:12)
  for(j in 1:nrow(plant_end_dt)){
    if(!is.na(plant_dt[j,5])){
      plant_end_dt[j,((as.numeric(plant_dt[j,5])):12)] <- 1:(12-as.numeric(plant_dt[j,5])+1)
      if(as.numeric(plant_dt[j,5])!=1){
        plant_end_dt[j,(1:(as.numeric(plant_dt[j,5])-1))] <- (12-as.numeric(plant_dt[j,5])+2):12
      }
    }
  }
  plant_end_dt <- as.data.table(plant_end_dt)
  
  
  harvest_dt <- crop_dt[,.(x,y,harvest_srt,harvest_mid,harvest_end)]
  colnames(harvest_dt)[3:5] <- paste0(crop,c("_harvest_srt","_harvest_mid","_harvest_end"),sep="")
  
  ## harvest calendars
  harvest_srt_dt <- matrix(0,nrow=nrow(harvest_dt),ncol=12)
  colnames(harvest_srt_dt) <- as.character(1:12)
  for(j in 1:nrow(harvest_srt_dt)){
    if(!is.na(harvest_dt[j,3])){
      harvest_srt_dt[j,((as.numeric(harvest_dt[j,3])):12)] <- 1:(12-as.numeric(harvest_dt[j,3])+1)
      if(as.numeric(harvest_dt[j,3])!=1){
        harvest_srt_dt[j,(1:(as.numeric(harvest_dt[j,3])-1))] <- (12-as.numeric(harvest_dt[j,3])+2):12
      }
    }
  }
  harvest_srt_dt <- as.data.table(harvest_srt_dt)
  
  harvest_mid_dt <- matrix(0,nrow=nrow(harvest_dt),ncol=12)
  colnames(harvest_mid_dt) <- as.character(1:12)
  for(j in 1:nrow(harvest_mid_dt)){
    if(!is.na(harvest_dt[j,4])){
      harvest_mid_dt[j,((as.numeric(harvest_dt[j,4])):12)] <- 1:(12-as.numeric(harvest_dt[j,4])+1)
      if(as.numeric(harvest_dt[j,4])!=1){
        harvest_mid_dt[j,(1:(as.numeric(harvest_dt[j,4])-1))] <- (12-as.numeric(harvest_dt[j,4])+2):12
      }
    }
  }
  harvest_mid_dt <- as.data.table(harvest_mid_dt)
  
  harvest_end_dt <- matrix(0,nrow=nrow(harvest_dt),ncol=12)
  colnames(harvest_end_dt) <- as.character(1:12)
  for(j in 1:nrow(harvest_end_dt)){
    if(!is.na(harvest_dt[j,5])){
      harvest_end_dt[j,((as.numeric(harvest_dt[j,5])):12)] <- 1:(12-as.numeric(harvest_dt[j,5])+1)
      if(as.numeric(harvest_dt[j,5])!=1){
        harvest_end_dt[j,(1:(as.numeric(harvest_dt[j,5])-1))] <- (12-as.numeric(harvest_dt[j,5])+2):12
      }
    }
  }
  harvest_end_dt <- as.data.table(harvest_end_dt)
  
  
  plant_srt_comb_dt <- data.table(plant_dt[,1:2],plant_srt_dt)
  plant_srt_long_dt <- data.table::melt(plant_srt_comb_dt,id.vars=c("x","y"),variable.name="mo",value.name=paste(crop,"_plant_srt",sep=""))
  plant_srt_long_dt <- plant_srt_long_dt[order(x,y,mo)]
  
  plant_mid_comb_dt <- data.table(plant_dt[,1:2],plant_mid_dt)
  plant_mid_long_dt <- data.table::melt(plant_mid_comb_dt,id.vars=c("x","y"),variable.name="mo",value.name=paste(crop,"_plant_mid",sep=""))
  plant_mid_long_dt <- plant_mid_long_dt[order(x,y,mo)]
  
  plant_end_comb_dt <- data.table(plant_dt[,1:2],plant_end_dt)
  plant_end_long_dt <- data.table::melt(plant_end_comb_dt,id.vars=c("x","y"),variable.name="mo",value.name=paste(crop,"_plant_end",sep=""))
  plant_end_long_dt <- plant_end_long_dt[order(x,y,mo)]
  
  harvest_srt_comb_dt <- data.table(harvest_dt[,1:2],harvest_srt_dt)
  harvest_srt_long_dt <- data.table::melt(harvest_srt_comb_dt,id.vars=c("x","y"),variable.name="mo",value.name=paste(crop,"_harvest_srt",sep=""))
  harvest_srt_long_dt <- harvest_srt_long_dt[order(x,y,mo)]
  
  harvest_mid_comb_dt <- data.table(harvest_dt[,1:2],harvest_mid_dt)
  harvest_mid_long_dt <- data.table::melt(harvest_mid_comb_dt,id.vars=c("x","y"),variable.name="mo",value.name=paste(crop,"_harvest_mid",sep=""))
  harvest_mid_long_dt <- harvest_mid_long_dt[order(x,y,mo)]
  
  harvest_end_comb_dt <- data.table(harvest_dt[,1:2],harvest_end_dt)
  harvest_end_long_dt <- data.table::melt(harvest_end_comb_dt,id.vars=c("x","y"),variable.name="mo",value.name=paste(crop,"_harvest_end",sep=""))
  harvest_end_long_dt <- harvest_end_long_dt[order(x,y,mo)]
  
  
  plant_srt_list[[i]] <- plant_srt_long_dt
  plant_mid_list[[i]] <- plant_mid_long_dt
  plant_end_list[[i]] <- plant_end_long_dt
  harvest_srt_list[[i]] <- harvest_srt_long_dt
  harvest_mid_list[[i]] <- harvest_mid_long_dt
  harvest_end_list[[i]] <- harvest_end_long_dt
  
}



plant_srt_dt <- Reduce(function(...) merge(...,by=c("x","y","mo"),all=T),plant_srt_list)
plant_mid_dt <- Reduce(function(...) merge(...,by=c("x","y","mo"),all=T),plant_mid_list)
plant_end_dt <- Reduce(function(...) merge(...,by=c("x","y","mo"),all=T),plant_end_list)
harvest_srt_dt <- Reduce(function(...) merge(...,by=c("x","y","mo"),all=T),harvest_srt_list)
harvest_mid_dt <- Reduce(function(...) merge(...,by=c("x","y","mo"),all=T),harvest_mid_list)
harvest_end_dt <- Reduce(function(...) merge(...,by=c("x","y","mo"),all=T),harvest_end_list)


plant_srt_dt$Month <- month.abb[plant_srt_dt$mo]
plant_srt_dt$Month <- factor(plant_srt_dt$Month,levels=unique(plant_srt_dt$Month))

plant_mid_dt$Month <- month.abb[plant_mid_dt$mo]
plant_mid_dt$Month <- factor(plant_mid_dt$Month,levels=unique(plant_mid_dt$Month))

plant_end_dt$Month <- month.abb[plant_end_dt$mo]
plant_end_dt$Month <- factor(plant_end_dt$Month,levels=unique(plant_end_dt$Month))


harvest_srt_dt$Month <- month.abb[harvest_srt_dt$mo]
harvest_srt_dt$Month <- factor(harvest_srt_dt$Month,levels=unique(harvest_srt_dt$Month))

harvest_mid_dt$Month <- month.abb[harvest_mid_dt$mo]
harvest_mid_dt$Month <- factor(harvest_mid_dt$Month,levels=unique(harvest_mid_dt$Month))

harvest_end_dt$Month <- month.abb[harvest_end_dt$mo]
harvest_end_dt$Month <- factor(harvest_end_dt$Month,levels=unique(harvest_end_dt$Month))


gs_mid_dt <- merge(plant_mid_dt[],harvest_mid_dt,by=c("x","y","mo","Month"))
gs_mid_dt[,`:=`(Maize_gs_mid=ifelse(Maize_plant_mid-Maize_harvest_mid<0 | Maize_harvest_mid==1,1,0),Sorghum_gs_mid=ifelse(Sorghum_plant_mid-Sorghum_harvest_mid<0 | Sorghum_harvest_mid==1,1,0),Millet_gs_mid=ifelse(Millet_plant_mid-Millet_harvest_mid<0 | Millet_harvest_mid==1,1,0),Rice_gs_mid=ifelse(Rice_plant_mid-Rice_harvest_mid<0 | Rice_harvest_mid==1,1,0),Wheat.Winter_gs_mid=ifelse(Wheat.Winter_plant_mid-Wheat.Winter_harvest_mid<0 | Wheat.Winter_harvest_mid==1,1,0))]

gs_srt_dt <- merge(plant_end_dt[],harvest_srt_dt,by=c("x","y","mo","Month"))
gs_srt_dt[,`:=`(Maize_gs_srt=ifelse(Maize_plant_end-Maize_harvest_srt<0 | Maize_harvest_srt==1,1,0),Sorghum_gs_srt=ifelse(Sorghum_plant_end-Sorghum_harvest_srt<0 | Sorghum_harvest_srt==1,1,0),Millet_gs_srt=ifelse(Millet_plant_end-Millet_harvest_srt<0 | Millet_harvest_srt==1,1,0),Rice_gs_srt=ifelse(Rice_plant_end-Rice_harvest_srt<0 | Rice_harvest_srt==1,1,0),Wheat.Winter_gs_srt=ifelse(Wheat.Winter_plant_end-Wheat.Winter_harvest_srt<0 | Wheat.Winter_harvest_srt==1,1,0))]

gs_lng_dt <- merge(plant_srt_dt[],harvest_end_dt,by=c("x","y","mo","Month"))
gs_lng_dt[,`:=`(Maize_gs_lng=ifelse(Maize_plant_srt-Maize_harvest_end<0 | Maize_harvest_end==1,1,0),Sorghum_gs_lng=ifelse(Sorghum_plant_srt-Sorghum_harvest_end<0 | Sorghum_harvest_end==1,1,0),Millet_gs_lng=ifelse(Millet_plant_srt-Millet_harvest_end<0 | Millet_harvest_end==1,1,0),Rice_gs_lng=ifelse(Rice_plant_srt-Rice_harvest_end<0 | Rice_harvest_end==1,1,0),Wheat.Winter_gs_lng=ifelse(Wheat.Winter_plant_srt-Wheat.Winter_harvest_end<0 | Wheat.Winter_harvest_end==1,1,0))]


save(plant_srt_dt,plant_mid_dt,plant_end_dt,harvest_srt_dt,harvest_mid_dt,harvest_end_dt,gs_srt_dt,gs_mid_dt,gs_lng_dt,file="data/calendar.RData")



