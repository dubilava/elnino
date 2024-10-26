library(data.table)
library(sf)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(ncdf4)

rm(list=ls())
gc()

# you will need to have downloaded relevant .nc files from the
# CPC Global Unified Gauge-Based Analysis of Daily Precipitation
# website available at 
# https://psl.noaa.gov/data/gridded/data.cpc.globalprecip.html
# and stored these .nc files in a folder called prec that is 
# a subfolder within the same folder where the present R file is;
# then proceed as follows.

# load the Africa map
africa <- ne_countries(continent="Africa",returnclass="sf",scale="large")
africa <- st_set_crs(africa,"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

"%!in%" <- Negate("%in%")

yr_vec <- 1979:2024

yrs <- vector("list",length(yr_vec))

for(r in 1:length(yr_vec)){
  
  precipitation_nc <- brick(paste0("prec/precip.",yr_vec[r],".nc"))
  
  lst <- vector("list",nlayers(precipitation_nc))
  
  for(i in 1:nlayers(precipitation_nc)){
    
    p01 <- subset(precipitation_nc,i)
    p02 <- rotate(p01)
    p03 <- crop(p02,africa)
    dm_p <- data.table(rasterToPoints(p03))
    dm_p[,`:=`(year=substr(colnames(dm_p)[3],2,5),mo=substr(colnames(dm_p)[3],7,8),dy=substr(colnames(dm_p)[3],10,11))]
    colnames(dm_p)[3] <- "prec"
    dm_p <- dm_p[order(x,y)]
    lst[[i]] <- dm_p
    
  }
  
  precipitation_dt <- data.table(Reduce(rbind,lst))
  precipitation_dt <- precipitation_dt[,.(prec=mean(prec)),by=.(x,y,year,mo)]
  
  yrs[[r]] <- precipitation_dt
  
  print(r)
  
}

prec_dt <- data.table(Reduce(rbind,yrs))

save(prec_dt,file="data/precipitation.RData")

