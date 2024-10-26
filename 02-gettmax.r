library(data.table)
library(sf)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(ncdf4)

rm(list=ls())
gc()

# you will need to have downloaded relevant .nc files for the
# maximum temperature from the CPC Global Unified Temperature 
# website available at
# https://psl.noaa.gov/data/gridded/data.cpc.globaltemp.html
# and stored these .nc files in a folder called tmax that is 
# a subfolder within the same folder where the present R file is;
# then proceed as follows.

# load the Africa map
africa <- ne_countries(continent="Africa",returnclass="sf",scale="large")
africa <- st_set_crs(africa,"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

"%!in%" <- Negate("%in%")

yr_vec <- 1979:2024

yrs <- vector("list",length(yr_vec))

for(r in 1:length(yr_vec)){
  
  temperature_nc <- brick(paste0("tmax/tmax.",yr_vec[r],".nc"))
  
  lst <- vector("list",nlayers(temperature_nc))
  
  for(i in 1:nlayers(temperature_nc)){
    
    p01 <- subset(temperature_nc,i)
    p02 <- rotate(p01)
    p03 <- crop(p02,africa)
    dm_t <- data.table(rasterToPoints(p03))
    dm_t[,`:=`(year=substr(colnames(dm_t)[3],2,5),mo=substr(colnames(dm_t)[3],7,8),dy=substr(colnames(dm_t)[3],10,11))]
    colnames(dm_t)[3] <- "tmax"
    dm_t <- dm_t[order(x,y)]
    lst[[i]] <- dm_t
    
  }
  
  temperature_dt <- data.table(Reduce(rbind,lst))
  temperature_dt <- temperature_dt[,.(tmax=mean(tmax)),by=.(x,y,year,mo)]
  
  yrs[[r]] <- temperature_dt
  
  print(r)
  
}

tmax_dt <- data.table(Reduce(rbind,yrs))

save(tmax_dt,file="data/temperature.RData")


