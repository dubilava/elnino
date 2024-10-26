library(data.table)
library(ggplot2)
library(sp)
library(sf)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(ncdf4)
library(R.utils)
library(stars)

rm(list=ls())
gc()

# you will need to have downloaded the relevant zipped folder from 
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/PRFF8V
# and stored it in a folder called SPAM within the same folder as the present R file.
# unzip/extract all files from this zipped folder (called dataverse files); 
# then proceed as follows.

# load the Africa map
africa <- ne_countries(continent="Africa",returnclass="sf",scale="large")
africa <- st_set_crs(africa,"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

"%!in%" <- Negate("%in%")

maize_a <- raster("spam2010v2r0_global_phys_area.geotiff/spam2010V2r0_global_A_MAIZ_A.tif")
sorghum_a <- raster("spam2010v2r0_global_phys_area.geotiff/spam2010V2r0_global_A_SORG_A.tif")
milletp_a <- raster("spam2010v2r0_global_phys_area.geotiff/spam2010V2r0_global_A_PMIL_A.tif")
millets_a <- raster("spam2010v2r0_global_phys_area.geotiff/spam2010V2r0_global_A_SMIL_A.tif")
rice_a <- raster("spam2010v2r0_global_phys_area.geotiff/spam2010V2r0_global_A_RICE_A.tif")
wheat_a <- raster("spam2010v2r0_global_phys_area.geotiff/spam2010V2r0_global_A_WHEA_A.tif")


maize_a05 <- aggregate(maize_a,fact=6,fun=sum)
sorghum_a05 <- aggregate(sorghum_a,fact=6,fun=sum)
milletp_a05 <- aggregate(milletp_a,fact=6,fun=sum)
millets_a05 <- aggregate(millets_a,fact=6,fun=sum)
rice_a05 <- aggregate(rice_a,fact=6,fun=sum)
wheat_a05 <- aggregate(wheat_a,fact=6,fun=sum)

mp <- rasterToPoints(maize_a05)
maize_dt <- data.table(round(mp,2))
colnames(maize_dt) <- c("x","y","area_maize")

sp <- rasterToPoints(sorghum_a05)
sorghum_dt <- data.table(round(sp,2))
colnames(sorghum_dt) <- c("x","y","area_sorghum")

mpp <- rasterToPoints(milletp_a05)
milletp_dt <- data.table(round(mpp,2))
colnames(milletp_dt) <- c("x","y","area_milletp")

msp <- rasterToPoints(millets_a05)
millets_dt <- data.table(round(msp,2))
colnames(millets_dt) <- c("x","y","area_millets")

millet_dt <- merge(milletp_dt,millets_dt,by=c("x","y"))
millet_dt[,area_millet:=area_milletp+area_millets]
millet_dt$area_milletp <- NULL
millet_dt$area_millets <- NULL

rp <- rasterToPoints(rice_a05)
rice_dt <- data.table(round(rp,2))
colnames(rice_dt) <- c("x","y","area_rice")

wp <- rasterToPoints(wheat_a05)
wheat_dt <- data.table(round(wp,2))
colnames(wheat_dt) <- c("x","y","area_wheat")

spam_dt <- Reduce(function(...) merge(...,by=c("x","y"),all=T),list(maize_dt,sorghum_dt,millet_dt,rice_dt,wheat_dt))

save(spam_dt,file="data/spam.RData")
