library(data.table)
library(sf)
library(raster)
library(stringr)
library(zoo)
library(haven)
library(rnaturalearth)
library(rnaturalearthdata)

rm(list=ls())
gc()

world <- ne_countries(returnclass="sf",scale="large")
africa <- ne_countries(returnclass="sf",scale="large",continent="Africa")

crs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

world <- st_set_crs(world,crs)
africa <- st_set_crs(africa,crs)

"%!in%" <- Negate("%in%")


load("data/climatecrops.RData")

load("../data/conflicts/acled/acled_africa_2024.RData")


# conflict ----
acled_dt <- acled_dt[,.(x=longitude,y=latitude,country,date=event_date,year,type=disorder_type,event=event_type,civilian_targeting,fatalities,actor1,actor2)]

acled_dt[,`:=`(x=floor(as.numeric(x)*2)/2+.25,y=floor(as.numeric(y)*2)/2+.25,mo=month(date))]

acled_dt <- acled_dt[date<as.Date("2024-06-01")]


acled_dt <- acled_dt[,.(incidents=.N,fatalities=sum(fatalities),actors=uniqueN(c(actor1,actor2))),by=.(x,y,country,year,mo,type,civilian_targeting)]

acled_dt[civilian_targeting=="Civilian targeting",actors:=actors-1]


acled_dt[,xy:=paste(x,y,sep=",")]


## I am doing this to ensure that conflict and crop data geo-coordinates match -- super inefficient coding but gets the stuff done (I think)
xy_acled_dt <- acled_dt[,.(xy)]
xy_acled_dt <- unique(xy_acled_dt)

climatecrops_dt[,xy:=paste(x,y,sep=",")]

xy_crops_dt <- climatecrops_dt[,.(xy)]
xy_crops_dt <- unique(xy_crops_dt)

xy_crops_dt$longitude <- as.numeric(unlist(strsplit(as.character(xy_crops_dt$xy),","))[c(T,F)])
xy_crops_dt$latitude <- as.numeric(unlist(strsplit(as.character(xy_crops_dt$xy),","))[c(F,T)])

xy_acled_dt$longitude <- as.numeric(unlist(strsplit(as.character(xy_acled_dt$xy),","))[c(T,F)])
xy_acled_dt$latitude <- as.numeric(unlist(strsplit(as.character(xy_acled_dt$xy),","))[c(F,T)])

xy_acled_dt <- xy_acled_dt[xy %!in% xy_crops_dt$xy]

d <- pointDistance(xy_acled_dt[,.(longitude,latitude)],xy_crops_dt[,.(longitude,latitude)],lonlat=T)

r <- apply(d,1,which.min)

p <- data.table(acled=xy_acled_dt$xy,crops=xy_crops_dt$xy[r])

colnames(p) <- c("xy","xy_crops")

xy_acled_dt <- acled_dt[,.(xy)]
xy_acled_dt <- unique(xy_acled_dt)

acled_xy_dt <- merge(acled_dt,p,by="xy",all.x=T)
acled_xy_dt[!is.na(xy_crops)]$xy <- acled_xy_dt[!is.na(xy_crops)]$xy_crops

acled_xy_dt$longitude <- as.numeric(unlist(strsplit(as.character(acled_xy_dt$xy),","))[c(T,F)])
acled_xy_dt$latitude <- as.numeric(unlist(strsplit(as.character(acled_xy_dt$xy),","))[c(F,T)])

acled_xy_dt$xy_crops <- NULL

acled_dt <- acled_xy_dt[,.(incidents=sum(incidents),fatalities=sum(fatalities),actors=sum(actors)),by=.(xy,x,y,year,mo,longitude,latitude,type,civilian_targeting)]

rm(d)


# making sure the location-periods when no conflict happens are not missing in the dataset
xy_acled_dt <- acled_dt[,.(xy)]
xy_acled_dt <- unique(xy_acled_dt)

xy_crops_dt <- climatecrops_dt[,.(xy)]
xy_crops_dt <- unique(xy_crops_dt)

xy_dt <- merge(xy_acled_dt,xy_crops_dt,all=T)
xy_dt <- unique(xy_dt)
xy <- as.character(xy_dt$xy)

yearmo <- substr(as.character(seq(as.Date("1997-01-01"),as.Date("2024-05-31"),by="month")),1,7)

xy_yearmo <- CJ(xy,yearmo)
xy_yearmo_dt <- unique(xy_yearmo)

xy_yearmo <- merge(xy_yearmo_dt,xy_dt,by="xy")
xy_yearmo[,`:=`(year=as.numeric(substr(yearmo,1,4)),mo=as.numeric(substr(yearmo,6,7)))]


pv_trg_dt <- acled_dt[type=="Political violence" & civilian_targeting=="Civilian targeting"]

pv_not_dt <- acled_dt[type=="Political violence" & civilian_targeting!="Civilian targeting"]


pv_trg_dt <- merge(xy_yearmo,pv_trg_dt,by=c("xy","year","mo"),all.x=T)
pv_not_dt <- merge(xy_yearmo,pv_not_dt,by=c("xy","year","mo"),all.x=T)


pv_trg_dt[,`:=`(type="Political violence",civilian_targeting="Civilian_targeting")]
pv_trg_dt[,c("x","y"):=tstrsplit(xy,",",fixed=T)]

pv_not_dt[,`:=`(type="Political violence",civilian_targeting="")]
pv_not_dt[,c("x","y"):=tstrsplit(xy,",",fixed=T)]


pv_trg_dt$longitude <- NULL
pv_trg_dt$latitude <- NULL

pv_trg_dt[is.na(pv_trg_dt)] <- 0
pv_trg_dt <- pv_trg_dt[,.(incidents=sum(incidents),fatalities=sum(fatalities),actors=sum(actors)),by=.(xy,x,y,yearmo,year,mo,type,civilian_targeting)]


pv_not_dt$longitude <- NULL
pv_not_dt$latitude <- NULL

pv_not_dt[is.na(pv_not_dt)] <- 0
pv_not_dt <- pv_not_dt[,.(incidents=sum(incidents),fatalities=sum(fatalities),actors=sum(actors)),by=.(xy,x,y,yearmo,year,mo,type,civilian_targeting)]


conflict_dt <- rbind(pv_trg_dt,pv_not_dt)#,unrest_dt)

rm(pv_trg_dt)
rm(pv_not_dt)
rm(xy_yearmo)
rm(xy_yearmo_dt)
gc()

conflict_dt[is.na(conflict_dt)] <- 0


conflict_dt[,`:=`(x=as.numeric(x),y=as.numeric(y),mo=as.factor(as.numeric(as.character(mo))))]


# merge with crops ----
climatecrops_dt[,`:=`(x=as.numeric(x),y=as.numeric(y),mo=as.factor(as.numeric(as.character(mo))))]


climateconflict_dt <- merge(conflict_dt,climatecrops_dt,by=c("x","y","xy","year","mo"))

climateconflict_dt <- climateconflict_dt[order(type,civilian_targeting,country,x,y,year,mo)]

climatecropsconflict_dt <- climateconflict_dt[elnino_year!=1996]


# some country renaming for general consistency

climatecropsconflict_dt[country=="Bir Tawil"]$country <- "Sudan"
climatecropsconflict_dt[country=="Central African Rep."]$country <- "Central African Republic"
climatecropsconflict_dt[country=="Dem. Rep. Congo"]$country <- "Democratic Republic of Congo"
climatecropsconflict_dt[country=="Eq. Guinea"]$country <- "Equatorial Guinea"
climatecropsconflict_dt[country=="Gambia"]$country <- "The Gambia"
climatecropsconflict_dt[country=="S. Sudan"]$country <- "South Sudan"
climatecropsconflict_dt[country=="Somaliland"]$country <- "Somalia"
climatecropsconflict_dt[country=="W. Sahara"]$country <- "Morocco"


save(climatecropsconflict_dt,file="data/climatecropsconflict.RData")


