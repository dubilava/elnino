library(data.table)
library(fixest)
library(zoo)
library(ggplot2)
library(ggforce)
library(ggnewscale)
library(ggbreak)
library(cowplot)
library(viridis)
library(stringr)
library(sf)
library(haven)
library(splines)
library(rnaturalearth)
library(rnaturalearthdata)

rm(list=ls())
gc()

"%!in%" <- Negate("%in%")

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


# load the acled data
load("data/climatecropsconflict.RData")



# Fig.2 - conflict ----

conflict_dt <- climatecropsconflict_dt[,.(incidents=sum(incidents)),by=.(x,y,type,civilian_targeting)]

conflict_dt[,conflict:=ifelse(civilian_targeting=="Civilian_targeting","one-sided","two-sided")]

conf_dt <- climatecropsconflict_dt[,.(x,y,year,mo,date,conflict=ifelse(civilian_targeting=="Civilian_targeting","one-sided","two-sided"),incidents,incidence=ifelse(incidents>0,1,0))]

sum_dt <- conf_dt[incidence==1,.(incidents=sum(incidents),cells=.N),by=.(date,conflict)]
sum_dt <- sum_dt[,cells:=as.numeric(cells)]
sum_dt <- sum_dt[,`incidents/cell`:=incidents/cells]

sum_dt <- sum_dt[order(conflict,date)]

sum_lg <- melt(sum_dt,id.vars=c("conflict","date"))


# generate the maps
gg_conflict1 <- ggplot(data = africa) +
  geom_sf(color=NA,fill="white",linewidth=.1)+
  geom_point(data=conflict_dt[conflict=="one-sided" & incidents>0],aes(x=x,y=y,size=incidents),color="black",fill="black",shape=22,stroke=.1)+
  scale_size_continuous(range=c(.1,.6),breaks=c(5,50,500,5000),trans="log",guide="none")+
  geom_sf(data=lakes,fill="gray",color="gray",linewidth=.1)+
  geom_sf(data=rivers,fill="gray",color="gray",linewidth=.1)+
  geom_sf(color="gray",fill=NA,linewidth=.2)+
  labs(subtitle="(a) One-sided")+
  coord_sf(ylim=c(-35,35),xlim=c(-17,51))+
  theme_void()+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank(),axis.title = element_blank(),axis.text = element_blank(),legend.title=element_text(size=9),legend.text=element_text(size=8,hjust=0.5),legend.position="inside",legend.position.inside=c(.2,.3),legend.key.height=unit(.15,"in"),panel.background=element_rect(fill="white",color=NA),plot.background=element_rect(fill="white",color=NA),plot.title.position="plot")

gg_conflict2 <- ggplot(data = africa) +
  geom_sf(color=NA,fill="white",linewidth=.1)+
  geom_point(data=conflict_dt[conflict=="two-sided" & incidents>0],aes(x=x,y=y,size=incidents),color="black",fill="black",shape=22,stroke=.1)+
  scale_size_continuous(range=c(.1,.6),breaks=c(5,50,500,5000),trans="log",guide="none")+
  geom_sf(data=lakes,fill="gray",color="gray",linewidth=.1)+
  geom_sf(data=rivers,fill="gray",color="gray",linewidth=.1)+
  geom_sf(color="gray",fill=NA,linewidth=.2)+
  labs(subtitle="(b) Two-sided")+
  coord_sf(ylim=c(-35,35),xlim=c(-17,51))+
  theme_void()+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank(),axis.title = element_blank(),axis.text = element_blank(),legend.title=element_text(size=9),legend.text=element_text(size=8,hjust=0.5),legend.position="inside",legend.position.inside=c(.2,.3),legend.key.height=unit(.15,"in"),panel.background=element_rect(fill="white",color=NA),plot.background=element_rect(fill="white",color=NA),plot.title.position="plot")


gg_tsconflict1 <- ggplot(sum_dt[conflict=="one-sided"],aes(x=date,y=incidents))+
  geom_line(linewidth=.3,color="black")+
  scale_y_continuous(breaks=seq(0,1200,400))+
  coord_cartesian(xlim=c(as.Date("1997-06-01"),as.Date("2024-05-31")),ylim=c(0,1200))+
  labs(x="Year",y="",subtitle="Incidents")+
  theme_minimal_hgrid()+
  theme(plot.title.position="plot",panel.background = element_rect(color=NA,fill="white"),plot.background = element_rect(color=NA,fill="white"),axis.text=element_text(size=10,color="dimgray"),panel.grid.major.y = element_line(linewidth=.4,linetype=3,color="dimgray"),panel.spacing=unit(1,"lines"),axis.ticks = element_blank(),axis.line.x=element_blank(),axis.title=element_text(size=10),plot.subtitle=element_text(size=10))

gg_tsconflict2 <- ggplot(sum_dt[conflict=="two-sided"],aes(x=date,y=incidents))+
  geom_line(linewidth=.3,color="black")+
  scale_y_continuous(breaks=seq(0,1500,500))+
  coord_cartesian(xlim=c(as.Date("1997-06-01"),as.Date("2024-05-31")),ylim=c(0,1500))+
  labs(x="Year",y="",subtitle="Incidents")+
  theme_minimal_hgrid()+
  theme(plot.title.position="plot",panel.background = element_rect(color=NA,fill="white"),plot.background = element_rect(color=NA,fill="white"),axis.text=element_text(size=10,color="dimgray"),panel.grid.major.y = element_line(linewidth=.4,linetype=3,color="dimgray"),panel.spacing=unit(1,"lines"),axis.ticks = element_blank(),axis.line.x=element_blank(),axis.title=element_text(size=10),plot.subtitle=element_text(size=10))


gg_conflict_combined <- plot_grid(gg_conflict1,gg_conflict2,gg_tsconflict1,gg_tsconflict2,align="none",axis="none",ncol=2,rel_heights=c(19,10))

gg_conflict_combined

ggsave("figures/fig2_conflict.png",gg_conflict_combined,width=6.5,height=6.5*9/16*1.5,dpi="retina")



# Fig.C1 - conflict series ----

ggc <- ggplot(sum_lg,aes(x=date,y=value,group=variable))+
  geom_line(color="black",linewidth=.4)+
  facet_wrap(variable~conflict,ncol=2,scales="free",labeller = label_wrap_gen(multi_line=F))+
  labs(x="",y="")+
  coord_cartesian(xlim=c(as.Date("1997-06-01"),as.Date("2024-05-31")))+
  theme_minimal_hgrid()+
  theme(panel.grid.major=element_line(linetype=3),plot.subtitle=element_text(size=12),axis.ticks=element_blank(),axis.title = element_blank(),axis.text=element_text(size=10),panel.background=element_rect(fill="white",color=NA),plot.background=element_rect(fill="white",color=NA),plot.title.position="plot",strip.text = element_text(size=10,hjust=0))

ggsave("figures/figc1_conflict_series.png",ggc,width=6.5,height=6.5*9/16*1.5,dpi="retina")




# Fig.3 - croplands ----

sub_dt <- unique(climatecropsconflict_dt[civilian_targeting=="Civilian_targeting" & elnino_year==2020,.(tc=mean(tc),gs_tc=mean(gs_tc),gs_tcp=mean(gs_tcp),gs_tct=mean(gs_tct),tc_prec=mean(tc_prec),gs_tc_prec=mean(gs_tc_prec),tc_tmax=mean(tc_tmax),gs_tc_tmax=mean(gs_tc_tmax),area10=area/10000,gs=sum(gs)),by=.(x,y,crop)])

sub0_dt <- sub_dt[,area10:=ifelse(area10>5,5,area10)]
sub0_dt <- sub0_dt[area10>0,.(x,y,area10)]

gg1 <- ggplot(data=africa) +
  geom_sf(color=NA,fill="white",linewidth=.1)+
  geom_point(data=sub_dt[area10>0],aes(x=x,y=y,size=area10),color="black",fill="black",shape=22,stroke=.1,na.rm=T)+
  scale_size_continuous(range=c(.1,0.6),breaks=c(1:5),guide="none")+
  geom_sf(color="gray",fill=NA,linewidth=.2)+
  geom_sf(data=lakes,fill="gray",color="gray",linewidth=.1)+
  geom_sf(data=rivers,fill="gray",color="gray",linewidth=.1)+
  geom_sf(color="gray",fill=NA,linewidth=.2)+
  coord_sf(ylim=c(-35,35),xlim=c(-17,51))+
  labs(subtitle = "(a) Croplands")+
  theme_void()+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank(),axis.title = element_blank(),axis.text = element_blank(),legend.title=element_text(size=9),legend.text=element_text(size=8,hjust=0.5),legend.position="inside",legend.position.inside=c(.2,.3),legend.key.height=unit(.15,"in"),panel.background=element_rect(fill="white",color=NA),plot.background=element_rect(fill="white",color=NA),plot.title.position="plot")

gg2 <- ggplot(data=africa) +
  geom_sf(color=NA,fill="white",linewidth=.1)+
  geom_point(data=sub_dt[gs_tc>0],aes(x=x,y=y,size=gs_tc),color="black",fill="black",shape=22,stroke=.1,na.rm=T)+
  scale_size_continuous(range=c(.1,0.6),breaks=seq(0,1,.2),guide="none")+
  geom_sf(color="gray",fill=NA,linewidth=.2)+
  geom_sf(data=lakes,fill="gray",color="gray",linewidth=.1)+
  geom_sf(data=rivers,fill="gray",color="gray",linewidth=.1)+
  geom_sf(color="gray",fill=NA,linewidth=.2)+
  coord_sf(ylim=c(-35,35),xlim=c(-17,51))+
  labs(subtitle = "(b) Teleconnections")+
  theme_void()+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank(),axis.title = element_blank(),axis.text = element_blank(),legend.title=element_text(size=9),legend.text=element_text(size=8,hjust=0.5),legend.position="inside",legend.position.inside=c(.2,.3),legend.key.height=unit(.15,"in"),panel.background=element_rect(fill="white",color=NA),plot.background=element_rect(fill="white",color=NA),plot.title.position="plot")



gg1h <- ggplot(sub0_dt,aes(x=area10)) +
  geom_histogram(binwidth=.25,boundary=0,color="white",fill="dimgray",linewidth=.2,na.rm=T)+
  scale_x_continuous(labels=c("0","1","2","3","4","5+"))+
  labs(y="",x="Cropland area (10,000 ha)",subtitle="Cells")+
  theme_minimal_hgrid()+
  theme(plot.title.position="plot",panel.background = element_rect(color=NA,fill="white"),plot.background = element_rect(color=NA,fill="white"),axis.text=element_text(size=10,color="dimgray"),panel.grid.major.y = element_line(linewidth=.4,linetype=3,color="dimgray"),panel.spacing=unit(1,"lines"),axis.ticks = element_blank(),axis.line.x=element_blank(),axis.title=element_text(size=10),plot.subtitle=element_text(size=10))

gg2h <- ggplot(sub_dt[gs_tc>0],aes(x=gs_tc)) +
  geom_histogram(binwidth=.05,boundary=0,color="white",fill="dimgray",linewidth=.2,na.rm=T)+
  scale_x_continuous(breaks=seq(0,1,.2))+
  labs(y="",x="Proportion of months in a growing season",subtitle="Cells")+
  theme_minimal_hgrid()+
  theme(plot.title.position="plot",panel.background = element_rect(color=NA,fill="white"),plot.background = element_rect(color=NA,fill="white"),axis.text=element_text(size=10,color="dimgray"),panel.grid.major.y = element_line(linewidth=.4,linetype=3,color="dimgray"),panel.spacing=unit(1,"lines"),axis.ticks = element_blank(),axis.line.x=element_blank(),axis.title=element_text(size=10),plot.subtitle=element_text(size=10))


gg_croplands_combined <- plot_grid(gg1,gg2,gg1h,gg2h,align="none",axis="none",ncol=2,rel_heights=c(19,10))

ggsave("figures/fig3_croplands.png",gg_croplands_combined,width=6.5,height=6.5*9/16*1.5,dpi="retina")







# Fig.C2 - crops and harvest ----


har_dt <- unique(climatecropsconflict_dt[civilian_targeting=="Civilian_targeting" & elnino_year==2020 & harv==1,.(area10=area/10000),by=.(x,y,month)])

har_dt <- merge(sub_dt[,.(x,y,crop)],har_dt,by=c("x","y"))
har_dt$crop <- factor(har_dt$crop,levels=c("rice","maize","wheat","millet","sorghum"))
har_dt$month <- factor(har_dt$month,levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

gg1 <- ggplot(data=africa) +
  geom_sf(color=NA,fill="white",linewidth=.1)+
  geom_point(data=har_dt,aes(x=x,y=y,color=crop,fill=crop,size=area10),shape=22,stroke=.1,na.rm=T)+
  scale_color_viridis_d(option="turbo",limits=c("rice","maize","wheat","millet","sorghum"),begin=.2,guide="none")+
  scale_fill_viridis_d(option="turbo",limits=c("rice","maize","wheat","millet","sorghum"),begin=.2)+
  scale_size_continuous(range=c(.2,.8),guide="none")+
  geom_sf(data=lakes,fill="gray",color="gray",linewidth=.1)+
  geom_sf(data=rivers,fill="gray",color="gray",linewidth=.1)+
  geom_sf(color="gray",fill=NA,linewidth=.2)+
  labs(subtitle="(a) Crop prevalence")+
  coord_sf(ylim=c(-35,35),xlim=c(-17,51))+
  theme_void()+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank(),axis.title = element_blank(),axis.text = element_blank(),legend.position="none",panel.background=element_rect(fill="white",color=NA),plot.background=element_rect(fill="white",color=NA),plot.title.position="plot")

gg2 <- ggplot(data=africa) +
  geom_sf(color=NA,fill="white",linewidth=.1)+
  geom_point(data=har_dt,aes(x=x,y=y,color=month,fill=month,size=area10),shape=22,stroke=.1,na.rm=T)+
  scale_color_viridis_d(option="turbo",limits=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),begin=.2,guide="none")+
  scale_fill_viridis_d(option="turbo",limits=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),begin=.2)+
  scale_size_continuous(range=c(.2,.8),guide="none")+
  geom_sf(data=lakes,fill="gray",color="gray",linewidth=.1)+
  geom_sf(data=rivers,fill="gray",color="gray",linewidth=.1)+
  geom_sf(color="gray",fill=NA,linewidth=.2)+
  labs(subtitle="(b) Harvest month")+
  coord_sf(ylim=c(-35,35),xlim=c(-17,51))+
  theme_void()+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank(),axis.title = element_blank(),axis.text = element_blank(),legend.position="none",panel.background=element_rect(fill="white",color=NA),plot.background=element_rect(fill="white",color=NA),plot.title.position="plot")


gg1h <- ggplot(har_dt,aes(x=crop,color=crop,fill=crop)) +
  geom_bar()+
  scale_color_viridis_d(option="turbo",limits=c("rice","maize","wheat","millet","sorghum"),begin=.2,guide="none")+
  scale_fill_viridis_d(option="turbo",limits=c("rice","maize","wheat","millet","sorghum"),begin=.2,guide="none")+
  labs(subtitle="Cells")+
  theme(axis.ticks=element_blank(),axis.title=element_blank(),axis.text.x=element_text(color="black",size=8,hjust=0.5),axis.text.y=element_text(color="black",size=8,hjust=1),plot.subtitle=element_text(size=10),panel.background=element_rect(fill="white",color=NA),plot.background=element_rect(fill="white",color=NA),plot.margin=unit(c(2,8,0,5),"pt"),plot.title.position = "plot")

gg2h <- ggplot(har_dt,aes(x=month,color=month,fill=month)) +
  geom_bar()+
  scale_color_viridis_d(option="turbo",limits=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),guide="none")+
  scale_fill_viridis_d(option="turbo",limits=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),guide="none")+
  scale_x_discrete(drop=F)+
  labs(subtitle="Cells")+
  theme(axis.ticks=element_blank(),axis.title=element_blank(),axis.text.x=element_text(color="black",size=8,hjust=0.5),axis.text.y=element_text(color="black",size=8,hjust=1),plot.subtitle=element_text(size=10),panel.background=element_rect(fill="white",color=NA),plot.background=element_rect(fill="white",color=NA),plot.margin=unit(c(2,8,0,5),"pt"),plot.title.position = "plot")


gg_harvest_combined <- plot_grid(gg1,gg2,gg1h,gg2h,align="none",axis="none",ncol=2,rel_heights=c(19,10))

ggsave("figures/figc2_cropsharvest.png",gg_harvest_combined,width=6.5,height=6.5*9/16*1.5,dpi="retina")




# Fig.4 - enso/oni ----

oni_dt <- fread("https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt")

oni_dt <- oni_dt[YR%in%1989:2024]
oni_dt$MO <- rep(1:12,times=length(unique(oni_dt$YR)))[1:nrow(oni_dt)]
oni_dt[,DATE:=as.Date(paste0(YR,"-",MO,"-01"))]
oni_dt <- oni_dt[DATE>="1997-06-01" & DATE<="2024-05-01"]

oni_dt[,EVENT:=ifelse(ANOM>=1,"El Nino",ifelse(ANOM<=-1,"La Nina","Neutral"))]

oni_dt$EVENT <- factor(oni_dt$EVENT,levels=unique(oni_dt$EVENT))

gg_enso <- ggplot()+
  geom_line(data=oni_dt,aes(x=DATE,y=ANOM),linewidth=.6,color="black")+
  geom_point(data=oni_dt[SEAS=="NDJ"],aes(x=DATE,y=ANOM),color=ifelse(oni_dt[SEAS=="NDJ"]$EVENT=="Neutral","black","white"),fill=ifelse(oni_dt[SEAS=="NDJ"]$EVENT=="Neutral","white","black"),shape=21,size=2,stroke=.4)+
  scale_color_viridis_d(option="turbo",begin=.2,end=.9,direction=-1,guide="none")+
  scale_fill_viridis_d(option="turbo",begin=.2,end=.9,direction=-1,guide="none")+
  coord_cartesian(ylim=c(-2,3))+
  labs(x="Year",y="",subtitle=expression(degree*"C"))+
  theme_minimal_hgrid()+
  theme(plot.title.position="plot",panel.background = element_rect(color=NA,fill="white"),plot.background = element_rect(color=NA,fill="white"),axis.text=element_text(size=10,color="dimgray"),panel.grid.major.y = element_line(linewidth=.4,linetype=3,color="dimgray"),panel.spacing=unit(1,"lines"),axis.ticks = element_blank(),axis.line.x=element_blank(),axis.title=element_text(size=12),plot.subtitle=element_text(size=12))

ggsave("figures/fig4_enso.png",gg_enso,width=6.5,height=6.5*9/16,dpi="retina")



# Fig.A1 - teleconnections ----

sub1_dt <- unique(climatecropsconflict_dt[type=="Political violence" & civilian_targeting=="Civilian_targeting" & elnino_year==2020,.(gs_tcp=mean(gs_tcp),gs_tct=mean(gs_tct),area10=area/10000),by=.(x,y)])

gg1 <- ggplot(data=africa) +
  geom_sf(color=NA,fill="white",linewidth=.1)+
  geom_point(data=sub1_dt[gs_tcp>0 & area10>0],aes(x=x,y=y,color=gs_tcp,fill=gs_tcp,size=area10),shape=22,stroke=.1,na.rm=T)+
  geom_point(data=sub1_dt[gs_tcp==0 & area10>0],aes(x=x,y=y,size=area10),color="darkgray",fill="white",shape=22,stroke=.1,na.rm=T)+
  scale_color_viridis(option="B",direction=-1,begin=.3,end=.9,breaks=seq(0,1,.2),limits=c(0,1))+
  scale_fill_viridis(option="B",direction=-1,begin=.3,end=.9,breaks=seq(0,1,.2),limits=c(0,1),guide="none")+
  scale_size_continuous(range=c(.2,.8),guide="none")+
  geom_sf(color="gray",fill=NA,linewidth=.2)+
  geom_sf(data=lakes,fill="gray",color="gray",linewidth=.1)+
  geom_sf(data=rivers,fill="gray",color="gray",linewidth=.1)+
  geom_sf(color="gray",fill=NA,linewidth=.2)+
  coord_sf(ylim=c(-35,35),xlim=c(-17,51))+
  labs(subtitle = "(a) Precipitation")+
  theme_void()+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank(),axis.title = element_blank(),axis.text = element_blank(),legend.title=element_text(size=9),legend.text=element_text(size=8,hjust=0.5),legend.position="bottom",legend.key.height=unit(.1,"in"),legend.key.width=unit(.55,"in"),panel.background=element_rect(fill="white",color=NA),plot.background=element_rect(fill="white",color=NA))+
  guides(color=guide_colorbar(title="Proportion of months in a growing season",title.position="top",label.position="bottom",nrow=1))

gg2 <- ggplot(data=africa) +
  geom_sf(color=NA,fill="white",linewidth=.1)+
  geom_point(data=sub1_dt[gs_tct>0 & area10>0],aes(x=x,y=y,color=gs_tct,fill=gs_tct,size=area10),shape=22,stroke=.1,na.rm=T)+
  geom_point(data=sub1_dt[gs_tct==0 & area10>0],aes(x=x,y=y,size=area10),color="darkgray",fill="white",shape=22,stroke=.1,na.rm=T)+
  scale_color_viridis(option="B",direction=-1,begin=.3,end=.9,breaks=seq(0,1,.2),limits=c(0,1))+
  scale_fill_viridis(option="B",direction=-1,begin=.3,end=.9,breaks=seq(0,1,.2),limits=c(0,1),guide="none")+
  scale_size_continuous(range=c(.2,.8),guide="none")+
  geom_sf(color="gray",fill=NA,linewidth=.2)+
  geom_sf(data=lakes,fill="gray",color="gray",linewidth=.1)+
  geom_sf(data=rivers,fill="gray",color="gray",linewidth=.1)+
  geom_sf(color="gray",fill=NA,linewidth=.2)+
  coord_sf(ylim=c(-35,35),xlim=c(-17,51))+
  labs(subtitle = "(b) Temperature")+
  theme_void()+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank(),axis.title = element_blank(),axis.text = element_blank(),legend.title=element_text(size=9),legend.text=element_text(size=8,hjust=0.5),legend.position="bottom",legend.key.height=unit(.1,"in"),legend.key.width=unit(.55,"in"),panel.background=element_rect(fill="white",color=NA),plot.background=element_rect(fill="white",color=NA))+
  guides(color=guide_colorbar(title="Proportion of months in a growing season",title.position="top",label.position="bottom",nrow=1))


gg1h <- ggplot(sub1_dt[gs_tcp>0],aes(x=gs_tcp)) +
  geom_histogram(binwidth=.05,color="white",fill="darkgray",na.rm=T,linewidth=.2,boundary=0.05)+
  scale_x_continuous(breaks=seq(0,1,.2),limits=c(0,1))+
  theme(axis.ticks=element_blank(),axis.title=element_blank(),axis.text.x=element_text(color="black",size=8,hjust=0.5),axis.text.y=element_blank(),panel.background=element_rect(fill="white",color=NA),title=element_text(color="black",size=9),plot.background=element_rect(fill="white",color=NA),plot.margin=unit(c(2,8,0,5),"pt"))

gg2h <- ggplot(sub1_dt[gs_tct>0],aes(x=gs_tct)) +
  geom_histogram(binwidth=.05,color="white",fill="darkgray",na.rm=T,linewidth=.2,boundary=0.05)+
  scale_x_continuous(breaks=seq(0,1,.2),limits=c(0,1))+
  theme(axis.ticks=element_blank(),axis.title=element_blank(),axis.text.x=element_text(color="black",size=8,hjust=0.5),axis.text.y=element_blank(),panel.background=element_rect(fill="white",color=NA),title=element_text(color="black",size=9),plot.background=element_rect(fill="white",color=NA),plot.margin=unit(c(2,8,0,5),"pt"))

gg_tc <- plot_grid(gg1,gg2,gg1h,gg2h,align="v",axis="none",ncol=2,rel_heights=c(5,2))

ggsave("figures/figa1_teleconnections.png",gg_tc,width=6.5,height=6.5*9/16*1.5,dpi="retina")



# Fig C3 - weather impact ----

# Fig.4 - teleconnections ----

gg1 <- ggplot(data=africa) +
  geom_sf(color=NA,fill="white",linewidth=.1)+
  geom_point(data=sub_dt[gs_tc>0 & area10>0],aes(x=x,y=y,color=gs_tc_prec,fill=gs_tc_prec,size=area10),shape=22,stroke=.1,na.rm=T)+
  geom_point(data=sub_dt[gs_tc==0 & area10>0],aes(x=x,y=y,size=area10),color="darkgray",fill="white",shape=22,stroke=.1,na.rm=T)+
  scale_color_viridis(option="turbo",direction=-1,begin=.2,breaks=seq(-.8,.8,.2),limits=c(-.8,.8))+
  scale_fill_viridis(option="turbo",direction=-1,begin=.2,breaks=seq(-.8,.8,.2),limits=c(-.8,.8),guide="none")+
  scale_size_continuous(range=c(.2,.8),guide="none")+
  geom_sf(color="gray",fill=NA,linewidth=.2)+
  geom_sf(data=lakes,fill="gray",color="gray",linewidth=.1)+
  geom_sf(data=rivers,fill="gray",color="gray",linewidth=.1)+
  geom_sf(color="gray",fill=NA,linewidth=.2)+
  coord_sf(ylim=c(-35,35),xlim=c(-17,51))+
  labs(subtitle = "(a) Precipitation impact")+
  theme_void()+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank(),axis.title = element_blank(),axis.text = element_blank(),legend.title=element_text(size=9),legend.text=element_text(size=8,hjust=0.5),legend.position="bottom",legend.key.height=unit(.1,"in"),legend.key.width=unit(.55,"in"),panel.background=element_rect(fill="white",color=NA),plot.background=element_rect(fill="white",color=NA))+
  guides(color=guide_colorbar(title="Growing season change in precipitation (mm/day)",title.position="top",label.position="bottom",nrow=1))


gg2 <- ggplot(data=africa) +
  geom_sf(color=NA,fill="white",linewidth=.1)+
  geom_point(data=sub_dt[gs_tc>0 & area10>0],aes(x=x,y=y,color=gs_tc_tmax,fill=gs_tc_tmax,size=area10),shape=22,stroke=.1,na.rm=T)+
  geom_point(data=sub_dt[gs_tc==0 & area10>0],aes(x=x,y=y,size=area10),color="darkgray",fill="white",shape=22,stroke=.1,na.rm=T)+
  scale_color_viridis(option="turbo",direction=1,begin=.2,breaks=seq(-.8,.8,.2),limits=c(-.8,.8))+
  scale_fill_viridis(option="turbo",direction=1,begin=.2,breaks=seq(-.8,.8,.2),limits=c(-.8,.8),guide="none")+
  scale_size_continuous(range=c(.2,.8),guide="none")+
  geom_sf(color="gray",fill=NA,linewidth=.2)+
  geom_sf(data=lakes,fill="gray",color="gray",linewidth=.1)+
  geom_sf(data=rivers,fill="gray",color="gray",linewidth=.1)+
  geom_sf(color="gray",fill=NA,linewidth=.2)+
  coord_sf(ylim=c(-35,35),xlim=c(-17,51))+
  labs(subtitle = "(b) Temperature impact")+
  theme_void()+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank(),axis.title = element_blank(),axis.text = element_blank(),legend.title=element_text(size=9),legend.text=element_text(size=8,hjust=0.5),legend.position="bottom",legend.key.height=unit(.1,"in"),legend.key.width=unit(.55,"in"),panel.background=element_rect(fill="white",color=NA),plot.background=element_rect(fill="white",color=NA))+
  guides(color=guide_colorbar(title="Growing season change in max temperature (C)",title.position="top",label.position="bottom",nrow=1))


gg1h <- ggplot(sub_dt[gs_tc>0],aes(x=gs_tc_prec)) +
  geom_histogram(binwidth=.1,color="white",fill="darkgray",na.rm=T,linewidth=.2,boundary=.1)+
  scale_x_continuous(breaks=seq(-.8,.8,.2),limits=c(-.8,.8))+
  theme(axis.ticks=element_blank(),axis.title=element_blank(),axis.text.x=element_text(color="black",size=8,hjust=0.5),axis.text.y=element_blank(),panel.background=element_rect(fill="white",color=NA),title=element_text(color="black",size=9),plot.background=element_rect(fill="white",color=NA),plot.margin=unit(c(2,8,0,5),"pt"))

gg2h <- ggplot(sub_dt[gs_tc>0],aes(x=gs_tc_tmax)) +
  geom_histogram(binwidth=.1,color="white",fill="darkgray",na.rm=T,linewidth=.2,boundary=.1)+
  scale_x_continuous(breaks=seq(-.8,.8,.2),limits=c(-.8,.8))+
  theme(axis.ticks=element_blank(),axis.title=element_blank(),axis.text.x=element_text(color="black",size=8,hjust=0.5),axis.text.y=element_blank(),panel.background=element_rect(fill="white",color=NA),title=element_text(color="black",size=9),plot.background=element_rect(fill="white",color=NA),plot.margin=unit(c(2,8,0,5),"pt"))


gg_tc <- plot_grid(gg1,gg2,gg1h,gg2h,align="v",axis="none",ncol=2,rel_heights=c(5,2))

ggsave("figures/figc3_weatherimpact.png",gg_tc,width=6.5,height=6.5*9/16*1.5,dpi="retina")



