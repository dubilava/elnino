library(data.table)
library(acled.api)

rm(list=ls())
gc()

# you will need to have registered to source the data using acled.api package in R 
# refer to the api user guide available at 
# https://www.acleddata.com/wp-content/uploads/dlm_uploads/2017/10/API-User-Guide.pdf

Sys.setenv(ACLED_EMAIL_ADDRESS="<your_email_address>")
Sys.setenv(ACLED_ACCESS_KEY="<your_access_key>")

acled_dt <- as.data.table(acled.api(region=c(1:5),all.variables=T))

save(acled_dt,file="conflict.RData")


# subset of incidents with "agrarian" term in the notes

acled_dt[,check:=str_detect(notes,paste0("\\b",c("farm","farms","farmer","farmers","peasant","peasants","producer","producers","agriculture","agricultural","livestock","animal","animals","cattle","cow","cows","sheep","goat","goats","crop","crops","cereal","cereals","grain","grains","maize","millet","rice","sorghum","wheat","produce","food","harvest","harvesting"),"\\b",collapse="|"))]

subset_dt <- acled_dt[check==T]
subset_dt$check <- NULL

save(subset_dt,file="conflict_agri.RData")

