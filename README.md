# Replication Material

This repository hosts data and R scripts for the paper "[Climate, Crops, and Postharvest Conflict](https://doi.org/10.1111/ajae.12504)" published in the American Journal of Agricultural Economics. 


## Data

### Climate

- The *precipitation* and *temperature* data are the NOAA/CPC Global Unified Gauge-Based Analysis of Daily Precipitation available at [https://psl.noaa.gov/data/gridded/data.cpc.globalprecip.html](https://psl.noaa.gov/data/gridded/data.cpc.globalprecip.html) and Global Unified Temperature available at [https://psl.noaa.gov/data/gridded/data.cpc.globaltemp.html](https://psl.noaa.gov/data/gridded/data.cpc.globaltemp.html). 

- The *sea surface temperature* data, which I use as the proxy for El Nino Southern Oscillation, are from the National Oceanic and Atmospheric Administration (NOAA) Climate Prediction Center available at [https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt](https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt).

### Crops

- The *harvest area* data are from the International Food and Policy Research Institute's MapSPAM 2010, obtained from Harvard Dataverse available at: [https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/PRFF8V](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/PRFF8V). For details, see: Yu, Q., You, L., Wood-Sichra, U., Ru, Y., Joglekar, A.K.B., Fritz, S., Xiong, W., Lu, M., Wu, W., and P., Yang. (2020). [A Cultivated Planet in 2010-Part 2: The Global Gridded Agricultural-Production Maps](https://essd.copernicus.org/articles/12/3545/2020/). Earth System Science Data, 12(4), 3545-3572.

- The *harvest calendar* data are from Sacks, W.J., D. Deryng, J.A. Foley, and N. Ramankutty (2010). Crop planting dates: an analysis of global patterns. Global Ecology and Biogeography 19, 607-620, available at: [https://sage.nelson.wisc.edu/data-and-models/datasets/crop-calendar-dataset/netcdf-0-5-degree/](https://sage.nelson.wisc.edu/data-and-models/datasets/crop-calendar-dataset/netcdf-0-5-degree/)
  
### Conflict

- The *conflict* data are from the Armed Conflict Location & Event Data Project (ACLED) available at: [http://acleddata.com](www.acleddata.com). For details, see: Raleigh, C., A. Linke, H. Hegre, and J. Karlsen (2010). Introducing ACLED: "An Armed Conflict Location and Event Dataset: Special Data Feature." Journal of Peace Research 47(5): 651–660.


## R scripts

### To compile the data

This replication material does not include the original source data. However, scripts 01 through 05 provide brief instructions on obtaining these data. The compiled data generated by these scripts is included in the replication materials and can be used to create the final dataset using scripts 06 through 08. The datasets compiled by scripts 06, 07, and 08 are not included here due to their large size. However, they can be replicated in just a few minutes by sequentially running scripts 06, 07, and 08.

- [01-getprec.r](01-getprec.r): download and store precipitation data
- [02-gettmax.r](02-gettmax.r): download and store temperature data
- [03-getcalendar.r](03-getcalendar.r): download and store harvest calendar data
- [04-getspam.r](04-getspam.r): download and store harvest area data
- [05-getacled.r](05-getacled.r): download and store conflict data
- [06-climatecrops.r](06-climatecrops.r): combine climate and crops data
  * the compiled data not stored in the replication material for its large size, but it can be replicated by running this script.
- [07-climatecropsconflict.r](07-climatecropsconflict.r): add conflict data to the combined climate and crops data
  * the compiled data not stored in the replication material for its large size, but it can be replicated by running this script.
- [08-climatecropsconflict_agri.r](07-climatecropsconflict.r): add agrarian conflict data to the combined climate and crops data
  * the compiled data not stored in the replication material for its large size, but it can be replicated by running this script.

### To replicate the results

To replicate the results, first run scripts 06, 07, and 08 to compile the dataset required for the analyses.

- [11-descriptive.r](11-descriptive.r): Figures 2-4 and Online Appendix Figures
- [12-results.r](12-results.r): Tables 1-3, Figure 5
- [13-results_agri.r](13-results_agri.r): Table 4
- [14-robust.r](14-robust.r): Online Appendix Tables


## License

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

