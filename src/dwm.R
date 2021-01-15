## -------------------------------------------------------
## Estimate down woody material carbon/ biomass/ volume
## stocks across HUC6 regions in the CONUS, using the most
## recent data from each state
## -------------------------------------------------------

## Note: data downloaded herein are extremely large (~50GB).
##       Using rFIA's 'out-of-memory' methods, we can still
##       process these data on a standard desktop computer.
##       Here we use a machine w/ 4 cores and ~32GB RAM


## Load packages --------------------------------------------------------------------
library(rFIA)
library(rgdal)
library(here)


## Download/ save the FIA data -------------------------------------------------------
## State abbs
allStates <- c('AL', 'AZ', 'AR', 'CA', 'CO', 'CT', 'DE', 'FL', 'GA', 'ID', 'IL', 'IN',
               'IA', 'KS', 'KY', 'LA', 'ME', 'MD', 'MA', 'MI', 'MN', 'MS', 'MO', 'MT',
               'NE', 'NV', 'NH', 'NJ', 'NM', 'NY', 'NC', 'ND', 'OH', 'OK', 'OR', 'PA',
               'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 'VT', 'VA', 'WA', 'WV', 'WI', 'WY')


## Downloading from the FIA Datamart, saving each stat subset in the same directory
## but NOT LOADING it (load = FALSE). Too large to load all at once
## Only need to run this once, i.e., no need to download again for future use
getFIA(allStates, dir = here('data/FIA'), load = FALSE)


## Note: If you haven't see the "here" package before, check out 
##       this link: https://here.r-lib.org/. 
##       Well worth the short learning curve, and makes sharing
##       code much easier. 



## Set up our 'remote' FIA.Database --------------------------------------------------
fia <- readFIA(here('data/FIA'), inMemory = FALSE, nCores = 3)
## Again, data is about 50GB total, so we can't fit it all in RAM at once
## Instead, we tell readFIA to 'remember' where the data are on our
## machine, and then other rFIA functions will read in the data they need
## state-by-state when they need to


## Take a most recent subset ---------------------------------------------------------
fiaMR <- clipFIA(fia)
## Retains only the data necessary to produce estimates for the most recent
## reporting years in each state. Speeds up processing, and allows us to
## 'merge' most recent (but different) reporting years for populations that
## span state boundaries (as HUC units often do)


## Bring in our geospatial data ------------------------------------------------------
## HUC6 clipped to CONUS boundary
huc <- st_read(here('data/GIS/HUC6/HUC6.shp'))



## Compute estimates of dwm stocks ---------------------------------------------------
dwmHUC <- dwm(fiaMR, polys = huc, returnSpatial = TRUE, nCores = 4)
## Here we say compute estimates of dwm stocks within each spatial unit in 'huc',
## treating each spatial unit as a unique population. As we hand 'dwm' the 'fiaMR'
## Remote.FIA.Database, dwm will cleverly read in the most recent data from each state
## one state at a time, process results to the estimation unit level, and then combine
## those results (much smaller than the full dataset) in memory at the end. By
## stating 'returnSpatial = TRUE', we tell 'dwm' to simply append results to the 'huc'
## polygon data, making it easy to produce spatial or spatio-temporal visualizations



## Save the results ---------------------------------------------------
write.csv(dwmHUC, here('results/dwm.csv', row.names = FALSE))
