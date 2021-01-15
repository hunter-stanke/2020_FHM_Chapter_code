
## -------------------------------------------------------
## Examine temporal trends in lodgepole pine mortality
## in Colorado using various design-based estimators
## -------------------------------------------------------


## Would you like to use parallel processing to speed things up?
## If so, check how many cores are available with
parallel::detectCores(logical = FALSE)

## I don't recommend using all of them, e.g., if 4 available, take 3
cores = 3

## Load packages --------------------------------------------------------------------
library(rFIA)
library(here)



## Download/ save the FIA data -------------------------------------------------------
getFIA('CO', dir = here('data/FIA'), load = FALSE)



## Read the data into memory --------------------------------------------------------
co <- readFIA(here('data/FIA'), states = 'CO', cores = 3)



## Estimate mortality with 5 different design-based estimators -----------------------
## Temporally Indifferent (EVALIDator method)
ti <- growMort(co, treeDomain = SPCD == 108, method = 'TI', nCores = cores)

## Simple moving average
sma <- growMort(co, treeDomain = SPCD == 108, method = 'SMA', nCores = cores)

## Linear moving average
lma <- growMort(co, treeDomain = SPCD == 108, method = 'LMA', nCores = cores)

## Exponential moving average (lambda = .5)
ema <- growMort(co, treeDomain = SPCD == 108, method = 'EMA', lambda = .5, nCores = cores)

## Annual panels, no combination
annual <- growMort(co, treeDomain = SPCD == 108, method = 'ANNUAL', nCores = cores)




## Estimate mortality with the exponential moving average estimators ------------------
## variying the values of lambda
emaLambda <- growMort(db = co,
                      treeDomain = SPCD == 108,
                      method = 'EMA', 
                      lambda = seq(0.05, 0.95, by = 0.05), 
                      nCores = cores)




## Save our results ------------------------------------------------------------------
write.csv(ti, here('results/mort_ti.csv'))
write.csv(sma, here('results/mort_sma.csv'))
write.csv(lma, here('results/mort_lma.csv'))
write.csv(ema, here('results/mort_ema.csv'))
write.csv(annual, here('results/mort_annual.csv'))
write.csv(emaLambda, here('results/mort_ema_ribbon.csv'))

