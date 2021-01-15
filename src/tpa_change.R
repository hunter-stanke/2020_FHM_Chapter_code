
## -------------------------------------------------------
## Estimate plot-level live tree abundance for every plot
## visit under the annual FIA design with rFIA.
## -------------------------------------------------------

## Note: data downloaded herein are extremely large (~50GB).
##       Using rFIA's 'out-of-memory' methods, we can still
##       process these data on a standard desktop computer.
##       Here we use a machine w/ 4 cores and ~32GB RAM

## Load packages --------------------------------------------------------------------
library(rFIA)
library(dplyr)
library(tidyr)
library(R2jags)
library(coda)
library(stringr)



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




## Set up our 'remote' FIA.Database --------------------------------------------------
fia <- readFIA(here('data/FIA'), inMemory = FALSE)
## Again, data is about 50GB total, so we can't fit it all in RAM at once
## Instead, we tell readFIA to 'remember' where the data are on our
## machine, and then other rFIA functions will read in the data they need
## state-by-state when they need to




## Compute estimates of plot-level TPA -----------------------------------------------
## Live tree abundance
pltTPA <- tpa(fia, grpBy = c(ECOSUBCD), byPlot = TRUE, treeType = 'live')




## Prep data for model ---------------------------------------------------------------
## Filter out nonforested plot/ times
pltTPA <- filter(pltTPA, PLOT_STATUS_CD == 1)

## Filter out visits w/ 0 TPA - causes computational issues
pltTPA <- filter(pltTPA, BAA > 0)

## ECOSUBCD often has whitespace that causes issues
pltTPA$ECOSUBCD <- str_trim(pltTPA$ECOSUBCD)



## Count the number of times each plot was measured
## as a forested plot, and remove those which have
## yet to be remeasured
nMeas <- pltTPA %>%
  group_by(pltID) %>%
  summarize(n = n())
pltTPA <- pltTPA %>%
  left_join(nMeas, by = 'pltID') %>%
  filter(n > 1)

## We want time to be the same everywhere
## Want to minimize uncertainty around the alpha
## So pick an approximate midpoint to center time
pltTPA <- pltTPA %>%
  group_by(pltID) %>%
  mutate(t = YEAR - 2010,
         ## Add a small constant to account for zero values
         y = log(BAA))

## A unique number identifying each plot location
pltTPA$g <- as.numeric(as.factor(pltTPA$pltID))


## A unique number identifying the spatial unit
## each plot location belongs to
spID <- pltTPA %>%
  ungroup() %>%
  select(g, ECOSUBCD) %>%
  distinct() %>%
  arrange(g) %>%
  mutate(h = as.numeric(as.factor(ECOSUBCD)))

## Set up parameters in a list
data <- list(N = nrow(pltTPA), ## number of obs
             I = nrow(spID), ## Number of unique plot locations
             J = length(unique(spID$h)), ## Number of unique spatial units
             y = pltTPA$y, ## TPA at each measurement
             t = pltTPA$t, ## time since intial measurement
             g = pltTPA$g, ## value linking observation to plot
             h = spID$h) ## value linking plot to spatial unit



### Fit model -------------------------------------------------------------------

# Parameters to estimate
params <- c("unit_beta", "global_beta")

# MCMC settings
ni <- 1000
nc <- 4

# Start Gibbs sampling
jags_mod <- jags(data,
                 parameters.to.save=params,
                 model.file=here("src/densityChangeModel.jag"),
                 n.chains=nc,
                 n.iter=ni)
jags_mod <- autojags(jags_mod, iter = 500)

## Convert to mcmc list
jags_mcmc <- as.mcmc(jags_mod)

chains <- jags_mcmc
## Convert to data.frame
for (i in 1:length(jags_mcmc)){
  chains[[i]] <- as.data.frame(jags_mcmc[[i]])
  names(chains)[i] <- i
}
class(chains) <- 'list'

## Tidy it up
post <- bind_rows(chains) %>%
  pivot_longer(cols = everything(), names_to = 'var', values_to = 'estimate') %>%
  filter(str_detect(var, 'beta')) %>%
  mutate(term = if_else(str_detect(var, 'global'), 'global', 'subsection')) %>%
  mutate(h = readr::parse_number(var),
         h = as.numeric(h)) %>%
  left_join(distinct(spID, h, ECOSUBCD), by = 'h') %>%
  left_join(distinct(ungroup(pltTPA), ECOSUBCD), by = 'ECOSUBCD') %>%
  select(ECOSUBCD, term, estimate) %>%
  ## We log scaled eariler, now exponentiate to get
  ## our original units back
  mutate(estimate = estimate)


## Summarize the posterior
post_sum <- post %>%
  group_by(ECOSUBCD, term) %>%
  summarise(m = median(estimate),
            variance = var(estimate),
            prec = 1 / variance,
            sd = sqrt(variance),
            upper = quantile(estimate, probs = .975),
            lower = quantile(estimate, probs = .025),
            sig = if_else(upper * lower > 0, 1, 0),
            pd = if_else(m > 0,
                         length(estimate[estimate > 0]) / length(estimate),
                         length(estimate[estimate < 0]) / length(estimate)))





## Save results ---------------------------------------------------------------
write.csv(pltTPA, here('results/tpa_change_observations.csv'))
write.csv(post, here('results/tpa_change_posterior_samples.csv'))
write.csv(post_sum, here('results/tpa_change_posterior_summary.csv'))
