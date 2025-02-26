### phytoSTOICH
### Seston C:N:P validation
### DG, September 2024

# set wd
rm(list = ls())
setwd("C:/Users/DanielGschwentner/Documents/GitHub/getting_droopy_with_it_woodstoich24")
###############################################################################
# Setup

#load packages
pck <- c("deSolve", "tidyverse", "cowplot", "ggthemes","ggpubr")
lapply(pck, require, character.only = T)
theme_set(theme_pubr() + theme(legend.position = "bottom"))

# Load algae parameters
source("algae_param_vctrs.R")

# Load models
# static model with fixed stoichiometry
source("models/static_liebig_zmix.R") # model 1 in Carly's framework
# Droop model
source("models/dynamic_liebig_zmix.R") 
# set timesteps
times <- 1:2000 # for troubleshooting, initial runs

# load data
NEON_Flat <- read_csv("NEON_Flathead_forcing.csv")
NEON_Flat <- as.data.frame(NEON_Flat)

################################################################################

### base predictions w. median values

# dynamic model
(start <- Sys.time())
NEON.dynamic <-  lapply(list(dynamic.algae, dynamic.diatoms, dynamic.greens, dynamic.cyanos), function(x) {
  params <- x
  lapply(1:nrow(NEON_Flat), function(i) {
    params["Pin"] = NEON_Flat[i, "Pin_mgm-3"]
    params["Nin"] = NEON_Flat[i, "Nin_mgm-3"]
    params["DOC"] = NEON_Flat[i, "DOC_mgL-1"]
    params["z"] = NEON_Flat[i, "mean depth_m"]
    params["SA"] = NEON_Flat[i, "SA_km2"]
    params["HRT"] = NEON_Flat[i, "HRT_years_mixedlayer"]
    y <- c("A1" = 100, "P" =NEON_Flat[i, "Pin_mgm-3"], "N" = NEON_Flat[i, "Nin_mgm-3"], 
           "QP1" = 0.015, "QN1" = 0.1)
    run <- ode(y, times, parms = params, func = dynamic.stoich.zmix)
    return(run[max(times),])
  })
})
(end <- Sys.time())
time.elapsed <- (end - start)/60/60
print(paste0("Time elapsed = ", time.elapsed, " hours!"))

# extract from list and convert to data-frame
NEON.dynamic.average <- as_data_frame(do.call(rbind, NEON.dynamic[[1]]))
NEON.dynamic.diatoms <- as_data_frame(do.call(rbind, NEON.dynamic[[2]]))
NEON.dynamic.greens <- as_data_frame(do.call(rbind, NEON.dynamic[[3]]))
NEON.dynamic.cyanos <- as_data_frame(do.call(rbind, NEON.dynamic[[4]]))

################################################################################

### bind everything together and clean data
NEON.CNP <- bind_rows(NEON.dynamic.average, 
                      NEON.dynamic.diatoms,
                      NEON.dynamic.greens,
                      NEON.dynamic.cyanos)
# add-in information
NEON.CNP$Lake <- rep(NEON_Flat$lake, 4)
NEON.CNP$species <- rep(c("average", "diatoms", "greens", "cyanos"), each = nrow(NEON_Flat))

# some mutating and cleaning
NEON.CNP <- NEON.CNP %>%
  mutate(CN_mass_model = 1/QN1, 
         CP_mass_model = 1/QP1, 
         NP_mass_model = QN1/QP1) %>%
  rename(P_mod_ugL = P, N_mod_ugL = N) %>%
  mutate(species = factor(species, levels = c("average", "diatoms", "greens", "cyanos")))

#############################################################################################
# save for plotting in SigmaPlot
#write_csv(NEON.CNP, "NEON_Flathead_model_seston_estimates.csv")


                   