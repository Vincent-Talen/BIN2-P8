## Copyright (c) 2023 Vincent Talen.
## Licensed under GPLv3. See LICENSE file.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: scenarios.R
##
## Purpose of script: Simulate scenarios using functions from the functions.R module
##
## Author: Vincent Talen
##
## Date Created: 08 Jan 2023
## Date Created: 09 Jan 2023
##
## Email: v.k.talen@st.hanze.nl
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Notes:
##   - Energetic efficiency plot needs to be added (back from creating the model section)
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
#setwd("BIN2-P8/endreport")
#setwd("../..")


# ######### #
#   Libs    #
# ######### #
source("src/functions.R")
source("src/simulateScenario.R")


# ######### #
#   Code    #
# ######### #
# Temperatures to do simulations of
temperatures <- c(5, 10, 15, 20, 25)

# SCENARIO 0: STANDARD SCENARIO ###########################################

# Gammarus mean body mass = 4.26 mgDM
# Annual leaf fall        = 300 gC/m2/an = 300 000 mgC/m2/an
# Gammarus density        = 30 mgDM/m2   = 15 mgC/m2

# Duration of the leaf fall in days
fall_duration_in_days <- 15

# Define scenario values
gamm_indv_mass <- 4.26
leaf_fall <- 300000 / fall_duration_in_days
gamm_start_biomass <- 15

# Get data for temperatures with values of current scenario
scen_df_list <- getScenarioDataList(gamm_indv_mass, leaf_fall, gamm_start_biomass)

# Create plots in an arranged grid
file_out <- "figures/Population Dynamics Standard Scenario.tiff"
image_title <- "Standard Scenario: Population Dynamics over 7 years"
plotScenarioDynamics(scen_df_list, image_title, file_out)

# Combine dataframes into one and add temperature, year, population metabolism- and ingestion columns
TestSD <- createLongDataFrame(scen_df_list)

# Simulate scenario and get final dataframes for both types of masses
LeafSD <- simulateLeafDynamics(TestSD, "SD")
GammSD <- simulateGammarusDynamics(TestSD, "SD")

############################################################################################################# #

