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


# ######### #
#   Libs    #
# ######### #
source("src/functions.R")
source("src/simulateScenario.R")


# ######### #
#   Code    #
# ######### #
# SETTINGS FOR ALL SCENARIOS ##############################################
# Temperatures to do simulations of
temperatures <- c(5, 10, 15, 20, 25)

# Duration of the leaf fall in days
fall_duration_in_days <- 15

# Gammarus mean body mass = 4.26 mgDM
gamm_indv_mass <- 4.26
# Annual leaf fall = 300 gC/m2/an = 300 000 mgC/m2/an
leaf_fall <- 300000 / fall_duration_in_days
# Gammarus density = 30 mgDM/m2 = 15 mgC/m2
gamm_start_biomass <- 15


# SCENARIO 0: STANDARD SCENARIO ###########################################
# Get data for temperatures with values of current scenario
scen_df_list <- getScenarioDataList(gamm_indv_mass, leaf_fall, gamm_start_biomass, NULL)

# Create plots in an arranged grid
file_out <- "Population Dynamics Standard Scenario.tiff"
image_title <- "Standard Scenario: Population Dynamics over 7 years"
plotScenarioDynamics(scen_df_list, image_title, file_out)

# Combine dataframes into one and add temperature, year, population metabolism- and ingestion columns
TestSD <- createLongDataFrame(scen_df_list, NULL)
# Simulate scenario and get final dataframes for both types of masses
DataSD <- simulateScenario(TestSD, "SD")


# SCENARIO 1: AVERAGE TSR RESPONSE ########################################
# Get data for temperatures with values of current scenario
scen_df_list2 <- getScenarioDataList(gamm_indv_mass, leaf_fall, gamm_start_biomass, calcTSR.Avg)

# Create plots in an arranged grid
file_out <- "Population Dynamics Average TSR.tiff"
image_title <- "Average Temperature-Size Rule Response: Population Dynamics over 7 years"
plotScenarioDynamics(scen_df_list2, image_title, file_out)

# Combine dataframes into one and add temperature, year, population metabolism- and ingestion columns
TestTSRA <- createLongDataFrame(scen_df_list2, calcTSR.Avg)
# Simulate scenario and get final dataframes for both types of masses
DataTSRA <- simulateScenario(TestTSRA, "TSRA")


# SCENARIO 2: MAXIMUM TSR RESPONSE ########################################
# Get data for temperatures with values of current scenario
scen_df_list3 <- getScenarioDataList(gamm_indv_mass, leaf_fall, gamm_start_biomass, calcTSR.Max)

# Create plots in an arranged grid
file_out <- "Population Dynamics Maximum TSR.tiff"
image_title <- "Maximum Temperature-Size Rule Response: Population Dynamics over 7 years"
plotScenarioDynamics(scen_df_list3, image_title, file_out)

# Combine dataframes into one and add temperature, year, population metabolism- and ingestion columns
TestTSRM <- createLongDataFrame(scen_df_list3, calcTSR.Max)
# Simulate scenario and get final dataframes for both types of masses
DataTSRM <- simulateScenario(TestTSRM, "TSRM")

# COMBINING SCENARIO OUTPUTS ##############################################
Data <- rbind(DataSD, DataTSRA, DataTSRM)
Pred <- createPredictionData(Data, c("SD", "TSRA", "TSRM"))
