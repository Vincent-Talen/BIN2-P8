## ---------------------------
##
## Script name: scenarios.R
##
## Purpose of script: Simulate scenarios using functions from the functions.R module
##
## Author: Vincent Talen
##
## Date Created: 17 June 2022
##
## Email: v.k.talen@st.hanze.nl
##
## ---------------------------
##
## Notes:
##   - x
##
## ---------------------------


#############
#   Libs    #
#############
source("lib/functions.R")


#############
#   Code    #
#############
#####################################
### SCENARIO 0: STANDARD SCENARIO ###
#####################################
# Gammarus mean body mass = 4.26 mgDM
# Annual leaf fall        = 300 gC/m2/an = 300 000 mgC/m2/an
# Gammarus density        = 30 mgDM/m2   = 15 mgC/m2

# Duration of the leaves fall
Days <- 15

# Scenario values
gamm_indv_mass <- 4.26
leaf_fall <- 300000 / Days
gamm_start_biomass <- 15

# Perform simulation of scenario
file_out <- "output/Population Dynamics Standard Scenario.tiff"
image_title <- "Standard Scenario: Population Dynamics over 7 years"
simulateScenario(file_out, image_title, gamm_indv_mass, leaf_fall, gamm_start_biomass)
