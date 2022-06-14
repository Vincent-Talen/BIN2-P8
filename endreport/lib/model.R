## ---------------------------
##
## Script name: model.R
##
## Purpose of script: Implementing the biological model for consumer-resource dynamics
##
## Author: Vincent Talen
##
## Date Created: 14 June 2022
##
## Email: v.k.talen@st.hanze.nl
##
## ---------------------------
##
## Notes:
##   - Still needs to be improved later to get 
##     an even more clear and continuous flow working towards the model
##   - Because of the new naming scheme when using supplied code keep in mind to change variable and function names
##
## ---------------------------


#############
#   Libs    #
#############
library(deSolve)
library(dplyr)


#############
#   Code    #
#############
## ---- loading in data ----
# Dataframe with respiration rate and ingestion rate (μg C/hour)
input_file <- "Data_Mismatch.txt"
data <- read.table(input_file, header = TRUE, dec = ".") 
data <- na.omit(data)

# Consider individuals as a qualitative variable
data$Indiv <- as.factor(data$Indiv)

# Suppress outlier individuals and rows with negative values
data <- data[!data$Indiv %in% c("12", "30", "76", "78"), ]
data <- data[data$Respi > 0 & data$Nutri > 0, ]


## ---- mte formulations ----
# Convert temperature into the Boltzmann term (°K)
boltz_const <- 8.62 * 10^-5
inversed_temps <- 1 / ((data$Temp + 273.15) * boltz_const)
mean_inverse_temp <- mean(inversed_temps)

# Quadratic functions for metabolic and ingestion rates (μg C/day)
calcMetabolicRate <- function(temp, mass) {
  exp(2.41599) * 
  mass^0.62308 * 
  exp(-0.66731 * (1 / ((temp + 273.15) * boltz_const) - mean_inverse_temp)) * 
  exp(-0.21153 * (1 / ((temp + 273.15) * boltz_const) - mean_inverse_temp)^2)
}
calcIngestionRate <- function(temp, mass) {
  exp(5.26814) * 
  mass^0.81654 * 
  exp(-0.31876 * (1 / ((temp + 273.15) * boltz_const) - mean_inverse_temp)) *
  exp(-0.18909 * (1 / ((temp + 273.15) * boltz_const) - mean_inverse_temp)^2)
}


## ---- assimilation efficiency ----
# Assimilation efficiency function based on exponential decay (quadratic model)
# Is a logistic equation with the MTE equation both at the numerator and the denominator
calcAssimEff <- function(temp) {
  mte_equation <- (exp(-0.84730) * exp(0.16400 * ((temp + 273.15) - 285.65) / (boltz_const * 285.65 * (temp + 273.15))))
  return(mte_equation / (1 + mte_equation))
}


## ---- attack rate parameter ----
leaf_mass <- 10.25 * 1000      # Mean initial leaf discs mass: 10.25 ± 0.68 mg
leaf_mass <- leaf_mass * 0.45  # Converted from mg to in μgC with a factor: 0.45

# Attack rate function based on exponential decay (quadratic model)
calcAttackRate <- function(temp, mass) {
  expt_duration <- 2
  decay_rate <- -log(1 - calcIngestionRate(temp, mass) * expt_duration / leaf_mass)
  attack_rate <- decay_rate / expt_duration
  return(attack_rate)
}


## ---- handling time parameter ----
# Handling time function based on exponential decay (quadratic model)
calcHandlingTime <- function(temp, mass) {1 / (calcIngestionRate(temp, mass) / 1000) }

## ---- leaf decomposition- and respiration rate ----
# Reference temperature
ref_temp <- 10

# Function for the leaf litter microbial decomposition rate
calcLeafDecomp <- function(temp) { 0.00956 * exp(-0.37000 * (1 / ((temp + 273.15) * boltz_const) - ref_temp)) } 

# Function for the leaf litter respiration rate
calcLeafRespi <- function(temp) { 2.5 / 0.45 * 10^-3 * exp(-0.65000 * (1 / ((temp + 273.15) * boltz_const) - ref_temp)) } 


## ---- consumer-resource model ----
GammLeafModel <- function(temp, gamm_indv_mass, leaf_biomass, gamm_pop_biomass) {
  Nutri <- function(time, state, parms) {
    with(as.list(c(state, parms)), {
      fL <- a * L / (1 + a * h * L)         # Holling type II functional response 
      dL <- -fL * G - k * L                 # Biomass changes of leaf litter stock
      dG <- G * (fL * A - M)                # Biomass changes of Gammarus population 
      list(c(dL, dG))
    })
  }
  
  # Leaf litter fall event function
  litterFallEvent <- function(time, state, parms) {
    with(as.list(c(state, parms)), {
      return(c(L + leaf_biomass, G))
    })
  }
  
  # Get time points to trigger litter fall event (first 15 days of the year)
  getFallTimesYearX <- function(year) {
    start_day_year <- year * 365
    return(seq(start_day_year + 1, start_day_year + 15))
  }
  litter_fall_times <- unlist(lapply(seq(0, 6), getFallTimesYearX))
  
  # Model parameters
  parameters <- c(
    M = calcMetabolicRate(temp, gamm_indv_mass) / 1000,  # Gammarus metabolic rate (in mgC/day)
    a = calcAttackRate(temp, gamm_indv_mass),            # Gammarus attack rate (in mgC/day)
    h = calcHandlingTime(temp, gamm_indv_mass),          # Gammarus handling time (in 1/day)
    A = calcAssimEff(temp),                              # Gammarus assimilation efficiency
    k = calcLeafDecomp(temp),                            # Leaf microbial decomposition (in 1/day)
    Rl = calcLeafRespi(temp)                             # [UNUSED???] Leaf respiration (in mgC/mgleaf/day)
  )
  
  # Times and starting conditions
  times <- seq(0, 365 * 7, by = 1)                      # Times in days for 7 years
  state <- c(L = leaf_biomass, G = gamm_pop_biomass)    # Starting biomasses (in g/m2)
  
  # Model output
  out <- as.data.frame(ode(time = times, func = Nutri, y = state, parms = parameters,
                           events = list(func = litterFallEvent, time = litter_fall_times)))
  return(out)
}


## ---- temperature-size rule models ----
# Average TSR response
calcTSR.Avg <- function(temp, mass) {
  conv_fact <- 6.5                                      # Avg. conversion factor from dry to fresh mass
  change_slope <- -3.90 - 0.53 * log10(mass)            # Slope of change in mass per carbon          
  change_prop <- log(1 + change_slope / 100)            # Proportion of change in mass per C
  change_const <- exp(log(mass) - 12.5 * change_prop)   # Constant of change in mass at 12.5°C
  
  dry_mass <- change_const * exp(change_prop * (temp))  # Dry body mass (mg)
  fresh_mass <- dry_mass / conv_fact                    # Fresh body mass (mg)
  return(dry_mass)
}

# Maximum TSR response
calcTSR.Max <- function(temp, mass) {
  conv_fact <- 6.5                                      # Avg. conversion factor from dry to fresh mass
  change_slope <- -8.0                                  # Slope of change in mass per carbon          
  change_prop <- log(1 + change_slope / 100)            # Proportion of change in mass per C
  change_const <- exp(log(mass) - 12.5 * change_prop)   # Constant of change in mass at 12.5°C
  
  dry_mass <- change_const * exp(change_prop * (temp))  # Dry body mass (mg)
  fresh_mass <- dry_mass / conv_fact                    # Fresh body mass (mg)
  return(dry_mass)
}


