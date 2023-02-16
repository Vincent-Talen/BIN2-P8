## Copyright (c) 2023 Vincent Talen.
## Licensed under GPLv3. See LICENSE file.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: model.R
##
## Purpose of script: Implements the biological model of consumer-resource dynamics
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
##   - Goal: formulate formula functions in a more readable/recognizable manner
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ######### #
#   Libs    #
# ######### #
library(deSolve)


# ######### #
#   Code    #
# ######### #
# Define Boltzmann term (°K)
boltz_const <- 8.62 * 10^-5
# Define mean inverse temperature (calculated from Data_Mismatch.txt)
mean_inverse_temp <- 40.5941593143742


## ---- mte formulations ----
# Quadratic function for metabolic rate (μg C/day)
calcMetabolicRate <- function(T.C, M) {
  alpha <- exp(2.41599)       # metabolic expression level at reference temperature
  b <- 0.62308                # body mass-scaling exponent
  p <- 0.66731                # curve steepness (of the relationship)
  q <- 0.21153                # quadratic term
  T.K <- T.C + 273.15         # convert temperature from Celsius to Kelvin
  
  # Repeating part with temperatures
  temp_dependancy_part <- (1 / (T.K * boltz_const)) - mean_inverse_temp
  
  # Calculate metabolic rate with full formula
  metabolic_rate <- alpha * M^b * exp(-p * temp_dependancy_part) * exp(-q * temp_dependancy_part^2)
  return(metabolic_rate)
}

# Quadratic function for ingestion rate (μg C/day)
calcIngestionRate <- function(T.C, M) {
  alpha <- exp(5.26814)       # ingestion expression level at reference temperature
  b <- 0.81654                # body mass-scaling exponent
  p <- 0.31876                # curve steepness (of the relationship)
  q <- 0.18909                # quadratic term
  T.K <- T.C + 273.15         # convert temperature from Celsius to Kelvin
  
  # Repeating part with temperatures
  temp_dependancy_part <- (1 / (T.K * boltz_const)) - mean_inverse_temp
  
  # Calculate ingestion rate with full formula
  ingestion_rate <- alpha * M^b * exp(-p * temp_dependancy_part) * exp(-q * temp_dependancy_part^2)
  return(ingestion_rate)
}


## ---- assimilation efficiency ----
# Assimilation efficiency function based on exponential decay (quadratic model)
# Is a logistic equation with the MTE equation both at the numerator and the denominator
calcAssimEff <- function(T.C) {
  alpha <- exp(-0.84730)    # normalization constant of assimilation efficiency
  Ea <- 0.16400             # activation energy
  T.0 <- 285.65             # reference temperature of 12.5 degrees Celsius in Kelvin
  T.K <- T.C + 273.15       # convert temperature from Celsius to Kelvin
  
  mte_equation <- alpha * exp( Ea * (T.K - T.0) / (boltz_const * T.0 * T.K) )
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
GammLeafModel <- function(temp, gamm_indv_mass, leaf_fall, gamm_start_biomass, tsr_model) {
  # Apply TSR Model to Gammarus individual mass
  if (!is.null(tsr_model)) {
    gamm_indv_mass <- tsr_model(temp, gamm_indv_mass)
  }
  
  Nutri <- function(time, state, parms) {
    with(as.list(c(state, parms)), {
      fL <- a * L / (1 + a * h * L)         # Holling type II functional response 
      dL <- -fL * G - k * L                 # Biomass changes of leaf litter stock
      dG <- G * (fL * A - M)                # Biomass changes of Gammarus population 
      list(c(dL, dG))
    })
  }
  
  # Leaf litter fall event function
  leafFallEvent <- function(time, state, parms) {
    with(as.list(c(state, parms)), {
      return(c(L + leaf_fall, G))
    })
  }
  
  # Get time points to trigger litter fall event (first 15 days of the year)
  getFallTimesYearX <- function(year) { seq(year * 365 + 1, year * 365 + 15) }
  leaf_fall_times <- unlist(lapply(seq(0, 6), getFallTimesYearX))
  
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
  times <- seq(0, 365 * 7, by = 1)                       # Times in days for 7 years
  state <- c(L = leaf_fall, G = gamm_start_biomass)      # Starting biomasses (in g/m2)
  
  # Model output
  out <- ode(time = times, func = Nutri, y = state, parms = parameters,
             events = list(func = leafFallEvent, time = leaf_fall_times))
  
  # Turn deSolve class object into dataframe and change very low and negative values to 0
  data_table <- as.data.table(out) %>% mutate(across(c(L, G), ~ fifelse(.x < 10^-3, 0, .x)))
  return(data_table)
}


## ---- temperature-size rule models ----
# Average TSR response
calcTSR.Avg <- function(temp, mass) {
  conv_fact <- 6.5                                       # Avg. conversion factor from dry to fresh mass
  change_slope <- -3.90 - 0.53 * log10(mass)             # Slope of change in mass per carbon          
  change_prop <- log(1 + change_slope / 100)             # Proportion of change in mass per C
  change_const <- exp(log(mass) - 12.5 * change_prop)    # Constant of change in mass at 12.5°C
  
  dry_mass <- change_const * exp(change_prop * (temp))   # Dry body mass (mg)
  fresh_mass <- dry_mass / conv_fact                     # Fresh body mass (mg)
  return(dry_mass)
}

# Maximum TSR response
calcTSR.Max <- function(temp, mass) {
  conv_fact <- 6.5                                       # Avg. conversion factor from dry to fresh mass
  change_slope <- -8.0                                   # Slope of change in mass per carbon          
  change_prop <- log(1 + change_slope / 100)             # Proportion of change in mass per C
  change_const <- exp(log(mass) - 12.5 * change_prop)    # Constant of change in mass at 12.5°C
  
  dry_mass <- change_const * exp(change_prop * (temp))   # Dry body mass (mg)
  fresh_mass <- dry_mass / conv_fact                     # Fresh body mass (mg)
  return(dry_mass)
}
