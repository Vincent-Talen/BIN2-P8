## ---------------------------
##
## Script name: model.R
##
## Purpose of script: Implementing the biological model for consumer-resource dynamics
##
## Author: Vincent Talen
##
## Date Created: 10 jun 2022
##
## Email: v.k.talen@st.hanze.nl
##
## ---------------------------
##
## Notes:
##   ?
##
## ---------------------------

## ---- loading in data ----
input_file <- "Data_Mismatch.txt"
# Dataframe with respiration rate and ingestion rate (?gC/hour)
data <- read.table(input_file, header = TRUE, dec = ".") 
data <- na.omit(data)

# Considerate individuals as a qualitative variable
data$Indiv <- as.factor(data$Indiv)

# Suppress rows with negative values
data <- data[data$Respi > 0,]
data <- data[data$Nutri > 0,]

# Suppress outlier individuals
data <- data[!data$Indiv == "12",]
data <- data[!data$Indiv == "30",]
data <- data[!data$Indiv == "76",]
data <- data[!data$Indiv == "78",]

# Convert temperature into the Boltzmann term (?K)
Boltz <- 8.62 * 10^(-5)
InvTem <- 1 / ((data$Temp + 273.15) * Boltz)
InvTemC <- (InvTem - mean(InvTem))


## ---- mte formulations ----
# Values for metabolic and ingestion rates (?gC/day)
Temp <- seq(0, 30, 0.1)            # [?UNUSED?]
MeanInvTem <- mean(InvTem)         # Mean inverse temperature
MeanMass <- mean(data$Masse)       # Mean Gammarus body mass

# With quadratic function
MetaQuadra <- function(Temp, Mass) {
  exp(2.41599) * 
  Mass^0.62308 * 
  exp(-0.66731 * (1 / ((Temp + 273.15) * Boltz) - MeanInvTem)) * 
  exp(-0.21153 * (1 / ((Temp + 273.15) * Boltz) - MeanInvTem)^2)
}
IngQuadra <- function(Temp, Mass) {
  exp(5.26814) * 
  Mass^0.81654 * 
  exp(-0.31876 * (1 / ((Temp + 273.15) * Boltz) - MeanInvTem)) *
  exp(-0.18909 * (1 / ((Temp + 273.15) * Boltz) - MeanInvTem)^2)
}


## ---- assimilation efficiency ----
# Assimilation efficiency function based on exponential decay (quadratic model)
AssimQuadra <- function(Temp) {
  mte.equation <- (exp(-0.84730) * exp(0.16400 * ((Temp + 273.15) - 285.65) / (Boltz * 285.65 * (Temp + 273.15))))
  return(mte.equation / (1 + mte.equation))
}

# Calculate assimilation efficiency 
data$Assim <- AssimQuadra(data$Temp)


## ---- energetic efficiency ----
# Calculate energetic efficiency
data$Energ <- (data$Nutri / data$Respi) * data$Assim
data <- data[data$Energ < 30,]


## ---- attack rate parameter ----
Duration <- 2                # Duration of the feeding experiment (days)
MeanMass <- mean(data$Masse) # Mean Gammarus body mass (mg)
LeafMass <- 10.25 * 1000     # Mean initial leaf discs mass: 10.25 ? 0.68 mg
LeafMass <- LeafMass * 0.45  # Converted from mg to in ?gC with a factor: 0.45

# [UNUSED] Calculation of K-value and attack rate (quadratic model)
K <- -log(1 - IngQuadra(Temp, MeanMass) * 2 / LeafMass)
Attack <- K / 2

# Attack rate function based on exponential decay (quadratic model)
AttackQuadra <- function(Temp, Mass) { 
  K <- -log(1 - IngQuadra(Temp, Mass) * 2 / LeafMass)
  Attack <- K / 2
  return(Attack)
}


## ---- handling time parameter ----
# Handling time function based on exponential decay (quadratic model)
HandlingQuadra <- function(Temp, Mass) {1 / (IngQuadra(Temp, Mass) / 1000) }


## ---- [UNUSED/TESTING] energetic efficiency parameter ----
EnergQuadra <- function(Temp, MeanMass) {
  LeafMass <- 10.25 * 0.45 * 1000                   # Leaf litter mass of the feeding (in ?gC)
  Metabo <- MetaQuadra(Temp, MeanMass)              # Gammarus metabolic rate (in ?gC/day)
  a <- AttackQuadra(Temp, MeanMass)                 # Gammarus attack rate (in ?gC/day)
  h <- HandlingQuadra(Temp, MeanMass)               # Gammarus handling time (in 1/day)
  A <- AssimQuadra(Temp)                            # Gammarus assimilation efficiency
  Inges <- (a * LeafMass) / (1 + a * h * LeafMass)  # Gammarus consumption rate (in ?gC/day)
  Energ <- A * Inges / Metabo                       # Gammarus energetic efficiency
  return(Energ)
}


## ---- leaf decomposition- and respiration rate ----
# Reference temperature
TRef=10

# Function for the leaf litter microbial decomposition rate
DecompLeaf <- function(Temp) { 0.00956 * exp(-0.37000 * (1 / ((Temp + 273.15) * Boltz) - TRef)) } 

# Function for the leaf litter respiration rate
RespLeaf <- function(Temp) { 2.5 / 0.45 * 10^-3 * exp(-0.65000 * (1 / ((Temp + 273.15) * Boltz) - TRef)) } 


## ---- consumer-resource model ----
GammLeafModel <- function(Temp, GammMass, Leaf, Gamm) {
  Nutri <- function(Times, state, parameters) {
    with(as.list(c(state, parameters)), {
      Fr <- a * L / (1 + a * h * L)               # Holling type II functional response 
      dL <- -G * Fr - Dl * L                      # Biomass changes of leaf litter stock
      dG <- G * (A * Fr - M)                      # Biomass changes of Gammarus population 
      list(c(dL, dG))
    })
  }
  
  # Define the years   
  y1 <- 365; y2 <- 2 * y1; y3 <- 3 * y1; y4 <- 4 * y1; y5 <- 5 * y1; y6 <- 6 * y1
  
  # Time points to trigger litter fall
  FallTime <- c(seq(1, 15), 
                seq(y1 + 1, y1 + 15), 
                seq(y2 + 1, y2 + 15), 
                seq(y3 + 1, y3 + 15), 
                seq(y4 + 1, y4 + 15), 
                seq(y5 + 1, y5 + 15), 
                seq(y6 + 1, y6 + 15))
  
  # Leaf litter fall event function
  FallFunc <- function(time, y, parms) {
    with(as.list(c(y, parms)), {
      return(c(L + Leaf, G))
    })
  }
  
  # Model parameters
  parameters <- c(
    M = MetaQuadra(Temp, GammMass) / 1000,         # Gammarus metabolic rate (in mgC/day)
    a = AttackQuadra(Temp, GammMass),              # Gammarus attack rate (in mgC/day)
    h = HandlingQuadra(Temp, GammMass),            # Gammarus handling time (in 1/day)
    A = AssimQuadra(Temp),                         # Gammarus assimilation efficiency
    Dl = DecompLeaf(Temp),                         # Leaf microbial decomposition (in 1/day)
    Rl = RespLeaf(Temp)                            # Leaf respiration (in mgC/mgleaf/day)
  )
  
  # Time and starting conditions
  Times <- seq(0, 365 * 7, by = 1)        # Time in days for 7 years
  State <- c(L = Leaf, G = Gamm)          # Starting biomasses (in g/m2)
  
  
  # Model output
  out <- as.data.frame(ode(time = Times, func = Nutri, y = State, parms = parameters,
                           events = list(func = FallFunc, time = FallTime)))
  return(out)
}


## ---- temperature-size rule models ----
Cf=6.5                                 # Average conversion factor from dry to fresh mass

# Average TSR response
TSRA <- function(Temp, Mass) {
  Slope <- -3.90 - 0.53 * log10(Mass)     # Slope of change in mass per carbon          
  Dx <- log(1 + Slope / 100)              # Proportion of change in mass per C
  Const <- exp(log(Mass) - 12.5 * Dx)     # Constant of change in mass at 12.5?C
  DryMass <- Const * exp(Dx * (Temp))     # Dry body mass (mg)
  FreshMass <- DryMass / Cf               # Fresh body mass (mg)
  return(DryMass)
}

# Maximum TSR response
TSRM <- function(Temp, Mass) {
  Slope <- -8.0                           # Slope of change in mass per carbon          
  Dx <- log(1 + Slope/100)                # Proportion of change in mass per C
  Const <- exp(log(Mass) - 12.5 * Dx)     # Constant of change in mass at 12.5?C
  DryMass <- Const * exp(Dx * (Temp))     # Dry body mass (mg)
  FreshMass <- DryMass / Cf               # Fresh body mass (mg)
  return(DryMass)
}



