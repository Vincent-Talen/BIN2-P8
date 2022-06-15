## ---------------------------
##
## Script name: scenarios.R
##
## Purpose of script: Simulate scenarios using the model made in model.R script
##
## Author: Vincent Talen
##
## Date Created: 15 June 2022
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
library(dplyr)
library(ggplot2)
source("lib/model.R")
source("lib/functions.R")


#############
#   Code    #
#############
#####################################
### SCENARIO 0: STANDARD SCENARIO ###
#####################################
# Leaf fall = 300 gC/m2/an = 300 000 mgC/m2/an
# Gammarus density = 30 mgDM/m2 = 15 mgC/m2
# Gammarus mean body mass = 4.26 mgDM

## ---- data los ff ----
simulation.data <- GammLeafModel(temp = 5, gamm_indv_mass = 4.26, daily_leaf_in = 300000 / Days, start_gamm_biomass = 15) %>% 
  mutate(L = if_else(L < 10^-3, 0, L)) %>% 
  mutate(G = if_else(G < 10^-3, 0, G))


## ---- functietjes ----
mine.perTemperature <- function(temp, gamm_indv_mass, daily_leaf_in, start_gamm_biomass) {
  data.simulated <- GammLeafModel(temp, gamm_indv_mass, daily_leaf_in, start_gamm_biomass) %>% 
    mutate(L = if_else(L < 10^-3, 0, L)) %>% 
    mutate(G = if_else(G < 10^-3, 0, G))
  #rownames(data.simulated) <- data.simulated[ ,1]
  #data.simulated[ ,1] <- NULL
  
  ## Create plot for current state variable
  plot <- ggplot(data.simulated, aes(x = time)) +
    labs(title = as.character(temp), x = "", y = "") +
    ### Add the lines and data points
    geom_line(aes(color = "Leaf Litter Biomass"), data = data.simulated["L"]) +
    geom_line(aes(color = "Gammarus Fossarum Biomass"), data = data.simulated["G"]) +
    theme(legend.position = "bottom", plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm")) + 
    ### Use scale_color_manual to color the data and also set their names in the legend
    scale_color_manual(name = "", values = c("black","red"), limits = c("Leaf Litter Biomass", "Gammarus Fossarum Biomass")) +
    ### Correct the displayed line types in the legend
    #guides(color = guide_legend(override.aes = list(linetype = c(1, 1))))
  return(plot)
}
mine.mainFunction()

them.perTemperature <- function(temp, gamm_indv_mass, daily_leaf_in, start_gamm_biomass) {
  # hoe zij hebben
  Test5SD <- GammLeafModel(temp, gamm_indv_mass, daily_leaf_in, start_gamm_biomass)
  Test5SD <- Test5SD %>% mutate(L = if_else(L < 10^-3, 0, L))
  Test5SD <- Test5SD %>% mutate(G = if_else(G < 10^-3, 0, G))
  
  matplot(Test5SD[, -1] / 10^5, type = "l", ylim = c(0, 6), ylab = "", xlab = "", cex.axis = 2.0, las = 1, col = c("black", "tomato2"), lty = c(1, 1, 1), lwd = c(1.5, 2))
  mtext(sprintf("%s°C", temp), cex = 1.5, side = 3, line = -3, at = 300)
}
them.mainFunction()

## ---- uitvoeren ----
# Duration of the leaves fall
Days <- 15

mine.mainFunction <- function() {
  # Scenario values
  gamm_indv_mass <- 4.26
  daily_leaf_in <- 300000 / Days
  start_gamm_biomass <- 15
  
  ## Create plots
  plots <- lapply(c(5, 10, 15, 20, 25), mine.perTemperature, gamm_indv_mass, daily_leaf_in, start_gamm_biomass)
  ## Print the plots in an arranged grid with the legend at the bottom
  printAndArrangePlots(plot.list = plots, grid.title = "[Standard Scenario]: Temporal Dynamics over 7 years")
}

them.mainFunction <- function() {
  # Plot the dynamics for 5 temperatures
  tiff('output/Population Dynamics SD.tiff', units = "in", width = 15, height = 8, res = 1000)
  par(mfrow = c(2, 3))
  par(oma = c(4, 4, 1, 1))
  par(mar = c(3, 3, 0, 0))
  
  them.perTemperature(temp = 5,  gamm_indv_mass = 4.26, daily_leaf_in = 300000 / Days, start_gamm_biomass = 15)
  them.perTemperature(temp = 10, gamm_indv_mass = 4.26, daily_leaf_in = 300000 / Days, start_gamm_biomass = 15)
  them.perTemperature(temp = 15, gamm_indv_mass = 4.26, daily_leaf_in = 300000 / Days, start_gamm_biomass = 15)
  them.perTemperature(temp = 20, gamm_indv_mass = 4.26, daily_leaf_in = 300000 / Days, start_gamm_biomass = 15)
  them.perTemperature(temp = 25, gamm_indv_mass = 4.26, daily_leaf_in = 300000 / Days, start_gamm_biomass = 15)
  
  mtext("Time (d)", side = 1, outer = TRUE, line = 1, cex = 1.5)
  mtext(expression('Biomass'~'('*10^5~mg~C~m^-2*')'), side = 2, outer = TRUE, line = 1, cex = 1.5)
  par(mfrow=c(1, 1))
  dev.off()
}

## ---- idk ----
# Main function for assignment 1: If given a valid dose, plots will be made for it's data
a1.mainFunction <- function(current.dose) {
  ## Change D parameter to median of current dose
  cur.median.MPL_conc <- median(data$MPL_conc[data$dose == current.dose], na.rm = T)
  basic.parameters$D <- calc.D(cur.median.MPL_conc)
  
  ## Get subset of experiment data of current dose
  cur.data.points <- subset(data, dose == 0.0 | dose == current.dose)
  ## Get subset of experiment median data of current dose
  cur.data.medians <- subset(data.medians, dose == 0.0 | dose == current.dose)
  ## Perform ODE function with the model to get simulation data
  cur.data.simulated <- as.data.frame(ode(func = modelBasic, times = basic.times, 
                                          y = basic.zero.state, parms = basic.parameters))
  
  ## Create plots
  plots <- lapply(c("Rm", "R"), a1.createPlot, 
                  cur.data.simulated, cur.data.medians, cur.data.points)
  ## Print the plots in an arranged grid with the legend at the bottom
  printAndArrangePlots(plot.list = plots,
                       grid.title = sprintf("Graphs for dose = %s (mg drug/kg rat/h)", current.dose))
}

# Function that fully creates a plot for given state variable with given data
a1.createPlot <- function(state.var, data.simulated, data.medians, data.points) {
  ## Collect current state variable from data frame that contains name, unit and title
  cur.var <- subset(state.info.df, name == state.var)
  
  ## Create plot for current state variable
  plot <- ggplot(data.simulated, aes(x = time, y = !!sym(cur.var[["name"]]))) +
    labs(title = cur.var["title"], x = "time (h)", 
         y = sprintf("%s (%s)", cur.var["name"], cur.var["unit"])) +
    ### Add the lines and data points
    geom_line(aes(color = "Simulation Data")) +
    geom_line(aes(color = "Experiment Median"), data = data.medians) + 
    geom_point(aes(color = "Experiment Data"), data = data.points, shape = 1) + 
    theme(legend.position = "bottom", plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm")) + 
    ### Use scale_color_manual to color the data and also set their names in the legend
    scale_color_manual(name = "", 
                       values = c("black","red","black"),
                       limits = c("Simulation Data", "Experiment Median", "Experiment Data")) +
    ### Correct the displayed line types in the legend
    guides(color = guide_legend(
      override.aes = list(linetype = c(1, 1, NA), shape = c(NA, NA, 1))))
  return(plot)
}