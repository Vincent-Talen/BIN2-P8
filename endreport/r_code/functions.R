## ---------------------------
##
## Script name: functions.R
##
## Purpose of script: Functions
##
## Author: Vincent Talen
##
## Date Created: 20 June 2022
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
library(data.table)
library(ggpubr)
library(reshape2)
library(tidyverse)
source("r_code/model.R")


#############
# Functions #
#############
# ---- scenario data gathering and preparations ----
getScenarioDataList <- function(gamm_indv_mass, leaf_fall, gamm_start_biomass) {
  # Get data for given values for each temperature using the model function that performs an ode
  data_list <- lapply(temperatures, GammLeafModel, gamm_indv_mass, leaf_fall, gamm_start_biomass) %>% 
    setNames(temperatures)
  return(data_list)
}

prepareBigDataFrame <- function(df_list) {
  # Function to get the population metabolism for a temperature with the population biomass
  getPopMetabolism <- function(cur_temp, gamm_pop_biomass) {
    # Get metabolic rate for current temperature
    meta_rate   <- calcMetabolicRate(cur_temp, gamm_indv_mass) / 1000  # Gammarus metabolic rate (in mg C/day)
    # Calculate population metabolism
    pop_metabolism <- meta_rate * gamm_pop_biomass
    return(fifelse(pop_metabolism < 0, 0, pop_metabolism))
  }
  # Function to get the population ingestion for a temperature with the population- and leaf biomasses
  getPopIngestion <- function(cur_temp, leaf_biomass, gamm_pop_biomass) {
    # Get ingestion- and attack rates for current temperature
    ingest_rate <- calcIngestionRate(cur_temp, gamm_indv_mass) / 1000  # Gammarus ingestion rate (in mg C/day)
    attack_rate <- calcAttackRate(cur_temp, gamm_indv_mass)            # Gammarus attack rate (in mg C/day)
    # Calculate population leaf ingestion
    pop_ingestion <- (attack_rate * leaf_biomass / (1 + attack_rate * 1 / (ingest_rate / 1000) * leaf_biomass)) * 0.30 * gamm_pop_biomass
    return(fifelse(pop_ingestion < 0, 0, pop_ingestion))
  }
  
  # Bind all dataframes from list to single big one and 
  # drop last days to have 2555 days/rows left per temperature
  big_df <- rbindlist(df_list)[!time == 2555] %>%
    # Rename 'time' column to conform to naming scheme
    setnames("time", "Time") %>%
    # Add temperature and year columns to facilitate future calculations
    "$<-"(Temp, rep(temperatures, each = 2555)) %>%
    "$<-"(Year, rep(rep(1:7, each = 365), 5)) %>%
    # Add population metabolism and leaf ingestion columns
    "$<-"(M, getPopMetabolism(.$Temp, .$G)) %>%
    "$<-"(I, getPopIngestion(.$Temp, .$L, .$G)) %>%
    # Set column order to a nicer one
    setcolorder(c("Time", "L", "G", "M", "I", "Temp", "Year"))
  return(big_df)
}

# ---- plotting list of scenario dataframes ----
createPlotForTemp <- function(cur_temp, cur_data) {
  # Divide L & G values to create a better readable plot
  divided_data <- copy(cur_data)
  set(divided_data, i = NULL, "L", divided_data$L / 10^5)
  set(divided_data, i = NULL, "G", divided_data$G / 10^5)
  
  # Create plot
  plot <- ggplot(divided_data, aes(x = time, y = value)) + 
    # Set axis limits and step size
    scale_x_continuous(breaks = seq(0, 7 * 365, 365)) +
    ylim(NA, ceiling(max(divided_data[, -1])) + 1) +
    # Add the data (lines)
    geom_line(aes(y = L, color = "Leaf Litter Biomass")) + 
    geom_line(aes(y = G, color = "Gammarus Fossarum Biomass")) +
    # Add styling
    labs(title = sprintf("%s°C", cur_temp), x = "", y = "") +
    # theme(plot.title = element_text(hjust = 0.075, vjust = -11)) +
    #theme(plot.margin = margin(0.1, 0.25, 0, 0, "cm")) +
    scale_color_manual(name = "", values = c("black","tomato2"), 
                       limits = c("Leaf Litter Biomass", "Gammarus Fossarum Biomass"))
  return(plot)
}

plotScenario <- function(data, image_title, file_out) {
  # Use lapply to create plots for each temperature in the list and collect the legend from a plot
  plot_list <- lapply(seq_along(data), function(i) { createPlotForTemp(names(data)[i], data[[i]]) })
  plot_legend <- get_legend(plot_list[[1]])
  
  # Remove legends from the plots and add extracted legend to end of the list
  plot_list <- lapply(plot_list, function(cur_plot) { cur_plot + theme(legend.position = "none") }) %>%
    "[[<-"(length(plot_list) + 1, plot_legend)
  
  # Place plots and legend in an arranged grid 
  col_num <- 3
  my.grid <- ggarrange(plotlist = plot_list, ncol = col_num, nrow = ceiling(length(plot_list) / col_num)) %>%
    annotate_figure(top = text_grob(image_title), bottom = text_grob("Time (d)"), 
                    left = text_grob(bquote("Biomass "(10^5~ mg~ C~ m^-2)), rot = 90))
  
  # Save the created arranged grid with the lossless 'lzw' compression that greatly reduces file size
  ggsave(file_out, bg = "white", width=15, height=8, units="in", dpi=300, compression = "lzw")
  dev.off()
}
