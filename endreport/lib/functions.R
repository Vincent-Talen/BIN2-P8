## ---------------------------
##
## Script name: functions.R
##
## Purpose of script: Functions
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
library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
source("lib/model.R")


#############
# Functions #
#############
createPlotForTemperature <- function(temp, gamm_indv_mass, leaf_fall, gamm_start_biomass) {
  # Call model function that performs ode with current values
  data.simulated <- GammLeafModel(temp, gamm_indv_mass, leaf_fall, gamm_start_biomass)
  
  # Divide L & G values to create a better readable plot
  divided.data <- data.simulated %>% 
    mutate(across(c(L, G), ~ .x / 10^5))
  
  # Create plot
  plot <- ggplot(divided.data, aes(x = time, y = value)) + 
    # Set axis limits and step size
    scale_x_continuous(breaks = seq(0, 7 * 365, 365)) +
    ylim(NA, ceiling(max(divided.data[, -1])) + 1) +
    # Add the data (lines)
    geom_line(aes(y = L, color = "Leaf Litter Biomass")) + 
    geom_line(aes(y = G, color = "Gammarus Fossarum Biomass")) +
    # Add styling
    labs(title = sprintf("%s°C", temp), x = "", y = "") +
    # theme(plot.title = element_text(hjust = 0.075, vjust = -11)) +
    #theme(plot.margin = margin(0.1, 0.25, 0, 0, "cm")) +
    scale_color_manual(name = "", values = c("black","tomato2"), 
                       limits = c("Leaf Litter Biomass", "Gammarus Fossarum Biomass"))
  return(plot)
}


simulateScenario <- function(file_out, image_title, gamm_indv_mass, leaf_fall, gamm_start_biomass) {
  # Use lapply to create a list with a plot for each temperature and collect the legend from a plot
  plot_list <- lapply(c(5, 10, 15, 20, 25), createPlotForTemperature, gamm_indv_mass, leaf_fall, gamm_start_biomass)
  plot_legend <- get_legend(plot_list[[1]])
  
  # Remove legends from the plots and add extracted legend to end of the list
  plot_list <- lapply(plot_list, function(cur_plot) {cur_plot + theme(legend.position = "none")}) %>%
    "[[<-"(length(plot_list) + 1, plot_legend)
  
  # Place plots and legend in an arranged grid 
  col_num <- 3
  my.grid <- ggarrange(plotlist = plot_list, ncol = col_num, nrow = ceiling(length(plot_list) / col_num)) %>%
    annotate_figure(top = text_grob(image_title), bottom = text_grob("Time (d)"), 
                    left = text_grob(bquote("Biomass "(10^5~ mg~ C~ m^-2)), rot = 90))
  
  # Save the created arranged grid
  # The lossless 'lzw' compression GREATLY reduces file size BUT generating, opening and loading files takes significantly longer
  ggsave(file_out, bg = "white", width=15, height=8, units="in", dpi=1000)#, compression = "lzw")
  dev.off()
}