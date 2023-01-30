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

## Plotting biomasses #####################################################
lossePlot <- function(deze_data, deze_pred, col_names, plot_title) {
  # y-axis limits
  max_lim <- round(max(deze_data[,get(col_names["Mean"])])/10^4+1)
  break_by <- ifelse(max_lim <= 4, 1, 2)

  plot <- ggplot(deze_data, aes(x=Temperature, y=!!sym(col_names["Mean"])/10^4, group=Scenario)) +
    # Plot data
    geom_point(aes(color=Scenario), size=3, position=position_dodge(0.5)) +
    geom_line(data=deze_pred, aes(x=Temperature, y=!!sym(col_names["Pred"])/10^4, color=Scenario, linetype=Scenario), show.legend=F) +
    geom_errorbar(color="grey50", width=0.75, position=position_dodge(0.5),
                  aes(ymin=(!!sym(col_names["Mean"]) - !!sym(col_names["Sd"]))/10^4, 
                      ymax=(!!sym(col_names["Mean"]) + !!sym(col_names["Sd"]))/10^4)) +
    # Axis (text) styling
    labs(x = "Temperature (\u00B0C)", y = expression('Mean biomass'~'('*10^4~mg~C~m^-2*')'), title = plot_title) +
    theme(axis.text.y=element_text(size=10, colour="black"), axis.text.x=element_text(size=10, colour="black"),
          plot.title = element_text(size=12), legend.title=element_text(face="bold")) +
    scale_y_continuous(breaks=seq(0, max_lim, by=break_by), limits=c(0, max_lim)) +#, labels=function(x) sprintf("%.0f", x)) +
    # Color and types of plotted data elements
    scale_color_manual(values=c(SD="black", TSRA="steelblue2", TSRM="tomato2")) +
    scale_linetype_manual(values=c("SD"="dotted", "TSRA"="solid", "TSRM"="solid"))
  return(plot)
}

plotGridOfColumn <- function(deze_data, deze_pred, col_name) {
  # Create plots and put them in a list
  plot_list <- lapply(c("L", "G"), function(l_or_g) {
    plot_title <- ifelse(l_or_g == "L", "Leaf Litter", "Gammarus")
    lossePlot(deze_data, deze_pred, sapply(c("Mean", "Sd", "Pred"), paste, l_or_g, sep=col_name), plot_title)
  })
  
  # Return the plots in an arranged grid with the legend at the bottom
  arranged_plots <- ggpubr::ggarrange(plotlist=plot_list, ncol=2, nrow=1, common.legend=TRUE, legend="bottom")
  return(arranged_plots)
}


createPlotImageForColumn <- function(column_name, image_title, out_file) {
  plotGridOfColumn(Data, Pred, column_name) %>%
           annotate_figure(top = text_grob(image_title, size=16, face="bold"))
  ggsave(paste("figures/", out_file, sep=""), bg="white", width=3840, height=2160, units="px", dpi=300, compression="lzw")
  dev.off()
}

createPlotImageForColumn("Biom", "Mean Biomass Scenarios", "Mean Biomass Scenarios.tiff")

