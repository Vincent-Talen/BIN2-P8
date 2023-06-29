## Copyright (c) 2023 Vincent Talen.
## Licensed under GPLv3. See LICENSE file.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: mainAnalysis.R
##
## Purpose of script: Perform all scenario simulations for the main analysis
##
## Author: Vincent Talen
##
## Date Created: 29 Jun 2023
##
## Email: v.k.talen@st.hanze.nl
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Notes:
##   - x
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~


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


# SCENARIO 0: REFERENCE SCENARIO ##########################################
# Get data for temperatures with values of current scenario
df_list.TSRR <- getScenarioDataList(gamm_indv_mass, leaf_fall, gamm_start_biomass, NULL)

# Create plots in an arranged grid
file_out.TSRR <- "Population Dynamics Reference Scenario.png"
image_title.TSRR <- "Reference Scenario: Population Dynamics over 7 years"
plotScenarioDynamics(df_list.TSRR, image_title.TSRR, file_out.TSRR)

# Combine dataframes into one and add temperature, year, population metabolism- and ingestion columns
combined_df.TSRR <- createLongDataFrame(df_list.TSRR, NULL)
# Simulate scenario and get final dataframes for both types of masses
data.TSRR <- simulateScenario(combined_df.TSRR, "TSRR")


# SCENARIO 1: AVERAGE TSR RESPONSE ########################################
# Get data for temperatures with values of current scenario
df_list.TSRA <- getScenarioDataList(gamm_indv_mass, leaf_fall, gamm_start_biomass, calcTSR.Avg)

# Create plots in an arranged grid
file_out.TSRA <- "Population Dynamics Average TSR Scenario.png"
image_title.TSRA <- "Average Temperature-Size Rule Response: Population Dynamics over 7 years"
plotScenarioDynamics(df_list.TSRA, image_title.TSRA, file_out.TSRA)

# Combine dataframes into one and add temperature, year, population metabolism- and ingestion columns
combined_df.TSRA <- createLongDataFrame(df_list.TSRA, calcTSR.Avg)
# Simulate scenario and get final dataframes for both types of masses
data.TSRA <- simulateScenario(combined_df.TSRA, "TSRA")


# SCENARIO 2: MAXIMUM TSR RESPONSE ########################################
# Get data for temperatures with values of current scenario
df_list.TSRM <- getScenarioDataList(gamm_indv_mass, leaf_fall, gamm_start_biomass, calcTSR.Max)

# Create plots in an arranged grid
file_out.TSRM <- "Population Dynamics Maximum TSR Scenario.png"
image_title.TSRM <- "Maximum Temperature-Size Rule Response: Population Dynamics over 7 years"
plotScenarioDynamics(df_list.TSRM, image_title.TSRM, file_out.TSRM)

# Combine dataframes into one and add temperature, year, population metabolism- and ingestion columns
combined_df.TSRM <- createLongDataFrame(df_list.TSRM, calcTSR.Max)
# Simulate scenario and get final dataframes for both types of masses
data.TSRM <- simulateScenario(combined_df.TSRM, "TSRM")


## PLOTTING FUNCTIONS #####################################################
createPlotForStateVariable <- function(state_variable_name, statistic_info, data_column_names, all_scenario_data, prediction_data) {
  # Calculate the y-axis limits
  max_lim <- round( max(all_scenario_data[,get(data_column_names["Mean"])]) + statistic_info$correction )
  break_by <- ifelse(max_lim <= 4, 1, 2)
  
  plot <- ggplot(all_scenario_data, aes(x=Temperature, y=abs(!!sym(data_column_names["Mean"])), group=Scenario)) +
    # Add the actual scenario simulation points to the plot
    geom_point(
      aes(color=Scenario), 
      size=3, 
      position=position_dodge(0.5)
    ) +
    # Add the prediction data lines to the plot
    geom_line(
      data=prediction_data, 
      aes(x=Temperature, y=abs(!!sym(data_column_names["Pred"])), color=Scenario, linetype=Scenario), 
      show.legend=F
    ) +
    # Add the standard deviation error bar to the plot
    geom_errorbar(
      color="grey50", 
      width=0.75, 
      position=position_dodge(0.5),
      aes(ymin=(abs(!!sym(data_column_names["Mean"]) - !!sym(data_column_names["Sd"]))), 
          ymax=(abs(!!sym(data_column_names["Mean"]) + !!sym(data_column_names["Sd"]))))
    ) +
    # Add the plot title and axis labels
    labs(x = "Temperature (\u00B0C)", y=statistic_info$y_label, title=state_variable_name) +
    # Style the plot title, axis and legend texts
    theme(
      axis.text.y=element_text(size=10, colour="black"), 
      axis.text.x=element_text(size=10, colour="black"),
      plot.title = element_text(size=12), 
      legend.title=element_text(face="bold")
    ) +
    # Set the y-axis limits
    scale_y_continuous(labels=function(x) {x / statistic_info$correction}) + 
    expand_limits(y = 0) +
    # Color and types of plotted data elements
    scale_color_manual(values=c(TSRR="black", TSRA="steelblue2", TSRM="tomato2")) +
    scale_linetype_manual(values=c("TSRR"="dotted", "TSRA"="solid", "TSRM"="solid"))
  return(plot)
}

plotGridOfStatistic <- function(statistic_info, all_scenario_data, prediction_data) {
  # Define state variables
  state_variable_list <- list(
    list(id="L", full_name="Leaf Litter"), 
    list(id="G", full_name="Gammarus Fossarum")
  )
  # Create the plots for the state variables in a list using lapply
  state_variable_plots <- lapply(state_variable_list, function(state_variable) {
    # Get the full names of the data columns for the current statistic + state variable combination
    data_column_names <- sapply(
      c("Mean", "Sd", "Pred"),  # for each data column
      paste,                    # perform the pasta function
      state_variable$id,        # put the id of the state variable at the end
      sep=statistic_info$name   # and the name of the statistic in between
    )
    # Create the plot of the data columns for the current statistic + state variable combination
    createPlotForStateVariable(state_variable$full_name, statistic_info, data_column_names, all_scenario_data, prediction_data)
  })
  
  # Return the plots in an arranged grid with the legend at the bottom
  arranged_plots <- ggpubr::ggarrange(plotlist=state_variable_plots, ncol=2, nrow=1, common.legend=TRUE, legend="bottom")
  return(arranged_plots)
}

createFigureForStatistic <- function(statistic_info, all_scenario_data) {
  # Get prediction data that is used for the continuous lines in the plots
  prediction_data <- createPredictionData(all_scenario_data, c("TSRR", "TSRA", "TSRM"))
  # Create the arranged grid with the two plots of the state variables for the statistic and annotate the figure title
  plotGridOfStatistic(statistic_info, all_scenario_data, prediction_data) %>%
    annotate_figure(top=text_grob(statistic_info$image_title, size=16, face="bold"))
  # Save the final figure for the column
  ggsave(
    filename=paste("figures/reproduced_plots/", statistic_info$image_title, ".png", sep=""),
    bg="white",
    width=3840,
    height=2160,
    units="px",
    dpi=300
  )
  dev.off()
}


## CREATING PLOTS #########################################################
# Combine the data frames into a single one
combined_data <- rbind(data.TSRR, data.TSRA, data.TSRM)

# Create a list with the statistics that a figure needs to be made for
statistic_info_list <- list(
  list(name = "Biom", 
       correction = 10^4, 
       y_label = expression('Mean biomass'~'('*10^4~mg~C~m^-2*')'), 
       image_title = "Mean Biomass over Temperature"),
  list(name = "PersTime", 
       correction = 1, 
       y_label = "Persistance time (days)", 
       image_title = "Persistance Time over Temperature"),
  list(name = "Slope", 
       correction = 10^4, 
       y_label = expression('Mean biomass slope'~'('*10^4~mg~C~m^-2~day^-1*')'), 
       image_title = "Mean Biomass Slope over Temperature")
)

# Create a figures for each statistic defined above
for (statistic_info in statistic_info_list) {
  createFigureForStatistic(statistic_info, combined_data)
}
