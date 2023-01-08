## Copyright (c) 2023 Vincent Talen.
## Licensed under GPLv3. See LICENSE file.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: simulateScenario.R
##
## Purpose of script: Make simulating a scenario easy with this one function
##
## Author: Vincent Talen
##
## Date Created: 08 Jan 2023
##
## Email: v.k.talen@st.hanze.nl
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Notes:
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ######### #
#   Libs    #
# ######### #
library(quantmod)
library(lme4)


# ######### #
#   Code    #
# ######### #
simulateScenario <- function(scenario_data, scenario_name, col_name, threshold) {
  # Exclude the first year and the first 16 days with leaf fall for all years for each temperature
  cut_df <- scenario_data[!Year == "1", tail(.SD, -16), by = .(Year)] %>%
    setcolorder(colnames(scenario_data)) # Restore column order
  
  # Biomass mean for each temperature ####
  Biom <- cut_df %>%
    # Calculate mean annual biomass per temperature, per year
    "["(j = .(all_means = mean(get(col_name))), by = .(Temp, Year)) %>%
    # Calculate biomass mean and deviation over 6 years per temperature
    "["(j = .(MeanBiom = mean(all_means), SdBiom = sd(all_means)), by = Temp)
  
  # Calculate persistence time for each temperature with the threshold ####
  PersTime <- cut_df %>%
    # Only keep days above the threshold
    "["(i = .[[col_name]] > threshold) %>%
    # Count the days above the threshold per temperature, per year
    "["(j = .(days_above = .N), by = .(Temp, Year)) %>%
    # Calculate persistence time mean and deviation over 6 years per temperature
    "["(j = .(MeanPersTime = mean(days_above), SdPersTime = sd(days_above)), by = Temp)
  
  # Define the biomass cycles ####
  ## Find the maximums and minimums and then get all the cycle's times ----
  CycleXSD2 <- scenario_data[by = .(Temp), j = .(Max = findPeaks(get(col_name)),#-1,
                                                 # Set minimum whilst selecting the correct correction for L or G using a switch
                                                 #Min = switch(col_name, "L" = findValleys(L)[seq(2,14,2)]-1, "G" = c(findValleys(G)-1, 2555)))] %>%
                                                 Min = switch(col_name, "L" = findValleys(L)[seq(2,14,2)], "G" = c(findValleys(G), 2555)))] %>%
    # Create new column 'Indices' with sequences of all the times in the cycles
    "$<-"(Indices, apply(., 1, function(cur_row) seq(cur_row[[2]], cur_row[[3]])))
  
  # For each temperature, get the list with indices and cycle identifiers
  createPerTempLists <- function(ind_lists) {
    # Returns a named list containing INDICES and IDENTIFIERS, both in a single array, for the given temperature
    getIdentifiers <- function(ind_lists) {
      # For each cycle create a list repeating the identifying letter for the length of that cycle
      sapply(1:length(ind_lists), function(i) rep( LETTERS[i], length(ind_lists[[i]]) ))
    }
    return( list(Indices = unlist(ind_lists), Identifiers = unlist(getIdentifiers(ind_lists))) )
  }
  per_temp_lists <- tapply(CycleXSD2$Indices, CycleXSD2$Temp, createPerTempLists)
  
  ## Combine indices and identifiers of all temperatures a single vector ----
  all_indices <- sapply(1:length(per_temp_lists), function(i) per_temp_lists[[i]]$Indices + 2555 * (i-1))
  all_identifiers <- sapply(1:length(per_temp_lists), function(i) per_temp_lists[[i]]$Identifiers)
  
  ## Subset data using the previously created vector ----
  cut_dt <- scenario_data[unlist(all_indices)] %>%
    # Add a column from the vector with identifiers
    "$<-"( Cycle, unlist(all_identifiers) ) %>%
    # Drop the first cycle for each temperature
    "["(Cycle != "A")
  
  # Calculate litter and Gammarus mean maximum biomasses and decreasing slopes ----
  ## Get the first time each cycle per temperature where the biomass is below the threshold
  first_below_threshold <- cut_dt[get(col_name) < threshold, head(.SD, 1), by = .(Temp, Cycle)]
  ## Also get the highest biomqss (first) value for each cycle per temperature
  maximum_values <- cut_dt[, head(.SD, 1), by = .(Temp, Cycle)]
  
  ## Calculate the slope for each cycle using the maximum and threshold times
  SlopLSD3 <- rbind(first_below_threshold, maximum_values) %>%
    "["( j = .(Slope = coef(lm(get(col_name) ~ Time))[[2]]), by = .(Cycle, Temp) )
  
  Maxi <- maximum_values[, .(MeanMax = mean(get(col_name)), SdMax = sd(get(col_name))), by = Temp]
  Slop <- SlopLSD3[, .(MeanSlope = mean(Slope), SdSlope = sd(Slope)), by = Temp]
  
  # Put all the information in a single data.table ----
  final_df <- setDT(c(Biom, PersTime, Maxi, Slop)) %>%
    # Remove duplicate 'Temp' columns
    "["(j = which(duplicated(names(.))) := NULL) %>%
    # Add column with scenario name/identifier
    "$<-"(Scenario, rep(scenario_name, 5))
  return(final_df)
}
