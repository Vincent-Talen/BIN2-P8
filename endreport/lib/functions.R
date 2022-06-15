## ---------------------------
##
## Script name: functions.R
##
## Purpose of script: Functions
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
library(ggplot2)
library(ggpubr)


#############
# Functions #
#############
# Function to print list of plots in an arranged grid with common legend at the bottom
printAndArrangePlots <- function(plot.list, grid.title) {
  my.grid <- ggarrange(plotlist = plot.list, common.legend = TRUE, legend = "bottom", ncol = 2, nrow = ceiling(length(plot.list) / 2))
  print(annotate_figure(my.grid, top = text_grob(grid.title)))
}
