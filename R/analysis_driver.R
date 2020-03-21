##==============================================================================
## analysis_driver.R
##
## To be run after calibration_driver.R
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

rm(list=ls())

library(extRemes)
library(Bolstad)

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/surge_comparison/R')
} else if(Sys.info()['user']=='aewsma') {
  machine <- 'office'
  setwd('/Users/aewsma/codes/surge_comparison/R')
} else {
  # assume on another cluster of some kind...
  machine <- 'remote'
  setwd('~/work/codes/surge_comparison/R')
}

source("best_models.R") # yields parameters[[model]][site,parameter]
source("make_projections.R") # yields return level estimates


##==============================================================================
## End
##==============================================================================
