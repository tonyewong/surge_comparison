##==============================================================================
## process_data.R
##
## Process data to yield time series of annual block maxima, to fit a GEV
## distribution, and a series of peaks-over-thresholds (POT), to fit a GPD,
## for a set of tide gauge data formatted as in the UHSLC repository
## (http://uhslc.soest.hawaii.edu/data/).
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.
##==============================================================================



##==============================================================================
## read in the raw tide gauge information, for the sites that will be the focus
## of this work

filetype='csv'
dat.dir <- '../input_data/tidegaugedata/'
files.tg <- list.files(path=dat.dir, pattern=filetype)
##==============================================================================



##==============================================================================
## construct GEV annual block maxima data object

data_gev <- vector('list', length(files.tg))
for (dd in 1:length(files.tg)) {
  # print an update of progress to the screen
  print(paste('now reading in data set ',dd,' / ',length(files.tg),sep=''))
  dotloc <- regexpr(".csv", files.tg[dd])
  uscloc <- regexpr("_", files.tg[dd])
  names(data_gev)[dd] <- substr(files.tg[dd], start=uscloc+1, stop=dotloc-1)
  data_gev[[dd]] <- process_gev(files.tg[dd], dat.dir)
}
##==============================================================================



##==============================================================================
## write GEV calibration data RDS file

## only take the stations with at least min_years years of available data once processed
idx_drop <- c()
for (dd in 1:length(files.tg)) {
  num_years <- nrow(data_gev[[dd]])
  if (num_years < min_years) idx_drop <- c(idx_drop, dd)
}
for (ii in rev(idx_drop)) {data_gev[[ii]] <- NULL}

## remove the stations that we won't use to avoid wasting time processing
## the GPD data set for them
files.tg <- files.tg[-idx_drop]

today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.gev <- paste('../input_data/processeddata_gev_',today,'.rds', sep='')
saveRDS(data_gev, file=filename.gev)
##==============================================================================



##==============================================================================
## construct GPD peaks-over-thresholds data object

data_gpd <- vector('list', length(files.tg))
for (dd in 1:length(files.tg)) {
  # print an update of progress to the screen
  print(paste('now reading in data set ',dd,' / ',length(files.tg),sep=''))
  dotloc <- regexpr(".csv", files.tg[dd])
  uscloc <- regexpr("_", files.tg[dd])
  names(data_gpd)[dd] <- substr(files.tg[dd], start=uscloc+1, stop=dotloc-1)
  data_gpd[[dd]] <- process_gpd(files.tg[dd], dat.dir)
  saveRDS(data_gpd, file="process_gpd_inprogress.rds")
}
##==============================================================================



##==============================================================================
## write GPD calibration data RDS file

today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.gpd <- paste('../input_data/processeddata_gpd_',today,'.rds', sep='')
saveRDS(data_gpd, file=filename.gpd)
##==============================================================================



##==============================================================================
## read in the potential covariates

source("get_timeseries_covariates.R")

# save the covariate information
filename.covariates <- paste('../input_data/covariates_',today,'.rds', sep='')
saveRDS(covariates_proj, file=filename.covariates)
##==============================================================================



##==============================================================================
## End
##==============================================================================
