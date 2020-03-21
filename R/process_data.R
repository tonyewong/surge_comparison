##==============================================================================
## process_data.R
##
## Process data to yield time series of annual block maxima, to fit a GEV
## distribution, for a set of tide gauge data formatted as in the UHSLC
## repository (http://uhslc.soest.hawaii.edu/data/).
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================


library(date)


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
}
##==============================================================================



##==============================================================================
## write GEV calibration data RDS file

today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.gev <- paste('../input_data/processeddata_gev_',today,'.rds', sep='')
saveRDS(data_gev, file=filename.gev)
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
