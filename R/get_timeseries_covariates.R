##==============================================================================
## get_timeseries_covariates.R
##
## First, everything is normalized to first 20 years (1928-1937) of the tide
## gauge data. Then, the covariates are each normalized so that the min-max
## range for the hindcast period is 0-1.
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


library(ncdf4)


# check the GEV data sets (annual block maxima) to see what years we need the
# covariates in for calibration
data_gev <- readRDS('../input_data/processeddata_gev_22Mar2020.rds')
first_year <- rep(NA, length(data_gev))
last_year <- rep(NA, length(data_gev))
for (dd in 1:length(first_year)) {
  first_year[dd] <- min(data_gev[[dd]][,1])
  last_year[dd] <- max(data_gev[[dd]][,1])
}
years <- (min(first_year)-1):max(last_year)


# get NAO index ================================================================
nao_dat <- read.table('../input_data/nao_3dp.dat')
colnames(nao_dat) <- c('year','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec','ann')
ibeg <- which(nao_dat['year']==1850)
iend <- max(which(nao_dat[,'ann']!=-99.99))
time_hist <- nao_dat[ibeg:iend, 'year']

# get DJF means
nao_hist <- rep(-999, length(time_hist))
for (y in 1:length(time_hist)) {
  nao_hist[y] <- mean( c(nao_dat$dec[ibeg+y-1], nao_dat$jan[ibeg+y], nao_dat$feb[ibeg+y]) )
}
time_hist <- nao_dat[ibeg:iend, 'year']

# re-normalize (mean and stdev) relative to 2001-2016 mean/stdev, so it is
# consistent with the projections
ind_norm <- which(time_hist==2001):which(time_hist==2016)
nao_hist <- (nao_hist - mean(nao_hist[ind_norm]))/sd(nao_hist[ind_norm])

# normalize relative to first 20 years of storm surge data
# (doing for both historical and projection, then stitching together and renormalizing)
ibeg <- which(time_hist==years[1])
iend <- which(time_hist==years[20])
nao_hist <- nao_hist - mean(nao_hist[ibeg:iend])

nao <- cbind(time_hist, nao_hist)
colnames(nao) <- c('year','nao')


# get GMST =====================================================================
data.tmp <- read.table('../input_data/noaa_temperature_1880-2017.csv', header = TRUE, sep=',')
time_hist <- data.tmp$Year
temperature_hist <- data.tmp$Value
temperature <- cbind(time_hist, temperature_hist)
colnames(temperature) <- c('year','temperature')


# get GMSL =====================================================================
sl.dat.new = read.table("../input_data/GMSL_ChurchWhite2011_yr_2015.txt")     #Reconstructed sea-level
SL.time= sl.dat.new[,1]-0.5     # times are in half-year
SL.new = sl.dat.new[,2]/1000    # data are in mm
SL.err = sl.dat.new[,3]/1000
ibeg=which(SL.time==1961); iend=which(SL.time==1990);
SL.new = SL.new - mean(SL.new[ibeg:iend])   # SL data are relative to 1961-1990
sealevel <- cbind(SL.time, SL.new)
colnames(sealevel) <- c('year','sealevel')


# just use Time (year relative to first year?) =================================
time <- cbind(years, years)
colnames(time) <- c('year','year')


# trim all potential features to ts.years
min_year <- max( min(nao[,'year']),      min(temperature[,'year']),
                 min(sealevel[,'year']), min(time[,'year'])        )
max_year <- min( max(nao[,'year']),      max(temperature[,'year']),
                 max(sealevel[,'year']), max(time[,'year'])        )

nao <- nao[which(nao[,1]==min_year):which(nao[,1]==max_year), ]
temperature <- temperature[which(temperature[,1]==min_year):which(temperature[,1]==max_year), ]
sealevel <- sealevel[which(sealevel[,1]==min_year):which(sealevel[,1]==max_year), ]
time <- time[which(time[,1]==min_year):which(time[,1]==max_year), ]


#===============================================================================
# create a single array to hold all the possible covariates

names_covariates <- c('time','temp','sealevel','nao')
covariates <- cbind(time[,2], temperature[,2], sealevel[,2], nao[,2])
colnames(covariates) <- names_covariates


#===============================================================================
#===============================================================================
# get time series for forcing/making return level projections


# CESM temperature projections =================================================
# gridcell areas
ncin <- nc_open("../input_data/cesm/GridcellAreaAtmosphere-CESM2-CMIP6_accessed26Jan2022/areacella_fx_CESM2_ssp585_r4i1p1f1_gn_v20200528.nc")
areaa <- ncvar_get(ncin,"areacella")
nc_close(ncin)
areaa_wgts <- areaa/sum(areaa)

# historical (Jan 1850 - Dec 2014)
ncin <- nc_open("../input_data/cesm/SurfaceTemperatureHistorical-CESM2-CMIP6_accessed26Jan2022/ts_Amon_CESM2_historical_r1i1p1f1_gn_185001-201412_v20190308.nc")
time_ts_hist <- ncvar_get(ncin,"time")
ts_hist <- ncvar_get(ncin,"ts")
nc_close(ncin)
# convert K to deg C
ts_hist <- ts_hist - 273.15
# global mean
ts_hist_global_monthly <- rep(-99999, length(time_ts_hist))
areaa_wgts <- areaa/sum(areaa)
for (i in 1:length(time_ts_hist)) {
  ts_hist_global_monthly[i] <- sum(ts_hist[,,i]*areaa_wgts)
}
# annual mean
ts_hist_global_annual <- rep(-99999, length(time_ts_hist)/12)
for (i in 1:(length(time_ts_hist)/12)) {
  ts_hist_global_annual[i] <- mean(ts_hist_global_monthly[((i-1)*12+1):(i*12)])
}

# 2015-2064
ncin <- nc_open("../input_data/cesm/SurfaceTemperature-CESM2-CMIP6_accessed26Jan2022/ts_Amon_CESM2_ssp585_r4i1p1f1_gn_201501-206412_v20200528.nc")
time_2015_2064 <- ncvar_get(ncin,"time")
ts_2015_2064 <- ncvar_get(ncin,"ts")
nc_close(ncin) # units: ncin$var$ts$units
# 2065-2100
ncin <- nc_open("../input_data/cesm/SurfaceTemperature-CESM2-CMIP6_accessed26Jan2022/ts_Amon_CESM2_ssp585_r4i1p1f1_gn_206501-210012_v20200528.nc")
time_2065_2100 <- ncvar_get(ncin,"time")
ts_2065_2100 <- ncvar_get(ncin,"ts")
nc_close(ncin)

# stitch projections together
time_ts <- c(time_2015_2064, time_2065_2100)
ts <- abind(ts_2015_2064, ts_2065_2100)
# convert K to deg C
ts <- ts - 273.15
# global aggregation
ts_global_monthly <- rep(-99999, length(time_ts))
for (i in 1:length(time_ts)) {
  ts_global_monthly[i] <- sum(ts[,,i]*areaa_wgts)
}
# annual means
# matches well against Figure 1a from Meehl et al 2020 (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020EA001296)
ts_global_annual <- rep(-99999, length(time_ts)/12)
for (i in 1:(length(time_ts)/12)) {
  ts_global_annual[i] <- mean(ts_global_monthly[((i-1)*12+1):(i*12)])
}

# stitch together the CESM2 historical temperatures from the repository
temperature_proj <- c(ts_hist_global_annual, ts_global_annual)
year_proj <- 1850:2100
# normalize like covariates for calibration
ibeg <- which(year_proj==years[1])
iend <- which(year_proj==years[20])
temperature_proj <- temperature_proj - mean(temperature_proj[ibeg:iend])
temperature_proj <- cbind(year_proj, temperature_proj)
colnames(temperature_proj) <- c('year','temp')


# Time =========================================================================
timec_proj <- temperature_proj[,'year'] - mean(temperature_proj[ibeg:iend,'year'])
timec_proj <- cbind(temperature_proj[,'year'], timec_proj)
colnames(timec_proj) <- c('year','time')


# BRICK sea level projection (from Wong and Keller 2017) =======================
ncdata <- nc_open("../input_data/BRICK_GMSL_WongKeller2017.nc")
gmsl_rcp85 <- t(ncvar_get(ncdata, 'GlobalSeaLevel_RCP85'))
time_proj <- ncvar_get(ncdata, 'time_proj')
nc_close(ncdata)

# get projection that yields median in 2100 under RCP8.5
ibeg <- 1
iend <- nrow(gmsl_rcp85)
# need only odd number so median = element
if (iend%%2 == 0) {iend <- iend-1}
ind_median <- which(gmsl_rcp85[,length(time_proj)]==median(gmsl_rcp85[ibeg:iend,length(time_proj)]))
gmsl_proj <- cbind(time_proj, gmsl_rcp85[ind_median,])
colnames(gmsl_proj) <- c('year','sealevel')

# normalize as in hindcast
ibeg <- which(time_proj==years[1])
iend <- which(time_proj==years[20])
gmsl_proj[,2] <- gmsl_proj[,2] - mean(gmsl_proj[ibeg:iend,2])


# NAO index ====================================================================
# (using method of Stephenson et al 2006, as in Wong et al 2018)
ncin <- nc_open("../input_data/cesm/SeaLevelPressureHistorical-CESM2-CMIP6_accessed26Jan2022/psl_Amon_CESM2_historical_r1i1p1f1_gn_185001-201412_v20190308.nc")
lon <- ncvar_get(ncin, "lon")
lat <- ncvar_get(ncin, "lat")
time_1850_2014 <- ncvar_get(ncin,"time")
psl_1850_2014 <- ncvar_get(ncin,"psl")
nc_close(ncin)

ncin <- nc_open("../input_data/cesm/SeaLevelPressure-CESM2-CMIP6_accessed26Jan2022/psl_Amon_CESM2_ssp585_r4i1p1f1_gn_201501-206412_v20200528.nc")
time_2015_2064 <- ncvar_get(ncin,"time")
psl_2015_2064 <- ncvar_get(ncin,"psl")
nc_close(ncin)

ncin <- nc_open("../input_data/cesm/SeaLevelPressure-CESM2-CMIP6_accessed26Jan2022/psl_Amon_CESM2_ssp585_r4i1p1f1_gn_206501-210012_v20200528.nc")
time_2065_2100 <- ncvar_get(ncin,"time")
psl_2065_2100 <- ncvar_get(ncin,"psl")
nc_close(ncin)

time_psl <- c(time_1850_2014,time_2015_2064, time_2065_2100)
psl <- abind(psl_1850_2014,psl_2015_2064, psl_2065_2100)

n_month <- length(time_psl)
n_year <- n_month/12 -1 # -1 reflects number of complete DJF NOA indices we'll get

# As in Stephenson et al, 2006 (doi:  10.1007/s00382-006-0140-x)
# use regional, since there is disagreement over how exactly to do the PCA
# (which EOFs to rotate, e.g.) and a physical interpretation is direct this way
# From CESM, -90 < lat < 90 and 0 < lon < 360
ilat_azores <- which(lat >= 20 & lat <= 55)
ilon_azores <- which(lon >= (360-90) | lon <= 60)

ilat_iceland <- which(lat >= 55 & lat <= 90)
ilon_iceland <- which(lon >= (360-90) | lon <= 60)

psl_azores <- psl[ilon_azores, ilat_azores, ]
psl_iceland <- psl[ilon_iceland, ilat_iceland, ]

# Don't do area-weighting as per explanation from Jesse Nusbaumer 27 March 2018
# email. Jesse says:
#   Don't area-weight the results.  Instead, just take the average of the SLP
#   over the specified region, treating every grid box the same.  The reason is
#   because in GCM papers they will usually specifically state "area-weighted"
#   if the actual meters-squared area is being taken into account, otherwise it
#   is assumed to just be a non-weighted average.  Also when it comes to a lot
#   of these indices the actual physical values (e.g. the total amount of
#   atmospheric mass producing the surface pressure) isn't really important,
#   they are just using the values to produce a strong statistical relationship.

# take the bulk average over each spatial area, for each month
psl_azores_spavg <- apply(psl_azores, 3, mean)
psl_iceland_spavg <- apply(psl_iceland, 3, mean)

# Standardize each site separately (as discussed in Jones et al 1997). Do relative
# to 2001-2016 mean/stdev so we can be consistent between projections and the
# historical record
psl_azores_std <- rep(-99999, n_month)
psl_iceland_std <- rep(-99999, n_month)

# indices corresponding to 2001-2016?
# > ncin$dim$time$units
# [1] "days since 0001-01-01 00:00:00"
# 1850-2000 (inclusive) = 12 months/year * 151 years = 1812 months
# 2001-2016 (inclusive) = 12 months/year * 16 years = 192 months
idx_2001_2016 <- (12*151+1):(12*151+12*16)

for (m in 1:12) {
  ind_this_month <- seq(from=m, to=n_month, by=12) # gather up all of month m's data points
  idx_2001_2016_this_month <- intersect(ind_this_month, idx_2001_2016)
  psl_tmp <- psl_azores_spavg[ind_this_month]      # only the Azores avgs for this month
  psl_azores_std[ind_this_month] <- (psl_tmp - mean(psl_azores_spavg[idx_2001_2016_this_month]))/sd(psl_azores_spavg[idx_2001_2016_this_month]) # slot the normalized values in
  psl_tmp <- psl_iceland_spavg[ind_this_month]
  psl_iceland_std[ind_this_month] <- (psl_tmp - mean(psl_azores_spavg[idx_2001_2016_this_month]))/sd(psl_azores_spavg[idx_2001_2016_this_month])
}

# get SLP difference
nao_monthly <- psl_azores_std - psl_iceland_std

# get winter mean
nao_proj <- rep(-99999, n_year)
for (y in 1:n_year) {
  nao_proj[y] <- mean(nao_monthly[(y-1)*12 + 12:14])  # DJF
#  nao_proj[y] <- mean(nao_monthly[(y-1)*12 + 12:15])  # DJFM
#  nao_proj[y] <- mean(nao_monthly[(y-1)*12 + 1:12])   # annual
}
time_proj <- 1850:2099 # years of projection; data point DJF avg for year Y has Dec in year Y

# normalize to first 20 years of storm surge data
ibeg <- which(time_proj==years[1])
iend <- which(time_proj==years[20])
nao_proj <- nao_proj - mean(nao_proj[ibeg:iend])

# paste together with the nao (historical)
ibeg <- which(time_proj==(max(nao[,1])+1))
iend <- length(time_proj)
nao_proj <- cbind(c(nao[,1], time_proj[ibeg:iend]),
                  c(nao[,2], nao_proj[ibeg:iend]))

# and normalize all together
ibeg <- which(nao_proj[,1]==years[1])
iend <- which(nao_proj[,1]==years[20])
nao_proj[,2] <- nao_proj[,2] - mean(nao_proj[ibeg:iend,2])
colnames(nao_proj) <- c('year','nao')


# put all together as in covariates ============================================

# trim to only the period we have all of them
min_year <- max( min(nao_proj[,'year']),  min(temperature_proj[,'year']),
                 min(gmsl_proj[,'year']), min(timec_proj[,'year'])        )
max_year <- min( max(nao_proj[,'year']),  max(temperature_proj[,'year']),
                 max(gmsl_proj[,'year']), max(timec_proj[,'year'])        )

year_proj <- min_year:max_year

##
## checking normalization of covariate time series
##

ibeg <- which(timec_proj[,'year']==min_year)
iend <- which(timec_proj[,'year']==max_year)
timec_proj <- timec_proj[ibeg:iend,]
ibeg <- which(nao_proj[,'year']==min_year)
iend <- which(nao_proj[,'year']==max_year)
nao_proj <- nao_proj[ibeg:iend,]
ibeg <- which(temperature_proj[,'year']==min_year)
iend <- which(temperature_proj[,'year']==max_year)
temperature_proj <- temperature_proj[ibeg:iend,]
ibeg <- which(gmsl_proj[,'year']==min_year)
iend <- which(gmsl_proj[,'year']==max_year)
gmsl_proj <- gmsl_proj[ibeg:iend,]

covariates_proj <- cbind(year_proj, timec_proj[,2], temperature_proj[,2], gmsl_proj[,2], nao_proj[,2])
colnames(covariates_proj) <- c('year',names_covariates) # time temp sealevel nao


# first, normalize the covariates also relative to the first 20 years
# of available tide gauge data
years_norm <- covariates[1:20,"time"]
idx_norm_hind <- match(years_norm, covariates[,"time"])
idx_norm_proj <- match(years_norm, covariates_proj[,"year"])
for (cc in names_covariates) {
  covariates[,cc] <- covariates[,cc] - mean(covariates[idx_norm_hind, cc])
  covariates_proj[,cc] <- covariates_proj[,cc] - mean(covariates_proj[idx_norm_proj, cc])
}

# now paste in the historical data for the period it's available for all time
# series
ibeg <- which(covariates[,"time"]==covariates_proj[1,"time"])
iend <- nrow(covariates)
covariates_proj[1:(iend-ibeg+1), names_covariates] <- covariates[ibeg:iend, names_covariates]

# finally, normalize so the hindcast period (covariates) are on 0-1
norms <- array(NA, c(length(names_covariates),2))
rownames(norms) <- names_covariates
colnames(norms) <- c('min','max')

for (cc in names_covariates) {
  # then, between 0 and 1
  norms[cc,'min'] <- min(covariates[ibeg:iend,cc])
  norms[cc,'max'] <- max(covariates[ibeg:iend,cc])
  covariates[,cc] <- (covariates[,cc] - norms[cc,'min'])/(norms[cc,'max'] - norms[cc,'min'])
  covariates_proj[,cc] <- (covariates_proj[,cc] - norms[cc,'min'])/(norms[cc,'max'] - norms[cc,'min'])
}

covariates <- cbind(covariates_proj[1:nrow(covariates), "year"] ,covariates)
colnames(covariates)[1] <- "year"

##==============================================================================
## End
##==============================================================================
