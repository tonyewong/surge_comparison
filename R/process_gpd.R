##==============================================================================
## process_gpd.R
##
## Function to take in a tide gauge hourly data set and return a time series of
## peaks-over-threshold local sea levels, together with the corresponding year.
## Yields output object that will show count of 0 exceedances during years with
## missing data, but the time_length will be less than a year, so it will behave
## nicely in the likelihood_gpd.R functions.
##
## filename_in = string, file name
## data_dir = string, directory path to the file nae, relative to current directory
## threshold_missing_data = filter out days that are missing more data than this (%)
## gpd_threshold = percentile (0-1) to serve as the GPD threshold (only modeling exceedances of this threshold)
## dt_decluster = declustering time-scale (days)
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

process_gpd <- function(filename_in, data_dir, threshold_missing_data=0.9, gpd_threshold=0.99, dt_decluster=3) {

  tbeg_total <- proc.time()

  fillvalue <- -32767 # fill-value
  data.tg <- read.csv(paste(data_dir, filename_in, sep=''))
  names(data.tg) <- c('year','month','day','hour','sl')
  data.tg$time.days <- as.numeric(mdy.date(month=data.tg$month, day=data.tg$day,
                                  year=data.tg$year)) + data.tg$hour/24

  # difference between time stamps (in units of days)
  time.diff <- diff(data.tg$time.days)

  # check that they are in the proper order, ascending
  print(paste('Are there any times out of order? ',any(time.diff < 0), sep=''))

  # put everything back in order - make sure you do this for the sea levels and
  # other fields too, so do as a matrix. also recalculate the time differences,
  # which you will need for averaging
# Note -  data from PMSLC is all in order
#  data.tg <- data.tg[order(data.tg$time.days),]
#  time.diff <- diff(data.tg$time.days)

  # what is one hour? in units of days
  one.hours <- 1/24

  # where are there missing data?
  idx_missing <- which(data.tg$sl == fillvalue)
  data.tg[idx_missing,"sl"] <- NA

  # Detrend by either subtracting annual means (moving 1-year window)
  print(paste('Detrending by subtracting moving window 1-year average ...', sep=''))

  # get a placeholder
  data.tg$sl.detrended <- data.tg$sl
  time.days.beg <- min(data.tg$time.days)
  time.days.end <- max(data.tg$time.days)

# for gaps:
# 1) detrend entire data set with linear interpolation
# 2) fit a mean annual cycle to the detrended hourly data
# 3) add the linear interpolation back in
# 4) any gaps are filled in with this annual cycle + linear trend

# -> 1) detrend
idx_present <- which(!is.na(data.tg$sl))
fit <- lm(data.tg$sl[idx_present] ~ data.tg$time.days[idx_present])
tg_trend <- fit$coefficients[1] + fit$coefficients[2]*data.tg$time.days
tg_detrended <- data.tg$sl - tg_trend

# -> 2) fit mean annual cycle
tg_cycle <- rep(NA, n_idx_per_year)
pb <- txtProgressBar(min=1,max=n_idx_per_year,initial=0,style=3)
for (tt in 1:n_idx_per_year) {
    idx_this_hour <- which( (data.tg$hour==data.tg$hour[tt]) &
                            (data.tg$day==data.tg$day[tt]) &
                            (data.tg$month==data.tg$month[tt]) )
    tg_cycle[tt] <- mean(tg_detrended[idx_this_hour], na.rm=TRUE)
    setTxtProgressBar(pb, tt)
}
close(pb)

# -> 3) add linear interpolation back in
# repeat annual cycle to match the full data length
tg_cycle_full <- rep(tg_cycle, ceiling(length(data.tg$sl)/n_idx_per_year))
tg_cycle_full <- tg_cycle_full[1:length(data.tg$sl)]
tg_cycle_full <- tg_cycle_full + tg_trend

# -> 4) fill gaps with this annual cycle + linear trend
data.tg$sl[idx_missing] <- tg_cycle_full[idx_missing]

# detrending the first half-year and last half-year
idx_first_year <- which(data.tg$time.days - time.days.beg <= 365.25)
idx_last_year <- which(time.days.end - data.tg$time.days <= 365.25)
n_idx_per_year <- 365.25*24+1 # +1 for the current step, and half year in each direction

idx_first_halfyear <- which(data.tg$time.days - time.days.beg < (365.25*0.5))
idx_last_halfyear <- which(time.days.end - data.tg$time.days < (365.25*0.5))

data.tg$sl.detrended[idx_first_halfyear] <- data.tg$sl[idx_first_halfyear] - mean(data.tg$sl[idx_first_year])
data.tg$sl.detrended[idx_last_halfyear] <- data.tg$sl[idx_last_halfyear] - mean(data.tg$sl[idx_last_year])

# beginning and ending indices for the middle bit
ibeg <- max(idx_first_halfyear) + 1
iend <- min(idx_last_halfyear) - 1

tbeg <- proc.time()
sl <- data.tg$sl
sl.detrended <- data.tg$sl.detrended
idx_window <- (ibeg - floor(n_idx_per_year/2)):(ibeg + floor(n_idx_per_year/2))
pb <- txtProgressBar(min=ibeg,max=iend,initial=0,style=3)
for (tt in ibeg:iend) {
    #data.tg$sl.detrended[tt] <- data.tg$sl[tt] - mean(data.tg$sl[idx_window], na.rm=TRUE)
    sl.detrended[tt] <- sl[tt] - mean(sl[idx_window], na.rm=TRUE)
    idx_window <- idx_window + 1
    setTxtProgressBar(pb, tt)
}
close(pb)
data.tg$sl.detrended[ibeg:iend] <- sl.detrended[ibeg:iend] # and data.tg$sl was unchanged
tend <- proc.time()
print(paste("Total detrending time:",(tend-tbeg)[3]/60,"minutes"))
# original took 93.27 minutes out of 106.92 minutes total
# replacing the data.tg with placeholders then assigning later took 6.30 minutes
# Still expect about 20 minutes/station, times 36 stations --> about 12 hours to process

# after this, the points with missing data have values - remove these
data.tg$sl.detrended[idx_missing] <- NA
data.tg$sl[idx_missing] <- NA


  # daily block maxima; calculate 99% quantile as GPD threshold

  # how many days in each year have at least 90% of their values?
  days.all <- floor(data.tg$time.days)
  days.unique <- unique(days.all)
  ind.days.to.remove <- NULL
  print('... filtering down to do a daily maxima time series of only the days with at least 90% of data ...')

  # first day
  ind.today <- which(floor(data.tg$time.days) == days.unique[1])
  perc.today <- (length(ind.today) - length(intersect(ind.today, idx_missing)))/24
  if (perc.today < 0.9) {ind.days.to.remove <- c(ind.days.to.remove, 1)}

  # last day
  ind.today <- which(floor(data.tg$time.days) == max(days.unique))
  perc.today <- (length(ind.today) - length(intersect(ind.today, idx_missing)))/24
  if (perc.today < 0.9) {ind.days.to.remove <- c(ind.days.to.remove, length(days.unique))}

  # in-between days - take advantage of continuity of data set
  #  and idx_missing are the only ones we need to worry about
  # --> only check unique(floor(data.tg$time.days[idx_missing]))

  print('... removing days that are missing too many data points ...')

  # days that are missing any data:
  days.missing <- unique(floor(data.tg$time.days[idx_missing]))

  tbeg <- proc.time()
  pb <- txtProgressBar(min=1,max=length(days.missing),initial=0,style=3)
  for (tt in 1:length(days.missing)) {
    ind.today <- which(floor(data.tg$time.days) == days.missing[tt])
    perc.missing.today <- length(intersect(ind.today, idx_missing))/24
    if (perc.missing.today > 0.1) {ind.days.to.remove <- c(ind.days.to.remove, match(days.missing[tt], days.unique))}
    setTxtProgressBar(pb, tt)
  }
  close(pb)
  tend <- proc.time()
  print(paste("Total missing data removal time:",(tend-tbeg)[3]/60,"minutes"))
  days.daily.max <- days.unique[-ind.days.to.remove]
  n.days <- length(days.daily.max)

  print('... getting time series of daily maxima ...')

  # calculate the daily maximum sea levels on the days of 'days.daily.max'
  sl.daily.max <- rep(NA, n.days)
  years.daily.max <- rep(NA, n.days)
  tbeg <- proc.time()
  pb <- txtProgressBar(min=0,max=n.days,initial=0,style=3)
  for (day in days.daily.max) {
    cnt <- match(day,days.daily.max)
    ind.today <- which(days.all == day)
    sl.daily.max[cnt] <- max(data.tg$sl.detrended[ind.today])
    years.daily.max[cnt] <- data.tg$year[ind.today][1]
    setTxtProgressBar(pb, cnt)
  }
  close(pb)
  tend <- proc.time()
  print(paste("Total daily maxima calculation time:",(tend-tbeg)[3]/60,"minutes"))

  # find all the excesses, "declustering" = if two are within dt.decluster of
  # each other, take only the maximum of the two (so make sure you save the
  # times of each excess)

  print('... getting threshold excesses ...')

  tbeg <- proc.time()
  # Make a list object for output for this site with everything we need
  # to calibrate the PP-GPD model
  data_out <- vector('list', 7)
  names(data_out) <- c('counts','year','time_length','excesses','threshold','p.threshold','dt.decluster')

  # threshold is 99% quantile of tide gauge's observed values.
  # Buchanan et al (2016) use the 99% quantile of the daily maximum time series.
  gpd.threshold <- as.numeric(quantile(sl.daily.max, gpd_threshold, na.rm=TRUE))
  data_out$dt.decluster <- dt_decluster
  data_out$p.threshold <- gpd_threshold
  data_out$threshold <- gpd.threshold

  ind.exceed <- which(sl.daily.max > gpd.threshold)
  days.exceed <- days.daily.max[ind.exceed]
  sl.exceed <- sl.daily.max[ind.exceed]
  years.exceed <- years.daily.max[ind.exceed]

  declustered.exceed <- decluster_timeseries(time=days.exceed, year=years.exceed, time.series=sl.exceed, min.dt=dt_decluster)
  days.exceed.decl <- declustered.exceed$time
  years.exceed.decl <- declustered.exceed$year
  sl.exceed.decl <- declustered.exceed$time.series

  # initialize
  years.unique <- unique(data.tg$year)
  data_out$counts <- data_out$year <- data_out$time_length <- rep(NA, length(years.unique))
  data_out$excesses <- vector('list', length(years.unique))

  for (ind.year in 1:length(years.unique)) {
    ind.hits.this.year <- which(years.exceed.decl == years.unique[ind.year])
    data_out$counts[ind.year] <- length(ind.hits.this.year)
    data_out$year[ind.year]   <- years.unique[ind.year]
    data_out$time_length[ind.year] <- length(which(years.daily.max == years.unique[ind.year]))
    if(length(ind.hits.this.year) > 0) {data_out$excesses[[ind.year]] <- sl.exceed.decl[ind.hits.this.year]
    } else                             {data_out$excesses[[ind.year]] <- NA}
  }
  tend <- proc.time()
  print(paste("Total threshold exceedances calculation time:",(tend-tbeg)[3]/60,"minutes"))

  # alternatively, could bin em all together. but this won't allow for potential
  # non-stationary behavior in the poisson process
  data_out$excesses_all <- sl.exceed.decl
  data_out$counts_all <- length(sl.exceed.decl)
  data_out$time_length_all <- length(days.daily.max)

  # that takes some time, so save the workspace image after each data set
#  save.image(file=filename.saveprogress)

  tend_total <- proc.time()
  print(paste('  ... done. Took ', (tend_total[3]-tbeg_total[3])/60, ' minutes.',sep=''))

  return(data_out)
}

##==============================================================================
## End
##==============================================================================
