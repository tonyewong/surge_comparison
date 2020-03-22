##==============================================================================
## process_gpd.R
##


print("TODO!")


TODO   TODO   TODO   TODO   TODO   TODO
TODO   TODO   TODO   TODO   TODO   TODO
TODO   TODO   TODO   TODO   TODO   TODO


## Function to take in a tide gauge hourly data set and return a time series of
## peaks-over-threshold local sea levels, together with the corresponding year.
##
## filename_in = string, file name
## data_dir = string, directory path to the file nae, relative to current directory
## threshold_missing_data = filter out days that are missing more data than this (%)
## gpd_threshold = percentile (0-1) to serve as the GPD threshold (only modeling exceedances of this threshold)
## dt_decluster = declustering time-scale (days)
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

process_gpd <- function(filename_in, data_dir, threshold_missing_data=0.9, gpd_threshold=0.99, dt_decluster=3) {

  tbeg <- proc.time()

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

  # where are there gaps longer than one hour? (+10sec for precision)
  igap <- which(time.diff > (one.hours+10/(24*60*60)))

  # Detrend by either subtracting annual means (moving 1-year window)
  print(paste('Detrending by subtracting moving window 1-year average ...', sep=''))

  # get a placeholder
  data.tg$sl.detrended <- data.tg$sl
  time.days.beg <- min(data.tg$time.days)
  time.days.end <- max(data.tg$time.days)

  pb <- txtProgressBar(min=0,max=length(data.tg$time.days),initial=0,style=3)
  for (tt in 1:length(data.tg$time.days)) {
    # if within half a year of either end of the time series, include either the
    # entire first year or entire last year to get a full year's worth of data in
    # the subtracted mean
    if (data.tg$time.days[tt] - time.days.beg < (365.25*0.5)) {
      ind.close <- which(data.tg$time.days - time.days.beg <= 365.25)
    } else if(time.days.end - data.tg$time.days[tt] < (365.25*0.5)) {
      ind.close <- which(time.days.end - data.tg$time.days <= 365.25)
    } else {
      ind.close <- which(abs(data.tg$time.days-data.tg$time.days[tt]) <= (365.25*0.5) )
    }
    data.tg$sl.detrended[tt] <- data.tg$sl[tt] - mean(data.tg$sl[ind.close])
    setTxtProgressBar(pb, tt)
  }
  close(pb)

here

  # daily block maxima; calculate 99% quantile as GPD threshold

  # how many days in each year have at least 90% of their values?
  days.all <- floor(data.tg$time.days)
  days.unique <- unique(days.all)
  ind.days.to.remove <- NULL
  print('... filtering down to do a daily maxima time series of only the days with at least 90% of data ...')
  pb <- txtProgressBar(min=min(days.unique),max=max(days.unique),initial=0,style=3)
  for (day in days.unique) {
    ind.today <- which(floor(data.tg$time.days) == day)
    perc.data.today <- length(ind.today)/24
    if(perc.data.today < 0.9) {ind.days.to.remove <- c(ind.days.to.remove, match(day, days.unique))}
    setTxtProgressBar(pb, day)
  }
  close(pb)
  days.daily.max <- days.unique[-ind.days.to.remove]
  n.days <- length(days.daily.max)

  # calculate the daily maximum sea levels on the days of 'days.daily.max'
  sl.daily.max <- rep(NA, n.days)
  years.daily.max <- rep(NA, n.days)
  print('... calculating time series of daily maxima ...')
  pb <- txtProgressBar(min=0,max=n.days,initial=0,style=3)
  for (day in days.daily.max) {
    cnt <- match(day,days.daily.max)
    ind.today <- which(days.all == day)
    sl.daily.max[cnt] <- max(data.tg$sl.detrended[ind.today])
    years.daily.max[cnt] <- data.tg$year[ind.today][1]
    setTxtProgressBar(pb, cnt)
  }
  close(pb)

  # find all the excesses, "declustering" = if two are within a day of each
  # other, take only the maximum of the two (so make sure you save the times
  # of each excess)

  print('... getting threshold excesses ...')

  # Make a list object for output for this site with everything we need
  # to calibrate the PP-GPD model
  data_out <- vector('list', 6)
  names(data_out) <- c('counts','year','time_length','excesses','threshold','p.threshold','dt.decluster')

  # threshold is 99% quantile of tide gauge's observed values.
  # Buchanan et al (2016) use the 99% quantile of the daily maximum time series.
  gpd.threshold <- as.numeric(quantile(sl.daily.max, gpd_threshold, na.rm=TRUE))
  data_out$dt.decluster <- dt_decluster
  data_out$p.threhsold <- gpd_threshold
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

  # alternatively, could bin em all together. but this won't allow for potential
  # non-stationary behavior in the poisson process
  data_out$excesses_all <- sl.exceed.decl
  data_out$counts_all <- length(sl.exceed.decl)
  data_out$time_length_all <- length(days.daily.max)

  # that takes some time, so save the workspace image after each data set
#  save.image(file=filename.saveprogress)

  tend <- proc.time()
  print(paste('  ... done. Took ', (tend[3]-tbeg[3])/60, ' minutes.',sep=''))

  return(data_out)
}

##==============================================================================
## End
##==============================================================================
