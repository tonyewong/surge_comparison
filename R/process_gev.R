##==============================================================================
## process_gev.R
##
## Function to take in a tide gauge hourly data set and return a time series of
## annual maximum local sea levels, together with the corresponding year.
##
## filename_in = string, file name
## data_dir = string, directory path to the file nae, relative to current directory
## threshold_missing_data = filter out years that are missing more data than this (%)
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

process_gev <- function(filename_in, data_dir, threshold_missing_data=0.9) {

  fillvalue <- -32767 # fill-value
  data.tg <- read.csv(paste(data_dir, filename_in, sep=''))
  names(data.tg) <- c('year','month','day','hour','sl')

  # any years with less than 90% of the data?
  all_years <- sort(unique(data.tg$year))
  n_data_per_year <- 24*365 # 24 hours/day * 365 days/year (minimum)
  missing_data <- NULL
  for (year in all_years) {
    idx_this_year <- which(data.tg$year==year)
    n_data_this_year <- length(which(data.tg$sl[idx_this_year] > fillvalue))
    if (n_data_this_year/n_data_per_year < threshold_missing_data) {
      missing_data <- c(missing_data, year)
    }
  }
  # for NOLA, missing_data should be NULL
  if (!is.null(missing_data)) {print("WARNING: dropping some years for missing data")}

  # which years are good?
  good_years <- setdiff(all_years, missing_data)
  nyear <- length(good_years)

  # detrend by subtracting off annual means and get block maxima
  data.tg$sl.detrended <- data.tg$sl
  lsl_max <- rep(NA,nyear)  # annual block maxima

  for (year in good_years) {
    idx_this_year <- which(data.tg$year==year & data.tg$sl > fillvalue)
    mean_sl_this_year <- mean(data.tg$sl[idx_this_year])
    data.tg$sl.detrended[idx_this_year] <- data.tg$sl[idx_this_year] - mean_sl_this_year
    lsl_max[which(good_years==year)] <- max(data.tg$sl.detrended[idx_this_year])
  }

  output_matrix <- matrix(NA, nrow=nyear, ncol=2)
  colnames(output_matrix) <- c("year", "lsl_max")
  output_matrix[,"year"] <- good_years
  output_matrix[,"lsl_max"] <- lsl_max

  return(output_matrix)
}

##==============================================================================
## End
##==============================================================================
