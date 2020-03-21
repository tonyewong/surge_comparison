# Useful:

# function to trim forcing to fit TG record unique years
# there might be missing years in TG record, so need to match each year and not
# just plop down an evenly spaced sequence
trimmed_forcing <- function(years_tidegauge, years_forcing, forcing) {
  output <- vector('list', 2); names(output) <- c('time','forcing')
  # check the beginning
  if(years_forcing[1] > years_tidegauge[1]) {print('ERROR - tide gauge record starts before forcing; add support for this situation')}
  # check the end
  if(max(years_forcing) < max(years_tidegauge)) {print('ERROR - tide gauge record ends after forcing; add support for this situation')}
  # match the indices of years_tidegauge within years_forcing
  imatch <- match(years_tidegauge, years_forcing)
  output$time <- years_forcing[imatch]
  output$forcing <- forcing[imatch]
  return(output)
}
