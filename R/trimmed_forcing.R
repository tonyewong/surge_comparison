##==============================================================================
## trimmed_forcing.R
##
## This is a helper function to trim the auxiliary covariate forcing to fit the
## tide gauge record's unique years. There might be missing years in the tide
## gauge data set, so we need to match each year and not just plot down an
## evenly spaced sequence. Does it sound like I'm speaking from experience?
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

##==============================================================================
## End
##==============================================================================
