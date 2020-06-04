# surge_comparison

## Directory structure

todo

## Workflow

todo

1. process_data.R
  1. process_gev.R - yields `processeddata_gev_[date].rds`
  1. get_timeseries_covariates.R - yields `covariates_[date].rds`
1. calibration_driver.R
  1. trimmed_forcing.R
  1. likelihood_gev.R
  1. parameter_setup_gev.R
  1. yields `optim_[covariate name]-gev-[date].rds`, output from maximum likelihood optimization for GEV parameters in each of the 8 potentially nonstationary models.


## Input data

### Tide gauge stations

* Only take stations with at least 15 years of available data

### Covariate time series forcing

* time is simply the year
* global mean sea level is taken from [Wong and Keller (2017; doi: 10.1002/2017EF000607)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1002/2017EF000607)) for projections and from [Church and White (2013; doi: 10.1007/s10712-011-9119-1)](https://link.springer.com/article/10.1007/s10712-011-9119-1) for hindcast
* global mean surface temperature is taken from [todo](todo) for projections and from [todo](todo) for hindcast
* winter mean NAO index is taken from [todo](todo) for hindcast and computed from [todo](todo) for projections

If you want to use your own time series forcing,
* get it onto an annual mean time scale,
* normalize it to have a range of 0-1 for the hindcast period,
* and pop this into the `get_timeseries_covariates.R` script.
* You will probably also have to adjust some of the names and plot/axis legends in the plots in `analysis_driver.R`.

## Copyright

 This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

## Questions?

Please do not hesitate to contact me if there is anything I can help you with.

Sincerely, Tony Wong (aewsma@rit.edu)
