# Output files

Output files from the calibration and projections. The date-stamps on the files given below are from the
originals in Wong et al. (2020; in review), but can give you a sense of what is in each of them.

* `returnlevels_31Jan2022.RData` - return level and return period projections for all of the sites, for all of the covariates and for all of the candidate model structures (8 each of GEV and GPD)
* `comparison_metrics_31Jan2022.RData` - contain the goodness-of-fit metrics for comparing the performance of the alternative models (AIC, BIC, negative log-likelihood and negative log-posterior score)
* `optim_gpd_like_28Jan2022` - differential evolution optimization output for GPD models, to optimize based on the negative log-likelihood
* `optim_gpd_post_28Jan2022` - differential evolution optimization output for GPD models, to optimize based on the negative log-posterior score
* `optim_gev_like_28Jan2022` - differential evolution optimization output for GEV models, to optimize based on the negative log-likelihood
* `optim_gev_post_28Jan2022` - differential evolution optimization output for GEV models, to optimize based on the negative log-posterior score

## Copyright

 This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

## Questions?

Please do not hesitate to contact me if there is anything I can help you with.

Sincerely, Tony Wong (aewsma@rit.edu)
