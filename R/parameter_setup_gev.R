##==============================================================================
## parameter_setup_gev.R
##
## Sets up a list `gev_models` that has the parameters names and bounds for each
## of the 8 types of GEV distribution, varying from fully stationary (first list
## element) to fully nonstationary (last list element).
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

nmodel <- 8
gev_models <- vector("list", nmodel)

# set up parameter names for each model, and
# set up parameter bounds for each model (need for optimization for priors)

# all stationary
gev_models[[1]]$parnames <- c("mu", "sigma", "xi")
gev_models[[1]]$bound_lower <- c(0, 0, -3)
gev_models[[1]]$bound_upper <- c(5000, 1000, 3)

# mu nonstationary
gev_models[[2]]$parnames <- c("mu0", "mu1", "sigma", "xi")
gev_models[[2]]$bound_lower <- c(0, -500, 0, -3)
gev_models[[2]]$bound_upper <- c(5000, 500, 1000, 3)

# sigma nonstationary
gev_models[[3]]$parnames <- c("mu", "sigma0", "sigma1", "xi")
gev_models[[3]]$bound_lower <- c(0, 0, -200, -3)
gev_models[[3]]$bound_upper <- c(5000, 1000, 200, 3)

# xi nonstationary
gev_models[[4]]$parnames <- c("mu", "sigma", "xi0", "xi1")
gev_models[[4]]$bound_lower <- c(0, 0, -3, -3)
gev_models[[4]]$bound_upper <- c(5000, 1000, 3, 3)

# mu and sigma nonstationary
gev_models[[5]]$parnames <- c("mu0", "mu1", "sigma0", "sigma1", "xi")
gev_models[[5]]$bound_lower <- c(0, -500, 0, -200, -3)
gev_models[[5]]$bound_upper <- c(5000, 500, 1000, 200, 3)

# mu and xi nonstationary
gev_models[[6]]$parnames <- c("mu0", "mu1", "sigma", "xi0", "xi1")
gev_models[[6]]$bound_lower <- c(0, -500, 0, -3, -3)
gev_models[[6]]$bound_upper <- c(5000, 500, 1000, 3, 3)

# sigma and xi nonstationary
gev_models[[7]]$parnames <- c("mu", "sigma0", "sigma1", "xi0", "xi1")
gev_models[[7]]$bound_lower <- c(0, 0, -200, -3, -3)
gev_models[[7]]$bound_upper <- c(5000, 1000, 200, 3, 3)

# all nonstationary
gev_models[[8]]$parnames <- c("mu0", "mu1", "sigma0", "sigma1", "xi0", "xi1")
gev_models[[8]]$bound_lower <- c(0, -500, 0, -200, -3, -3)
gev_models[[8]]$bound_upper <- c(5000, 500, 1000, 200, 3, 3)

##==============================================================================
## End
##==============================================================================
