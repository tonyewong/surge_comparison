## parameter_setup_gpd.R

nmodel <- 8
gpd_models <- vector("list", nmodel)

# set up parameter names for each model, and
# set up parameter bounds for each model (need for optimization for priors)

# all stationary
gpd_models[[1]]$parnames <- c("lambda", "sigma", "xi")
gpd_models[[1]]$bound_lower <- c(0, 0, -3)
gpd_models[[1]]$bound_upper <- c(0.05, 1000, 3)

# lambda nonstationary
gpd_models[[2]]$parnames <- c("lambda0", "lambda1", "sigma", "xi")
gpd_models[[2]]$bound_lower <- c(0, -0.1, 0, -3)
gpd_models[[2]]$bound_upper <- c(0.05, 0.1, 1000, 3)

# sigma nonstationary
gpd_models[[3]]$parnames <- c("lambda", "sigma0", "sigma1", "xi")
gpd_models[[3]]$bound_lower <- c(0, 0, -200, -3)
gpd_models[[3]]$bound_upper <- c(0.05, 1000, 200, 3)

# xi nonstationary
gpd_models[[4]]$parnames <- c("lambda", "sigma", "xi0", "xi1")
gpd_models[[4]]$bound_lower <- c(0, 0, -3, -3)
gpd_models[[4]]$bound_upper <- c(0.05, 1000, 3, 3)

# lambda and sigma nonstationary
gpd_models[[5]]$parnames <- c("lambda0", "lambda1", "sigma0", "sigma1", "xi")
gpd_models[[5]]$bound_lower <- c(0, -0.1, 0, -200, -3)
gpd_models[[5]]$bound_upper <- c(0.05, 0.1, 1000, 200, 3)

# lambda and xi nonstationary
gpd_models[[6]]$parnames <- c("lambda0", "lambda1", "sigma", "xi0", "xi1")
gpd_models[[6]]$bound_lower <- c(0, -0.1, 0, -3, -3)
gpd_models[[6]]$bound_upper <- c(0.05, 0.1, 1000, 3, 3)

# sigma and xi nonstationary
gpd_models[[7]]$parnames <- c("lambda", "sigma0", "sigma1", "xi0", "xi1")
gpd_models[[7]]$bound_lower <- c(0, 0, -200, -3, -3)
gpd_models[[7]]$bound_upper <- c(0.05, 1000, 200, 3, 3)

# all nonstationary
gpd_models[[8]]$parnames <- c("lambda0", "lambda1", "sigma0", "sigma1", "xi0", "xi1")
gpd_models[[8]]$bound_lower <- c(0, -0.1, 0, -200, -3, -3)
gpd_models[[8]]$bound_upper <- c(0.05, 0.1, 1000, 200, 3, 3)
