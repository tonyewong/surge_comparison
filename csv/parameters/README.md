# Optimized parameter sets

## Directory structure

The first layer in the `parameters` directory here is the choice of objective function to maximize. There are two choices for this:
* `nll`, the negative log-likelihood (thus, you are looking at the maximum likelihood parameter estimates), or
* `nps`, the negative posterior score (thus, you are looking at the maximum a posteriori parameter estimates).

Once you decide which objective function to use, go into that directory. The second layer is which extreme value statistical model form to use:
* `gev` for generalized extreme value distributions and annual block maxima, or
* `gpd` for a Poisson process/generalized Pareto distribution model, for a time series of peaks-over-thresholds.

Decide which extreme value model to use and enter that directory. The third layer is which covariate time series to use to modulate the nonstationary parameters (if any) in the extreme value model. Within each of these choices, one of the models corresponds to a stationary model. The choices for covariate time series are:
* `time`, a linear increase in time,
* `temp`, global mean surface temperature,
* `sealevel`, global mean sea level, or
* `nao`, winter mean North Atlantic Oscillation index.

Within each of those covariate directories will be 8 CSV files. Each corresponds to one of the different model structures, shown in Table 2 of the Wong et al. (under review) manuscript that this repository accompanies. These models are:
* `model1`: stationary
* `model2`: $\mu$ (GEV) or $\lambda$ (GPD) nonstationary, but $\sigma$ and $\xi$ stationary
* `model3`: $\sigma$ nonstationary, but $\mu$ ($\lambda$) and $\xi$ stationary
* `model4`: $\xi$ nonstationary, but $\mu$ ($\lambda$) and $\sigma$ stationary
* `model5`: $\mu$ ($\lambda$) and $\sigma$ nonstationary, but $\xi$ stationary
* `model6`: $\mu$ ($\lambda$) and $\xi$ nonstationary, but $\sigma$ stationary
* `model7`: $\sigma$ and $\xi$ nonstationary, but $\mu$ ($\lambda$) stationary
* `model8`: all three parameters nonstationary

---

Questions? Tony Wong (aewsma@rit.edu)
