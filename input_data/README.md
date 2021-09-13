# Input data files

Input files for processing the data and performing the calibration. The date-stamps on the files given below are from the
originals in Wong et al. (2020; in review), but can give you a sense of what is in each of them.

## Related to the tide gauge sites

* `tidegaugedata` - (directory) in here you will find 40 CSV data sets, one for each of the tide gauge stations considered in this work. Four of them are removed from the analysis because they do not have sufficient data after processing (removing years with missing data). If you want to add new sites to the analysis, just stick the CSV files in here, make sure the columns match those in the existing data sets (and units match!), and use the same file naming convention. You will want to add the site information to the Site Information spreadsheets, described next
* `Sites Information Spreadsheet.xlsx` - tide gauge station information for all 40 of the original stations considered
* `SiteInformation.csv` - tide gauge station information for the 36 sites (out of the original 40) used in the analysis. This is the version of the file that is actually used by the R code

## Related to the covariate forcing time series

* `covariates_22Mar2020.rds` - historical and future projections covariates for forcing the nonstationary models
* `GMSL_ChurchWhite2011_yr_2015.txt` - global mean sea level historical data from [Church and White (2013; doi: 10.1007/s10712-011-9119-1)](https://link.springer.com/article/10.1007/s10712-011-9119-1)
* `nao_3dp.dat` - winter mean NAO index historical data from [Jones et al. (1997)](https://doi.org/10.1002/(SICI)1097-0088(19971115)17:13%3C1433::AID-JOC203%3E3.0.CO;2-P)
* `noaa_temperature_1880-2017.csv` - global mean surface temperature historical data from [National Centers for Environmental Information data portal](http://www.ncdc.noaa.gov/cag/)

The following three files must be downloaded using the commands in `R/install_packages.R` (if you do not have these files already). All of the other files should be in the Github repository already. If you downloaded the zipped folder from Zenodo, they will be in there.
* `BRICK_GMSL_WongKeller2017.nc` - global mean sea level projections from [Wong and Keller (2017; doi: 10.1002/2017EF000607)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1002/2017EF000607)
* `DMIEH5_SRA1B_4_MM_psl.1-1200.nc` - sea-level pressure projection from [MPI-ECHAM5](http://www.mpimet.mpg.de/fileadmin/models/echam/mpi_report_349.pdf)
* `global.tas.aann.CNRM-CM5.historical+rcp85.r1i1p1.18500101-21001231.nc` - global mean surface temperature projection from the [CNRM-CM5 simulation (member 1) under Representative Concentration Pathway 8.5 (RCP8.5)](http://cmip-pcmdi.llnl.gov/cmip5/)

## Related to fitting prior distributions for the model parameters

* `tidegauge_processed_manystations_decl3-pot99-annual_10Dec2017.rds` - this one is a head start on fitting prior distributions, because previous work ([Wong et al., 2018](https://doi.org/10.1088/1748-9326/aacb3d)) used this set of 28 tide gauge stations with long data records to fit priors. This RDS file contains the annual block maxima (`gev_year` element on the lists) and peaks-over-thresholds data (`gpd` list elements) for each site.
* `surge_MLEs_gev_15May2020.rds` - maximum likelihood parameter estimates for all of the GEV model structures for the set of long-record tide gauge stations around the world
* `surge_MLEs_gpd_15May2020.rds` - maximum likelihood parameter estimates for all of the GPD model structures for the set of long-record tide gauge stations around the world
* `surge_priors_uniform_gev_15May2020.rds` - uniform prior distributions fit for the GEV model parameters. Different priors are fit for each of the different types of model (so, for example, just because $\mu_0$ has prior $\pi(\mu_0)$ for one model structure, and that parameter is in another model structure, the prior distribution is not necessarily also $\pi$ in the second model). The same is true for the other prior distribution files too.
* `surge_priors_uniform_gpd_15May2020.rds` - uniform prior distributions fit for the GPD model parameters
* `surge_priors_normalgamma_gev_15May2020.rds` - normal and gamma prior distributions for the GEV model parameters
* `surge_priors_normalgamma_gpd_15May2020.rds` - normal and gamma prior distributions for the GPD model parameters

## Related to the final processed data sets, which serve as input to the calibration

* `processeddata_gev_22Mar2020.rds` - processed annual block maxima data for fitting GEV distributions (with 8 model structural forms and 4 covariates considered)
* `processeddata_gpd_23Mar2020.rds` - processed peaks-over-thresholds data for fitting GPD distributions (with 8 model structures and 4 covariates), using in the default results presented the 99th percentile of daily maxima as the GPD threshold and a 3-day declustering time-scale to separate exceedances from the same event

## Copyright

 This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

## Questions?

Please do not hesitate to contact me if there is anything I can help you with.

Sincerely, Tony Wong (aewsma@rit.edu)
