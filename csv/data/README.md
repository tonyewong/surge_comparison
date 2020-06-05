# Processed tide gauge data sets

## Directory structure

The first layer in the `data` directory here is the choice of which extreme value statistical model form to use:
* `gev` for generalized extreme value distributions and annual block maxima, or
* `gpd` for a Poisson process/generalized Pareto distribution model, for a time series of peaks-over-thresholds.

Decide which extreme value model to use and enter that directory. There will be a separate file for each of the tide gauge stations.

## GEV

For the GEV sets of data:
* first column corresponds to the year
* second column corresponds to the detrended annual maximum local sea level
* do not assume that every year will be included in the data sets! If a year is missing more than 90% of its hourly data, then that year is discarded from the analysis entirely.

## GPD

For the Poisson process/GPD sets of data, each data set is first detrended by subtracting off a moving window 1-year average. Where there are gaps, they are filled by fitting a mean annual cycle with a linear trend (described in the accompanying manuscript). But, the gap-filling is only used for the detrending; no portions of missing data are used for analysis.

We then compute the time series of daily maxima for each station. We discard any days that are missing 4 hours or more (90%) of the data for that day. We use the 99th percentile as the threshold to designate an event as "extreme" and count as an exceedance for the peaks-over-thresholds GPD model. Finally, we use a 3-day declustering time-scale to remove events that are likely part of the same storm. If two events are within 3 days of one another, we remove the event with the lower sea level signal.

The file `thresholds_declustering.csv` gives for each site:
* column 1: the percentile used for the peaks-over-thresholds threshold
* column 2: the sea level height corresponding to the extreme sea level threshold
* column 3: the declustering time-scale (days)

In the files for each site:
* column 1: year
* column 2: number of threshold exceedances in that year
* column 3: amount of time (days) of data present for that year
* column 4 and higher: sea level heights of each of the exceedances for that year.

---

Questions? Tony Wong (aewsma@rit.edu)
