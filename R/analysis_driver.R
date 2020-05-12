##==============================================================================
## analysis_driver.R
##
## To be run after calibration_driver.R
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

rm(list=ls())

library(extRemes)
library(Bolstad)

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/surge_comparison/R')
} else if(Sys.info()['user']=='aewsma') {
  machine <- 'office'
  setwd('/Users/aewsma/codes/surge_comparison/R')
} else {
  # assume on another cluster of some kind...
  machine <- 'remote'
  setwd('~/work/codes/surge_comparison/R')
}

## Set up. Make sure consistent with the
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename_covariates <- "../input_data/covariates_22Mar2020.rds"
filename_calibdata_gev  <- "../input_data/processeddata_gev_22Mar2020.rds"
filename_calibdata_gpd  <- "../input_data/processeddata_gpd_23Mar2020.rds"
filename_optim_gev <- "../output/optim_gev-22Mar2020.rds"
filename_optim_gpd <- "../output/optim_gpd-23Mar2020.rds"
names_evm <- c("gev", "gpd")
n_evm <- length(names_evm)

## do some analyzing
source("best_models.R") # yields parameters[[model]][site,parameter]
source("make_projections.R") # yields return level estimates
metrics <- vector("list", 3); names(metrics) <- c("NLL","AIC","BIC")
metrics$NLL <- nll
metrics$AIC <- aic
metrics$BIC <- bic

## number of years of data for each site
num_years <- rep(NA, nsite); names(num_years) <- site_names
for (dd in 1:nsite) {num_years[dd] <- nrow(data_calib$gev[[dd]])}


## get to the good stuff!


##==============================================================================
## Compare model choices based
## on covariate
##============================

##  comparison of model choice based on penalization of over-fitting
best_models_all <- vector("list", 2); names(best_models_all) <- names_evm
for (gg in names_evm) {best_models_all[[gg]] <- vector("list", 3); names(best_models_all[[gg]]) <- names(metrics)}
model_choices <- c(1,2,3,4,5,6,7,8)
covariate_labels <- c("Time","NAO","Temp","Sea level","Stat")
covariate_breaks <- sort(c(seq(from=0.75, to=4.75, by=1), seq(from=1.25, to=5.25, by=1)))

pdf('../figures/covariate_choice.pdf', width=5.5, height=6, pointsize=11, colormodel='cmyk')
par(mfrow=c(3,2), mai=c(.5,.45,.2,.08))
panel_labels <- c('a','b','c','d','e','f'); label_cnt <- 1
for (metric in names(metrics)) {
  for (gg in names_evm) {
    best_models_all[[gg]][[metric]] <- rep(NA, nsite)
    names(best_models_all[[gg]][[metric]]) <- site_names
    for (dd in site_names) {
        mat <- metrics[[metric]][[gg]][dd,,model_choices]
        idx_min <- which.min(mat)
        idx_row <- idx_min%%nrow(mat)
        if(idx_row==0) {idx_row <- nrow(mat)}
        idx_col <- which.min(as.vector(mat[idx_row,]))
        if(idx_col==1) {idx_row <- 5} # stationary model
        best_models_all[[gg]][[metric]][dd] <- idx_row
    }
    hist(best_models_all[[gg]][[metric]], breaks=covariate_breaks, freq=TRUE, xlim=c(0.5,5.5), ylim=c(0,25), main="", xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i')
    grid()
    hist(best_models_all[[gg]][[metric]], breaks=covariate_breaks, freq=TRUE, col="white", xlim=c(0.5,5.5), ylim=c(0,25), main="", xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i', add=TRUE)
    axis(1, at=1:5, labels=covariate_labels)
    axis(2, at=seq(from=0,to=30,by=5), labels=seq(from=0,to=30,by=5))
    mtext(side=1, text="Covariate", line=2.3, cex=0.85)
    mtext(side=2, text="# stations", line=2.3, cex=0.85)
    mtext(side=3, text=metric, cex=0.85)
    mtext(side=3, text=panel_labels[label_cnt], adj=0); label_cnt <- label_cnt + 1
  }
}
dev.off()

##==============================================================================



##==============================================================================
## Compare model choices based
## on covariate and
## amount of available data
##============================

# note: none of the stations (as of this writing) have num_years = median(num_years) (55 years)
idx_short <- which(num_years <= median(num_years))
idx_long <- which(num_years >= median(num_years))

##  comparison of model choice based on penalization of over-fitting
best_models_datlen <- vector("list", 2); names(best_models_datlen) <- c("short", "long")
for (dl in 1:2) {best_models_datlen[[dl]] <- vector("list", 2); names(best_models_datlen[[dl]]) <- names_evm}
for (gg in names_evm) {best_models_datlen[[dl]][[gg]] <- vector("list", 3); names(best_models_datlen[[dl]][[gg]]) <- names(metrics)}
model_choices <- c(1,2,3,4,5,6,7,8)
covariate_breaks_short <- sort(c(seq(from=0.75, to=4.75, by=1), seq(from=1.001, to=5.001, by=1)))
covariate_breaks_long <- sort(c(seq(from=0.999, to=4.999, by=1), seq(from=1.25, to=5.25, by=1)))

pdf('../figures/covariate_choice_datlen.pdf', width=5.5, height=6, pointsize=11, colormodel='cmyk')
par(mfrow=c(3,2), mai=c(.5,.45,.2,.08))
panel_labels <- c('a','b','c','d','e','f'); label_cnt <- 1
for (metric in names(metrics)) {
  for (gg in names_evm) {
    best_models_datlen$short[[gg]][[metric]] <- rep(NA, length(idx_short)); names(best_models_datlen$short[[gg]][[metric]]) <- site_names[idx_short]
    best_models_datlen$long[[gg]][[metric]] <- rep(NA, length(idx_long)); names(best_models_datlen$long[[gg]][[metric]]) <- site_names[idx_long]
    for (dd in site_names) {
      mat <- metrics[[metric]][[gg]][dd,,model_choices]
      idx_min <- which.min(mat)
      idx_row <- idx_min%%nrow(mat)
      if(idx_row==0) {idx_row <- nrow(mat)}
      idx_col <- which.min(as.vector(mat[idx_row,]))
      if(idx_col==1) {idx_row <- 5} # stationary model
      if (dd %in% names(idx_short)) {best_models_datlen$short[[gg]][[metric]][dd] <- idx_row} else {best_models_datlen$long[[gg]][[metric]][dd] <- idx_row}
    }
    hist(best_models_datlen$short[[gg]][[metric]], breaks=covariate_breaks_short, freq=TRUE, xlim=c(0.5,5.5), ylim=c(0,15), main="", xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i')
    grid()
    hist(best_models_datlen$short[[gg]][[metric]], breaks=covariate_breaks_short, freq=TRUE, col="white", xlim=c(0.5,5.5), ylim=c(0,15), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    hist(best_models_datlen$long[[gg]][[metric]], breaks=covariate_breaks_long, freq=TRUE, col="gray50", xlim=c(0.5,5.5), ylim=c(0,15), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    if (label_cnt==1) {legend(0.6,15,c("< 55 years","> 55 years"),pch=c(0,15),col=c("black","gray50"),pt.cex=2,bty='n')}
    axis(1, at=1:5, labels=covariate_labels)
    axis(2, at=seq(from=0,to=30,by=5), labels=seq(from=0,to=30,by=5))
    mtext(side=1, text="Covariate", line=2.3, cex=0.85)
    mtext(side=2, text="# stations", line=2.3, cex=0.85)
    mtext(side=3, text=metric, cex=0.85)
    mtext(side=3, text=panel_labels[label_cnt], adj=0); label_cnt <- label_cnt + 1
  }
}
dev.off()

##==============================================================================



##==============================================================================
## Compare model choices based
## on number of parameters
## NB: this version lumps all of the potential covariates in together
##============================

##  comparison of model choice based on penalization of over-fitting
best_models_all <- vector("list", 2); names(best_models_all) <- names_evm
for (gg in names_evm) {best_models_all[[gg]] <- vector("list", 3); names(best_models_all[[gg]]) <- names(metrics)}
model_choices <- c(1,2,3,4,5,6,7,8)
model_labels <- vector('list', 2); names(model_labels) <- names_evm
model_labels$gev <- c("none",expression(mu),expression(sigma),expression(xi),expression(mu * ', ' * sigma),expression(mu * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * mu * ', ' * sigma * ', ' * xi))
model_labels$gpd <- c("none",expression(lambda),expression(sigma),expression(xi),expression(lambda * ', ' * sigma),expression(lambda * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * lambda * ', ' * sigma * ', ' * xi))
model_breaks <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.4, to=8.4, by=1)))

pdf('../figures/model_choice.pdf', width=6, height=6, pointsize=11, colormodel='cmyk')
par(mfrow=c(3,2), mai=c(.5,.45,.2,.08))
panel_labels <- c('a','b','c','d','e','f'); label_cnt <- 1
for (metric in names(metrics)) {
  for (gg in names_evm) {
    best_models_all[[gg]][[metric]] <- rep(NA, nsite)
    names(best_models_all[[gg]][[metric]]) <- site_names
    for (dd in site_names) {
        mat <- metrics[[metric]][[gg]][dd,,model_choices]
        idx_min <- which.min(mat)
        idx_row <- idx_min%%nrow(mat)
        if(idx_row==0) {idx_row <- nrow(mat)}
        idx_col <- which.min(as.vector(mat[idx_row,]))
        best_models_all[[gg]][[metric]][dd] <- idx_col
    }
    hist(best_models_all[[gg]][[metric]], breaks=model_breaks, freq=TRUE, xlim=c(0.5,8.5), ylim=c(0,26), main="", xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i')
    grid()
    hist(best_models_all[[gg]][[metric]], breaks=model_breaks, freq=TRUE, col="white", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    axis(1, at=1:8, labels=model_labels[[gg]])
    axis(2, at=seq(from=0,to=30,by=5), labels=seq(from=0,to=30,by=5))
    mtext(side=1, text="Nonstationary parameters", line=2.3, cex=0.85)
    mtext(side=2, text="# stations", line=2.3, cex=0.85)
    mtext(side=3, text=metric, cex=0.85)
    mtext(side=3, text=panel_labels[label_cnt], adj=0); label_cnt <- label_cnt + 1
  }
}
dev.off()

##==============================================================================



##==============================================================================
## Compare model choices based
## on number of parameters and
## amount of available data
## NB: this version lumps all of the potential covariates in together
##============================

##  comparison of model choice based on penalization of over-fitting
best_models_datlen <- vector("list", 2); names(best_models_datlen) <- c("short", "long")
for (dl in 1:2) {best_models_datlen[[dl]] <- vector("list", 2); names(best_models_datlen[[dl]]) <- names_evm}
for (gg in names_evm) {best_models_datlen[[dl]][[gg]] <- vector("list", 3); names(best_models_datlen[[dl]][[gg]]) <- names(metrics)}
model_choices <- c(1,2,3,4,5,6,7,8)
model_breaks_short <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.001, to=8.001, by=1)))
model_breaks_long <- sort(c(seq(from=0.999, to=7.999, by=1), seq(from=1.4, to=8.4, by=1)))

pdf('../figures/model_choice_datlen.pdf', width=6, height=6, pointsize=11, colormodel='cmyk')
par(mfrow=c(3,2), mai=c(.5,.45,.2,.08))
panel_labels <- c('a','b','c','d','e','f'); label_cnt <- 1
for (metric in names(metrics)) {
  for (gg in names_evm) {
    best_models_datlen$short[[gg]][[metric]] <- rep(NA, length(idx_short)); names(best_models_datlen$short[[gg]][[metric]]) <- site_names[idx_short]
    best_models_datlen$long[[gg]][[metric]] <- rep(NA, length(idx_long)); names(best_models_datlen$long[[gg]][[metric]]) <- site_names[idx_long]
    for (dd in site_names) { # yeah, redoing this...
        mat <- metrics[[metric]][[gg]][dd,,model_choices]
        idx_min <- which.min(mat)
        idx_row <- idx_min%%nrow(mat)
        if(idx_row==0) {idx_row <- nrow(mat)}
        idx_col <- which.min(as.vector(mat[idx_row,]))
        if (dd %in% names(idx_short)) {best_models_datlen$short[[gg]][[metric]][dd] <- idx_col} else {best_models_datlen$long[[gg]][[metric]][dd] <- idx_col}
    }
    hist(best_models_datlen$short[[gg]][[metric]], breaks=model_breaks_short, freq=TRUE, xlim=c(0.5,8.5), ylim=c(0,15), main="", xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i')
    grid()
    hist(best_models_datlen$short[[gg]][[metric]], breaks=model_breaks_short, freq=TRUE, col="white", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    hist(best_models_datlen$long[[gg]][[metric]], breaks=model_breaks_long, freq=TRUE, col="gray50", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    if (label_cnt==1) {legend(0.5,15,c("< 55 years","> 55 years"),pch=c(0,15),col=c("black","gray50"),pt.cex=2,bty='n')}
    axis(1, at=1:8, labels=model_labels[[gg]])
    axis(2, at=seq(from=0,to=30,by=5), labels=seq(from=0,to=30,by=5))
    mtext(side=1, text="Nonstationary parameters", line=2.3, cex=0.85)
    mtext(side=2, text="# stations", line=2.3, cex=0.85)
    mtext(side=3, text=metric, cex=0.85)
    mtext(side=3, text=panel_labels[label_cnt], adj=0); label_cnt <- label_cnt + 1
  }
}
dev.off()

##==============================================================================




##==============================================================================
## Compare model choices based
## on number of parameters
## NB: this version uses only one of the covariates. Cycles through them.
##============================

for (cc in names_covariates) {
  ##  comparison of model choice based on penalization of over-fitting
  best_models_all <- vector("list", 2); names(best_models_all) <- names_evm
  for (gg in names_evm) {best_models_all[[gg]] <- vector("list", 3); names(best_models_all[[gg]]) <- names(metrics)}
  model_choices <- c(1,2,3,4,5,6,7,8)
  model_labels <- vector('list', 2); names(model_labels) <- names_evm
  model_labels$gev <- c("none",expression(mu),expression(sigma),expression(xi),expression(mu * ', ' * sigma),expression(mu * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * mu * ', ' * sigma * ', ' * xi))
  model_labels$gpd <- c("none",expression(lambda),expression(sigma),expression(xi),expression(lambda * ', ' * sigma),expression(lambda * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * lambda * ', ' * sigma * ', ' * xi))
  model_breaks <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.4, to=8.4, by=1)))

  pdf(paste('../figures/model_choice_',cc,'.pdf',sep=''), width=6, height=6, pointsize=11, colormodel='cmyk')
  par(mfrow=c(3,2), mai=c(.5,.45,.2,.08))
  panel_labels <- c('a','b','c','d','e','f'); label_cnt <- 1
  for (metric in names(metrics)) {
    for (gg in names_evm) {
      best_models_all[[gg]][[metric]] <- rep(NA, nsite)
      names(best_models_all[[gg]][[metric]]) <- site_names
      for (dd in site_names) {
          mat <- metrics[[metric]][[gg]][dd,cc,model_choices]
          idx_col <- which.min(mat)
          best_models_all[[gg]][[metric]][dd] <- idx_col
      }
      hist(best_models_all[[gg]][[metric]], breaks=model_breaks, freq=TRUE, xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i')
      grid()
      hist(best_models_all[[gg]][[metric]], breaks=model_breaks, freq=TRUE, col="white", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
      axis(1, at=1:8, labels=model_labels[[gg]])
      axis(2, at=seq(from=0,to=30,by=5), labels=seq(from=0,to=30,by=5))
      mtext(side=1, text="Nonstationary parameters", line=2.3, cex=0.85)
      mtext(side=2, text="# stations", line=2.3, cex=0.85)
      mtext(side=3, text=metric, cex=0.85)
      mtext(side=3, text=panel_labels[label_cnt], adj=0); label_cnt <- label_cnt + 1
    }
  }
  dev.off()
} # end covariates loop

##==============================================================================



##==============================================================================
## Compare model choices based
## on number of parameters and
## amount of available data
## NB: this version uses only NAO covariate
##============================

for (cc in names_covariates) {
  ##  comparison of model choice based on penalization of over-fitting
  best_models_datlen <- vector("list", 2); names(best_models_datlen) <- c("short", "long")
  for (dl in 1:2) {best_models_datlen[[dl]] <- vector("list", 2); names(best_models_datlen[[dl]]) <- names_evm}
  for (gg in names_evm) {best_models_datlen[[dl]][[gg]] <- vector("list", 3); names(best_models_datlen[[dl]][[gg]]) <- names(metrics)}
  model_choices <- c(1,2,3,4,5,6,7,8)
  model_breaks_short <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.001, to=8.001, by=1)))
  model_breaks_long <- sort(c(seq(from=0.999, to=7.999, by=1), seq(from=1.4, to=8.4, by=1)))

  pdf(paste('../figures/model_choice_datlen_',cc,'.pdf',sep=''), width=6, height=6, pointsize=11, colormodel='cmyk')
  par(mfrow=c(3,2), mai=c(.5,.45,.2,.08))
  panel_labels <- c('a','b','c','d','e','f'); label_cnt <- 1
  for (metric in names(metrics)) {
    for (gg in names_evm) {
      best_models_datlen$short[[gg]][[metric]] <- rep(NA, length(idx_short)); names(best_models_datlen$short[[gg]][[metric]]) <- site_names[idx_short]
      best_models_datlen$long[[gg]][[metric]] <- rep(NA, length(idx_long)); names(best_models_datlen$long[[gg]][[metric]]) <- site_names[idx_long]
      for (dd in site_names) { # yeah, redoing this...
        mat <- metrics[[metric]][[gg]][dd,cc,model_choices]
        idx_col <- which.min(mat)
        if (dd %in% names(idx_short)) {best_models_datlen$short[[gg]][[metric]][dd] <- idx_col} else {best_models_datlen$long[[gg]][[metric]][dd] <- idx_col}
      }
      hist(best_models_datlen$short[[gg]][[metric]], breaks=model_breaks_short, freq=TRUE, xlim=c(0.5,8.5), ylim=c(0,17), main="", xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i')
      grid()
      hist(best_models_datlen$short[[gg]][[metric]], breaks=model_breaks_short, freq=TRUE, col="white", xlim=c(0.5,8.5), ylim=c(0,17), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
      hist(best_models_datlen$long[[gg]][[metric]], breaks=model_breaks_long, freq=TRUE, col="gray50", xlim=c(0.5,8.5), ylim=c(0,17), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
      if (label_cnt==1) {legend(0.5,15,c("< 55 years","> 55 years"),pch=c(0,15),col=c("black","gray50"),pt.cex=2,bty='n')}
      axis(1, at=1:8, labels=model_labels[[gg]])
      axis(2, at=seq(from=0,to=30,by=5), labels=seq(from=0,to=30,by=5))
      mtext(side=1, text="Nonstationary parameters", line=2.3, cex=0.85)
      mtext(side=2, text="# stations", line=2.3, cex=0.85)
      mtext(side=3, text=metric, cex=0.85)
      mtext(side=3, text=panel_labels[label_cnt], adj=0); label_cnt <- label_cnt + 1
    }
  }
  dev.off()
}

##==============================================================================



##==============================================================================
## Compare return levels based
## on model choice
##============================

todo

##==============================================================================



##==============================================================================
## Map of covariate choice based on location
##==========================================

todo

site_dat <- read.csv("../input_data/SiteInformation.csv", header=TRUE)
idx_map <- match(site_names, site_dat[,"site_name"])
site_dat <- site_dat[idx_map,]
lats <- site_dat[,"lat"]
lons <- site_dat[,"lon"]

best_models_map <- vector("list", 2); names(best_models_all) <- names_evm
for (gg in names_evm) {best_models_all[[gg]] <- vector("list", length(metrics)); names(best_models_all[[gg]]) <- names(metrics)}
model_choices <- c(1,2,3,4,5,6,7,8)
for (metric in names(metrics)) {
  for (gg in names_evm) {
    best_models_map[[gg]][[metric]] <- rep(NA, nsite)
    names(best_models_map[[gg]][[metric]]) <- site_names
    for (dd in site_names) {
        mat <- metrics[[metric]][[gg]][dd,,model_choices]
        idx_min <- which.min(mat)
        idx_row <- idx_min%%nrow(mat)
        if(idx_row==0) {idx_row <- nrow(mat)}
        idx_col <- which.min(as.vector(mat[idx_row,]))
        if(idx_col==1) {idx_row <- 5} # stationary model
        best_models_map[[gg]][[metric]][dd] <- idx_row
    }
  }
}

## site choices  #covariate_labels
metric <- "AIC" # which metric to use for the plot?
covariate_choices_site <- vector("list", 5)
for (cc in 1:5) {covariate_choices_site[[cc]] <- which(best_models_map[[gg]][[metric]]==cc)}
covariate_colors <- vector("list", 5); names(covariate_colors) <- covariate_labels
covariate_colors$Time <- "darkorange3"
covariate_colors$NAO <- "seagreen"
covariate_colors$Temp <- "mediumslateblue"
covariate_colors$`Sea level` <- "mediumvioletred"
covariate_colors$Stat <- "black"

## Make a plot

library(maps)

png("../covariate_choice_map.png", width=600, height=600)
par(mfrow=c(1,1))
map("state", fill=TRUE, col="gray85", bg="white", xlim=c(-105, -64), ylim=c(22, 48), mar=c(8,12,0,0))
for (cc in 1:5) {points(lons[covariate_choices_site[[cc]]], lats[covariate_choices_site[[cc]]], col=covariate_colors[[cc]], cex=2, pch=16)}
mtext(side=1, text='Longitude', line=3.5, cex=1.5)
mtext(side=2, text='Latitude', line=4.5, cex=1.5)
axis(side=1, at=seq(from=-100, to=-64, by=5), labels=c("100 W", "", "90 W", "", "80 W", "", "70 W", ""), cex.axis=1.5)
axis(side=2, at=seq(from=20, to=50, by=5), labels=c("20 N", "25 N", "30 N", "35 N", "40 N", "45 N", "50 N"), las=1, cex.axis=1.5)
dev.off()

##==============================================================================



##==============================================================================
## Map of model choice based on location
## (all covariates together? or just one?)
##======================================

todo

##==============================================================================



##==============================================================================
## End
##==============================================================================
