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
# calibration using likelihood:
filename_maxlike_gev <- "../output/optim_gev_like_22Mar2020.rds"
filename_maxlike_gpd <- "../output/optim_gpd_like_23Mar2020.rds"
# calibration using posterior (score):
filename_maxpost_gev <- "../output/optim_gev_post_15May2020.rds"
filename_maxpost_gpd <- "../output/optim_gpd_post_15May2020.rds"
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
idx_short <- which(num_years <= median(num_years)) # note: none of the stations (as of this writing) have num_years = median(num_years) (55 years)
idx_long <- which(num_years >= median(num_years))


## get to the good stuff!


##==============================================================================
## Compare model choices based on
## covariate and amount of available data
##============================


##  comparison of model choice based on penalization of over-fitting
best_models_datlen <- vector("list", 2); names(best_models_datlen) <- c("short", "long")
best_models_all <- vector("list", 2); names(best_models_all) <- names_evm
for (dl in 1:2) {
  best_models_datlen[[dl]] <- vector("list", 2); names(best_models_datlen[[dl]]) <- names_evm
  for (gg in names_evm) {
    best_models_datlen[[dl]][[gg]] <- vector("list", 3); names(best_models_datlen[[dl]][[gg]]) <- names(metrics)
    best_models_all[[gg]] <- vector("list", 3); names(best_models_all[[gg]]) <- names(metrics) # yes, doing this multiple times
  }
}
model_choices <- c(1,2,3,4,5,6,7,8)
covariate_labels <- c("Time","Temp","Sea level","NAO","Stat") # needs to be in the same order as the second dimension of the `metrics` array: [1] "time" "temp" "sealevel" "nao"
covariate_breaks <- sort(c(seq(from=0.75, to=4.75, by=1), seq(from=1.25, to=5.25, by=1)))
covariate_breaks_short <- sort(c(seq(from=0.75, to=4.75, by=1), seq(from=1.001, to=5.001, by=1)))
covariate_breaks_long <- sort(c(seq(from=0.999, to=4.999, by=1), seq(from=1.25, to=5.25, by=1)))

pdf('../figures/covariate_choice.pdf', width=5.5, height=6, pointsize=11, colormodel='cmyk')
par(mfrow=c(3,2), mai=c(.5,.45,.2,.08))
panel_labels <- c(expression('(a)'),expression('(b)'),expression('(c)'),expression('(d)'),expression('(e)'),expression('(f)')); label_cnt <- 1
for (metric in names(metrics)) {
  for (gg in names_evm) {
    best_models_datlen$short[[gg]][[metric]] <- rep(NA, length(idx_short)); names(best_models_datlen$short[[gg]][[metric]]) <- site_names[idx_short]
    best_models_datlen$long[[gg]][[metric]] <- rep(NA, length(idx_long)); names(best_models_datlen$long[[gg]][[metric]]) <- site_names[idx_long]
    best_models_all[[gg]][[metric]] <- rep(NA, nsite); names(best_models_all[[gg]][[metric]]) <- site_names
    for (dd in site_names) {
      mat <- metrics[[metric]][[gg]][dd,,model_choices]
      idx_min <- which.min(mat)
      idx_row <- idx_min%%nrow(mat)
      if(idx_row==0) {idx_row <- nrow(mat)}
      idx_col <- which.min(as.vector(mat[idx_row,]))
      if(idx_col==1) {idx_row <- 5} # stationary model
      best_models_all[[gg]][[metric]][dd] <- idx_row
      if (dd %in% names(idx_short)) {best_models_datlen$short[[gg]][[metric]][dd] <- idx_row} else {best_models_datlen$long[[gg]][[metric]][dd] <- idx_row}
    }
    hist(best_models_all[[gg]][[metric]], breaks=covariate_breaks, freq=TRUE, xlim=c(0.5,5.5), ylim=c(0,25), main="", lty=5, xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i')
    grid()
    hist(best_models_datlen$short[[gg]][[metric]], breaks=covariate_breaks_short, freq=TRUE, col="gray80", xlim=c(0.5,5.5), ylim=c(0,15), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    hist(best_models_datlen$long[[gg]][[metric]], breaks=covariate_breaks_long, freq=TRUE, col="gray30", xlim=c(0.5,5.5), ylim=c(0,15), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    if (label_cnt==1) {legend(0.6,25,c("< 55 years","> 55 years","combined"),pch=c(15,15,NA),lty=c(NA,NA,5),col=c("gray80","gray30","black"),pt.cex=2,bty='n')}
    axis(1, at=1:5, labels=covariate_labels)
    axis(2, at=seq(from=0,to=30,by=5), labels=seq(from=0,to=30,by=5), las=1)
    mtext(side=1, text="Covariate", line=2.3, cex=0.85)
    mtext(side=2, text="# stations", line=2.3, cex=0.85)
    mtext(side=3, text=metric, cex=0.85)
    mtext(side=3, text=panel_labels[label_cnt], cex=0.85, adj=0); label_cnt <- label_cnt + 1
  }
}
dev.off()

##==============================================================================



##==============================================================================
## Compare model choices based
## on number of parameters and amount of available data
## NB: this version lumps all of the potential covariates in together
##============================

##  comparison of model choice based on penalization of over-fitting
best_models_datlen <- vector("list", 2); names(best_models_datlen) <- c("short", "long")
best_models_all <- vector("list", 2); names(best_models_all) <- names_evm
for (dl in 1:2) {
  best_models_datlen[[dl]] <- vector("list", 2); names(best_models_datlen[[dl]]) <- names_evm
  for (gg in names_evm) {
    best_models_datlen[[dl]][[gg]] <- vector("list", 3); names(best_models_datlen[[dl]][[gg]]) <- names(metrics)
    best_models_all[[gg]] <- vector("list", 3); names(best_models_all[[gg]]) <- names(metrics)
  }
}
model_choices <- c(1,2,3,4,5,6,7,8)
model_labels <- vector('list', 2); names(model_labels) <- names_evm
model_labels$gev <- c("none",expression(mu),expression(sigma),expression(xi),expression(mu * ', ' * sigma),expression(mu * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * mu * ', ' * sigma * ', ' * xi))
model_labels$gpd <- c("none",expression(lambda),expression(sigma),expression(xi),expression(lambda * ', ' * sigma),expression(lambda * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * lambda * ', ' * sigma * ', ' * xi))
model_breaks_short <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.001, to=8.001, by=1)))
model_breaks_long <- sort(c(seq(from=0.999, to=7.999, by=1), seq(from=1.4, to=8.4, by=1)))
model_breaks <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.4, to=8.4, by=1)))

pdf('../figures/model_choice.pdf', width=6, height=6, pointsize=11, colormodel='cmyk')
par(mfrow=c(3,2), mai=c(.5,.45,.2,.08))
panel_labels <- c(expression('(a)'),expression('(b)'),expression('(c)'),expression('(d)'),expression('(e)'),expression('(f)')); label_cnt <- 1
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
        best_models_all[[gg]][[metric]][dd] <- idx_col
        if (dd %in% names(idx_short)) {best_models_datlen$short[[gg]][[metric]][dd] <- idx_col} else {best_models_datlen$long[[gg]][[metric]][dd] <- idx_col}
    }
    hist(best_models_all[[gg]][[metric]], breaks=model_breaks, freq=TRUE, lty=5, xlim=c(0.5,8.5), ylim=c(0,26), main="", xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i')
    grid()
    hist(best_models_datlen$short[[gg]][[metric]], breaks=model_breaks_short, freq=TRUE, col="white", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    hist(best_models_datlen$long[[gg]][[metric]], breaks=model_breaks_long, freq=TRUE, col="gray50", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    if (label_cnt==1) {legend(0.6,25,c("< 55 years","> 55 years","combined"),pch=c(15,15,NA),lty=c(NA,NA,5),col=c("gray80","gray30","black"),pt.cex=2,bty='n')}
    axis(1, at=1:8, labels=model_labels[[gg]])
    axis(2, at=seq(from=0,to=30,by=5), labels=seq(from=0,to=30,by=5), las=1)
    mtext(side=1, text="Nonstationary parameters", line=2.3, cex=0.85)
    mtext(side=2, text="# stations", line=2.3, cex=0.85)
    mtext(side=3, text=metric, cex=0.85)
    mtext(side=3, text=panel_labels[label_cnt], cex=0.85, adj=0); label_cnt <- label_cnt + 1
  }
}
dev.off()

##==============================================================================



##==============================================================================
## Compare model choices based
## on number of parameters and amount of available data
## NB: this version uses one covariate at a time
##============================

## individual versions with NLL, AIC and BIC
for (cc in names_covariates) {
  ##  comparison of model choice based on penalization of over-fitting
  best_models_datlen <- vector("list", 2); names(best_models_datlen) <- c("short", "long")
  best_models_all <- vector("list", 2); names(best_models_all) <- names_evm
  for (dl in 1:2) {
    best_models_datlen[[dl]] <- vector("list", 2); names(best_models_datlen[[dl]]) <- names_evm
    for (gg in names_evm) {
      best_models_datlen[[dl]][[gg]] <- vector("list", 3); names(best_models_datlen[[dl]][[gg]]) <- names(metrics)
      best_models_all[[gg]] <- vector("list", 3); names(best_models_all[[gg]]) <- names(metrics)
    }
  }
  model_choices <- c(1,2,3,4,5,6,7,8)
  model_labels <- vector('list', 2); names(model_labels) <- names_evm
  model_labels$gev <- c("none",expression(mu),expression(sigma),expression(xi),expression(mu * ', ' * sigma),expression(mu * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * mu * ', ' * sigma * ', ' * xi))
  model_labels$gpd <- c("none",expression(lambda),expression(sigma),expression(xi),expression(lambda * ', ' * sigma),expression(lambda * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * lambda * ', ' * sigma * ', ' * xi))
  model_breaks_short <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.001, to=8.001, by=1)))
  model_breaks_long <- sort(c(seq(from=0.999, to=7.999, by=1), seq(from=1.4, to=8.4, by=1)))
  model_breaks <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.4, to=8.4, by=1)))

  pdf(paste('../figures/model_choice_',cc,'.pdf',sep=''), width=6, height=6, pointsize=11, colormodel='cmyk')
  par(mfrow=c(3,2), mai=c(.5,.45,.2,.08))
  panel_labels <- c(expression('(a)'),expression('(b)'),expression('(c)'),expression('(d)'),expression('(e)'),expression('(f)')); label_cnt <- 1
  for (metric in names(metrics)) {
    for (gg in names_evm) {
      best_models_datlen$short[[gg]][[metric]] <- rep(NA, length(idx_short)); names(best_models_datlen$short[[gg]][[metric]]) <- site_names[idx_short]
      best_models_datlen$long[[gg]][[metric]] <- rep(NA, length(idx_long)); names(best_models_datlen$long[[gg]][[metric]]) <- site_names[idx_long]
      best_models_all[[gg]][[metric]] <- rep(NA, nsite); names(best_models_all[[gg]][[metric]]) <- site_names
      for (dd in site_names) {
        mat <- metrics[[metric]][[gg]][dd,cc,model_choices]
        idx_col <- which.min(mat)
        best_models_all[[gg]][[metric]][dd] <- idx_col
        if (dd %in% names(idx_short)) {best_models_datlen$short[[gg]][[metric]][dd] <- idx_col} else {best_models_datlen$long[[gg]][[metric]][dd] <- idx_col}
      }
      hist(best_models_all[[gg]][[metric]], breaks=model_breaks, freq=TRUE, xlim=c(0.5,8.5), lty=5, ylim=c(0,32), main="", xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i')
      grid()
      hist(best_models_datlen$short[[gg]][[metric]], breaks=model_breaks_short, freq=TRUE, col="white", xlim=c(0.5,8.5), ylim=c(0,17), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
      hist(best_models_datlen$long[[gg]][[metric]], breaks=model_breaks_long, freq=TRUE, col="gray50", xlim=c(0.5,8.5), ylim=c(0,17), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
      if (label_cnt==1) {legend(0.6,31,c("< 55 years","> 55 years","combined"),pch=c(15,15,NA),lty=c(NA,NA,5),col=c("gray80","gray30","black"),pt.cex=2,bty='n')}
      axis(1, at=1:8, labels=model_labels[[gg]])
      axis(2, at=seq(from=0,to=30,by=5), labels=seq(from=0,to=30,by=5), las=1)
      mtext(side=1, text="Nonstationary parameters", line=2.3, cex=0.85)
      mtext(side=2, text="# stations", line=2.3, cex=0.85)
      mtext(side=3, text=metric, cex=0.85)
      mtext(side=3, text=panel_labels[label_cnt], cex=0.85, adj=0); label_cnt <- label_cnt + 1
    }
  }
  dev.off()
}

##==============================================================================



##==============================================================================
## Compare model choices based
## on number of parameters and amount of available data
## NB: this version uses one covariate at a time
##============================

## main text version with only AIC, but all covariates and the pooled

model_choices <- c(1,2,3,4,5,6,7,8)
model_labels <- vector('list', 2); names(model_labels) <- names_evm
model_labels$gev <- c("none",expression(mu),expression(sigma),expression(xi),expression(mu * ', ' * sigma),expression(mu * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * mu * ', ' * sigma * ', ' * xi))
model_labels$gpd <- c("none",expression(lambda),expression(sigma),expression(xi),expression(lambda * ', ' * sigma),expression(lambda * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * lambda * ', ' * sigma * ', ' * xi))
model_breaks_short <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.001, to=8.001, by=1)))
model_breaks_long <- sort(c(seq(from=0.999, to=7.999, by=1), seq(from=1.4, to=8.4, by=1)))
model_breaks <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.4, to=8.4, by=1)))

pdf('../figures/model_choice_aic.pdf', width=6.5, height=9, pointsize=11, colormodel='cmyk')
par(mfrow=c(5,2), mai=c(.43,.65,.22,.1))

panel_labels <- c(expression('(a)'),expression('(b)'),expression('(c)'),expression('(d)'),
                  expression('(e)'),expression('(f)'),expression('(g)'),expression('(h)'),
                  expression('(i)'),expression('(j)')); label_cnt <- 1
metric <- "AIC"
# first, do the pooled one
best_models_datlen <- vector("list", 3); names(best_models_datlen) <- c("short", "long", "all")
for (dl in 1:3) {best_models_datlen[[dl]] <- vector("list", length(names_evm)); names(best_models_datlen[[dl]]) <- names_evm}
for (gg in names_evm) {
  best_models_datlen$short[[gg]] <- rep(NA, length(idx_short)); names(best_models_datlen$short[[gg]]) <- site_names[idx_short]
  best_models_datlen$long[[gg]] <- rep(NA, length(idx_long)); names(best_models_datlen$long[[gg]]) <- site_names[idx_long]
  best_models_datlen$all[[gg]] <- rep(NA, nsite); names(best_models_datlen$all[[gg]]) <- site_names
  for (dd in site_names) {
      mat <- metrics[[metric]][[gg]][dd,,model_choices]
      idx_min <- which.min(mat)
      idx_row <- idx_min%%nrow(mat)
      if(idx_row==0) {idx_row <- nrow(mat)}
      idx_col <- which.min(as.vector(mat[idx_row,]))
      best_models_datlen$all[[gg]][dd] <- idx_col
      if (dd %in% names(idx_short)) {best_models_datlen$short[[gg]][dd] <- idx_col} else {best_models_datlen$long[[gg]][dd] <- idx_col}
  }
  hist(best_models_datlen$all[[gg]], breaks=model_breaks, freq=TRUE, lty=5, xlim=c(0.5,8.5), ylim=c(0,18), main="", xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i')
  grid()
  hist(best_models_datlen$short[[gg]], breaks=model_breaks_short, freq=TRUE, col="white", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
  hist(best_models_datlen$long[[gg]], breaks=model_breaks_long, freq=TRUE, col="gray50", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
  if (label_cnt==1) {legend(0.6,18.5,c("< 55 years","> 55 years","combined"),pch=c(15,15,NA),lty=c(NA,NA,5),col=c("gray80","gray30","black"),pt.cex=2,bty='n')}
  axis(1, at=1:8, labels=model_labels[[gg]])
  axis(2, at=seq(from=0,to=30,by=5), labels=seq(from=0,to=30,by=5), las=1)
  mtext(side=1, text="Nonstationary parameters", line=2.3, cex=0.85)
  mtext(side=2, text="# stations", line=2.3, cex=0.85)
  if (gg=="gev") {mtext(side=2, text="All covariates", line=4, cex=0.85)}
  if(gg=="gev") {mtext(side=3, text="GEV", line=.6, cex=0.85)} else {mtext(side=3, text="GPD", line=.6, cex=0.85)}
  mtext(side=3, text=panel_labels[label_cnt], cex=0.85, adj=0); label_cnt <- label_cnt + 1
}
# then cycle through all the covariates
for (cc in names_covariates) {
  for (gg in names_evm) {
    best_models_datlen$short[[gg]] <- rep(NA, length(idx_short)); names(best_models_datlen$short[[gg]]) <- site_names[idx_short]
    best_models_datlen$long[[gg]] <- rep(NA, length(idx_long)); names(best_models_datlen$long[[gg]]) <- site_names[idx_long]
    best_models_datlen$all[[gg]] <- rep(NA, nsite); names(best_models_datlen$all[[gg]]) <- site_names
    for (dd in site_names) {
      mat <- metrics[[metric]][[gg]][dd,cc,model_choices]
      idx_col <- which.min(mat)
      best_models_datlen$all[[gg]][dd] <- idx_col
      if (dd %in% names(idx_short)) {best_models_datlen$short[[gg]][dd] <- idx_col} else {best_models_datlen$long[[gg]][dd] <- idx_col}
    }
    hist(best_models_datlen$all[[gg]], breaks=model_breaks, freq=TRUE, lty=5, xlim=c(0.5,8.5), ylim=c(0,18), main="", xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i')
    grid()
    hist(best_models_datlen$short[[gg]], breaks=model_breaks_short, freq=TRUE, col="white", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    hist(best_models_datlen$long[[gg]], breaks=model_breaks_long, freq=TRUE, col="gray50", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    if (label_cnt==1) {legend(0.6,18.5,c("< 55 years","> 55 years","combined"),pch=c(15,15,NA),lty=c(NA,NA,5),col=c("gray80","gray30","black"),pt.cex=2,bty='n')}
    axis(1, at=1:8, labels=model_labels[[gg]])
    axis(2, at=seq(from=0,to=30,by=5), labels=seq(from=0,to=30,by=5), las=1)
    mtext(side=1, text="Nonstationary parameters", line=2.3, cex=0.85)
    mtext(side=2, text="# stations", line=2.3, cex=0.85)
    if (gg=="gev") {mtext(side=2, text=covariate_labels[match(cc,names_covariates)], line=4, cex=0.85)}
    mtext(side=3, text=panel_labels[label_cnt], cex=0.85, adj=0); label_cnt <- label_cnt + 1
  }
}
dev.off()

##==============================================================================



##==============================================================================
## Compare return levels based
## on model choice
##============================

#> dim(rl$gev$FernandinaBeach$time)
#[1]  8  8 80
# 1st = return period (2/5/10/20/50/100/200/500) ; 2nd = model structure ; 3rd = year (2020:2099)

todo

##==============================================================================



##==============================================================================
## Map of covariate choice based on location
##==========================================

site_dat <- read.csv("../input_data/SiteInformation.csv", header=TRUE)
idx_map <- match(site_names, site_dat[,"site_name"])
site_dat <- site_dat[idx_map,]
lats <- site_dat[,"lat"]
lons <- site_dat[,"lon"]

best_models_map <- vector("list", 2); names(best_models_map) <- names_evm
for (gg in names_evm) {best_models_map[[gg]] <- vector("list", length(metrics)); names(best_models_map[[gg]]) <- names(metrics)}
model_choices <- c(1,2,3,4,5,6,7,8)
for (metric in names(metrics)) {
  for (gg in names_evm) {
    best_models_map[[gg]][[metric]] <- vector('list', 2); names(best_models_map[[gg]][[metric]]) <- c('covar','model')
    for (ii in 1:2) {best_models_map[[gg]][[metric]][[ii]] <- rep(NA, nsite); names(best_models_map[[gg]][[metric]][[ii]]) <- site_names}
    for (dd in site_names) {
        mat <- metrics[[metric]][[gg]][dd,,model_choices]
        idx_min <- which.min(mat)
        idx_row <- idx_min%%nrow(mat)
        if(idx_row==0) {idx_row <- nrow(mat)}
        idx_col <- which.min(as.vector(mat[idx_row,]))
        if(idx_col==1) {idx_row <- 5} # stationary model
        best_models_map[[gg]][[metric]]$covar[dd] <- idx_row
        best_models_map[[gg]][[metric]]$model[dd] <- idx_col
    }
  }
}

## site choices  #covariate_labels
metric <- "AIC" # which metric to use for the plot?
covariate_choices_site <- model_choices_site <- vector("list", 3)
names(covariate_choices_site) <- names(model_choices_site) <- c("gev","gpd","both")
for (gg in names_evm) {
  covariate_choices_site[[gg]] <- vector("list", 5)
  model_choices_site[[gg]] <- vector("list", length(model_choices))
  for (cc in 1:5) {covariate_choices_site[[gg]][[cc]] <- which(best_models_map[[gg]][[metric]]$covar==cc)}
  for (mm in 1:length(model_choices)) {model_choices_site[[gg]][[mm]] <- which(best_models_map[[gg]][[metric]]$model==mm)}
}
covariate_colors <- vector("list", 5); names(covariate_colors) <- covariate_labels
covariate_colors$Time <- "darkorange3"
covariate_colors$NAO <- "seagreen"
covariate_colors$Temp <- "mediumslateblue"
covariate_colors$`Sea level` <- "mediumvioletred"
covariate_colors$Stat <- "black"

## Make a plot

library(maps)

png("../figures/covariate_choice_map.png", width=500, height=760)
par(mfrow=c(2,1))

# GEV model structures
gg <- "gev"
map("world", fill=TRUE, col="gray85", bg="white", xlim=c(-105, -66.5), ylim=c(23.5, 47), mar=c(4,10,0,0))
map("state", fill=TRUE, col="gray85", bg="white", xlim=c(-105, -66.5), ylim=c(23.5, 47), mar=c(4,10,0,0), add=TRUE)
for (cc in 1:5) {points(lons[covariate_choices_site[[gg]][[cc]]], lats[covariate_choices_site[[gg]][[cc]]], col=covariate_colors[[cc]], cex=2, pch=16)}
mtext(side=1, text='Longitude', line=2.0, cex=1)
mtext(side=2, text='Latitude', line=3.8, cex=1)
mtext(side=3, text=expression('(a)'), line=0.1, cex=1.3, adj=0); mtext(side=3, text='GEV', line=0.1, cex=1.3)
axis(side=1, at=seq(from=-110, to=-64, by=5), labels=c("110 °W", "", "100 °W", "", "90 °W", "", "80 °W", "", "70 °W", ""), cex.axis=1)
axis(side=2, at=seq(from=20, to=50, by=5), labels=c("20 °N", "25 °N", "30 °N", "35 °N", "40 °N", "45 °N", "50 °N"), las=1, cex.axis=1)
legend(-74, 32, c("Stat","Time","NAO","Temp","Sea level"), pch=16,
       col=c(covariate_colors$Stat,covariate_colors$Time,covariate_colors$NAO,covariate_colors$Temp,covariate_colors$`Sea level`),
       pt.cex=2, cex=1, bty='n')

# GPD model structures
gg <- "gpd"
map("world", fill=TRUE, col="gray85", bg="white", xlim=c(-105, -66.5), ylim=c(23.5, 47), mar=c(4,10,0,0))
map("state", fill=TRUE, col="gray85", bg="white", xlim=c(-105, -66.5), ylim=c(23.5, 47), mar=c(4,10,0,0), add=TRUE)
for (cc in 1:5) {points(lons[covariate_choices_site[[gg]][[cc]]], lats[covariate_choices_site[[gg]][[cc]]], col=covariate_colors[[cc]], cex=2, pch=16)}
mtext(side=1, text='Longitude', line=2.0, cex=1)
mtext(side=2, text='Latitude', line=3.8, cex=1)
mtext(side=3, text=expression('(b)'), line=0.1, cex=1.3, adj=0); mtext(side=3, text='GPD', line=0.1, cex=1.3)
axis(side=1, at=seq(from=-110, to=-64, by=5), labels=c("110 °W", "", "100 °W", "", "90 °W", "", "80 °W", "", "70 °W", ""), cex.axis=1)
axis(side=2, at=seq(from=20, to=50, by=5), labels=c("20 °N", "25 °N", "30 °N", "35 °N", "40 °N", "45 °N", "50 °N"), las=1, cex.axis=1)

dev.off()

## for the text, how many of each covariate choice are there?
for (gg in names_evm) {for (cc in 1:5) print(paste(gg,covariate_labels[cc],length(covariate_choices_site[[gg]][[cc]])))}
##==============================================================================



##==============================================================================
## Map of new return period at points in future,
## based on "best" model choice (covariate and structure)
## for each site, for each of GEV and GPD
##======================================

# HERE NOW



# revisit model choices, but eliminate all with xi nonstationary
best_models_map <- vector("list", 2); names(best_models_map) <- names_evm
for (gg in names_evm) {best_models_map[[gg]] <- vector("list", length(metrics)); names(best_models_map[[gg]]) <- names(metrics)}
model_choices <- c(1,2,3,5) # eliminate any xi-nonstationary choices
for (metric in names(metrics)) {
  for (gg in names_evm) {
    best_models_map[[gg]][[metric]] <- vector('list', 2); names(best_models_map[[gg]][[metric]]) <- c('covar','model')
    for (ii in 1:2) {best_models_map[[gg]][[metric]][[ii]] <- rep(NA, nsite); names(best_models_map[[gg]][[metric]][[ii]]) <- site_names}
    for (dd in site_names) {
        mat <- metrics[[metric]][[gg]][dd,,model_choices]
        idx_min <- which.min(mat)
        idx_row <- idx_min%%nrow(mat)
        if(idx_row==0) {idx_row <- nrow(mat)}
        idx_col <- which.min(as.vector(mat[idx_row,]))
        if(idx_col==1) {idx_row <- 5} # stationary model
        best_models_map[[gg]][[metric]]$covar[dd] <- idx_row
        best_models_map[[gg]][[metric]]$model[dd] <- model_choices[idx_col]
    }
  }
}

metric <- "AIC"
rp_old <- "50"
year <- 2050
yy <- match(year, years_proj)

library(plot3D)

png("../figures/return_periods_map.png", width=500, height=760)
par(mfrow=c(2,1))

rp_new <- vector("list",length(names_evm)); names(rp_new) <- names_evm
# GEV model structures
gg <- "gev"
rp_new[[gg]] <- rep(NA, length(site_names)) # needs to be the return periods in the same order as lats and lons
for (dd in 1:length(site_names)) {
  cc <- best_models_map[[gg]][[metric]]$covar[dd]
  mm <- best_models_map[[gg]][[metric]]$model[dd]
  if (cc==5) {rp_new[[gg]][dd] <- as.numeric(rp_old)
  } else     {rp_new[[gg]][dd] <- mean(rp[[gg]][[dd]][[cc]][rp_old,mm,(yy-10):(yy+10)])}
}
idx_long2 <- idx_long[-8] # because New London is all weird
gr <- .bincode(rp_new[[gg]][idx_long2], seq(min(rp_new[[gg]][idx_long2]), max(rp_new[[gg]][idx_long2]), len=length(idx_long2)), include.lowest = T)
cols <- c("firebrick", "orange", "yellow", "white")
col <- colorRampPalette(cols, bias=1.75)(length(idx_long))[gr]
map("world", fill=TRUE, col="gray85", bg="white", xlim=c(-105, -66.5), ylim=c(23.5, 47), mar=c(4,10,0,0))
map("state", fill=TRUE, col="gray85", bg="white", xlim=c(-105, -66.5), ylim=c(23.5, 47), mar=c(4,10,0,0), add=TRUE)
points(lons[idx_long2], lats[idx_long2], bg=col, col="black", cex=2, pch=21)
mtext(side=1, text='Longitude', line=2.0, cex=1)
mtext(side=2, text='Latitude', line=3.8, cex=1)
mtext(side=3, text=expression('(a)'), line=0.1, cex=1.3, adj=0); mtext(side=3, text='GEV', line=0.1, cex=1.3)
axis(side=1, at=seq(from=-110, to=-64, by=5), labels=c("110 °W", "", "100 °W", "", "90 °W", "", "80 °W", "", "70 °W", ""), cex.axis=1)
axis(side=2, at=seq(from=20, to=50, by=5), labels=c("20 °N", "25 °N", "30 °N", "35 °N", "40 °N", "45 °N", "50 °N"), las=1, cex.axis=1)
#legend(-74, 32, c("Stat","Time","NAO","Temp","Sea level"), pch=16,
#       col=c(covariate_colors$Stat,covariate_colors$Time,covariate_colors$NAO,covariate_colors$Temp,covariate_colors$`Sea level`),
#       pt.cex=2, cex=1, bty='n')
colkey(clim=c(quantile(rp_new[[gg]][idx_long2],c(0,1))), at=seq(from=0,to=as.numeric(rp_old),by=10), labels=seq(from=1, to=as.numeric(rp_old), by=10), cex.axis=1.5, side=4, dist=0, width=0.65, add=TRUE, col=colorRampPalette(cols, bias=1.75)(length(idx_long2)))

# COLOR BAR CLEARLY NOT LINED UP WITH THE FIGURE DOTS

# GPD model structures
gg <- "gpd"
rp_new[[gg]] <- rep(NA, length(site_names)) # needs to be the return periods in the same order as lats and lons
for (dd in 1:length(site_names)) {
  cc <- best_models_map[[gg]][[metric]]$covar[dd]
  mm <- best_models_map[[gg]][[metric]]$model[dd]
  if (cc==5) {rp_new[[gg]][dd] <- as.numeric(rp_old)
  } else     {rp_new[[gg]][dd] <- mean(rp[[gg]][[dd]][[cc]][rp_old,mm,(yy-10):(yy+10)])}
}
gr <- .bincode(rp_new[[gg]][idx_long2], seq(min(rp_new[[gg]][idx_long2]), max(rp_new[[gg]][idx_long2]), len=length(idx_long2)), include.lowest = T)
cols <- c("firebrick", "orange", "yellow", "white")
col <- colorRampPalette(cols, bias=1.75)(length(idx_long))[gr]
map("world", fill=TRUE, col="gray85", bg="white", xlim=c(-105, -66.5), ylim=c(23.5, 47), mar=c(4,10,0,0))
map("state", fill=TRUE, col="gray85", bg="white", xlim=c(-105, -66.5), ylim=c(23.5, 47), mar=c(4,10,0,0), add=TRUE)
points(lons[idx_long2], lats[idx_long2], bg=col, col="black", cex=2, pch=21)
mtext(side=1, text='Longitude', line=2.0, cex=1)
mtext(side=2, text='Latitude', line=3.8, cex=1)
mtext(side=3, text=expression('(b)'), line=0.1, cex=1.3, adj=0); mtext(side=3, text='GPD', line=0.1, cex=1.3)
axis(side=1, at=seq(from=-110, to=-64, by=5), labels=c("110 °W", "", "100 °W", "", "90 °W", "", "80 °W", "", "70 °W", ""), cex.axis=1)
axis(side=2, at=seq(from=20, to=50, by=5), labels=c("20 °N", "25 °N", "30 °N", "35 °N", "40 °N", "45 °N", "50 °N"), las=1, cex.axis=1)
colkey(clim=c(quantile(rp_new[[gg]][idx_long2],c(0,1))), at=seq(from=0,to=as.numeric(rp_old),by=10), labels=seq(from=1, to=as.numeric(rp_old), by=10), cex.axis=1.5, side=4, dist=0, width=0.65, add=TRUE, col=colorRampPalette(cols, bias=1.75)(length(idx_long2)))

dev.off()


#cbind(site_names, lats, lons, rp_new$gpd, covariate_labels[best_models_map[[gg]][[metric]]$covar], best_models_map[[gg]][[metric]]$model, num_years)[idx_long,]


# OLD CODE
gr <- .bincode(rp_new, seq(min(rp_new), max(rp_new), len=length(rp_new)), include.lowest = T)
cols <- c("paleturquoise4", "lightgoldenrod", "firebrick")
bias <- 4.5
col <- colorRampPalette(cols, bias=bias)(length(rp_new))[gr]

png("./figures/Netherlands_cost_rcp85p50.png", width = 600, height = 600)
par(mfrow=c(1,1), mai=c(6,6,.5,.5))
map("world", fill=TRUE, col="gray85", bg="white", xlim=c(3, 7), ylim=c(51, 54), mar=c(8,7,0,0))
points(cost[[site]][[rcp]]$p50[,"lon"], cost[[site]][[rcp]]$p50[,"lat"], col=col, cex=2, pch=16)
points(rot[1], rot[2], col="black", pch=15, lwd=3, cex=1.5); text(rot[1]+0.5, rot[2], pos=1, "Rotterdam", cex=1.3)
points(ams[1], ams[2], col="black", pch=8, lwd=3, cex=1.5); text(ams[1]+0.5, ams[2]-0.12, pos=1, "Amsterdam", cex=1.3)
points(hag[1], hag[2], col="black", pch=15, lwd=3, cex=1.5); text(hag[1]-0.62, hag[2], pos=3, "The Hague", cex=1.3)
points(del[1], del[2], col="black", pch=15, lwd=3, cex=1.5); text(del[1]-0.25, del[2], pos=1, "Delfzijl", cex=1.3)
mtext(side=1, text='Longitude', line=2.5, cex=1.5)
mtext(side=2, text='Latitude', line=4.5, cex=1.5)
mtext(side=3, text='[$B]', line=1, cex=1.5, adj=1.13)
axis(side=1, at=seq(from=3, to=7, by=0.5), labels=c("3 E", "", "4 E", "", "5 E", "", "6 E", "", "7 E"), cex.axis=1.5)
axis(side=2, at=seq(from=51, to=54, by=0.5), labels=c("51 N", "", "52 N", "", "53 N", "", "54 N"), las=1, cex.axis=1.5)
colkey(clim=c(0,200), at=seq(from=0,to=200,by=10), labels=c("0","","20","","40","","60","","80","","100","","120","","140","","160","","180","","200"), cex.axis=1.5, side=4, dist=0, width=0.65, add=TRUE, col=colorRampPalette(cols, bias=bias)(length(cost_norm)))
dev.off()

##==============================================================================



##==============================================================================
## Map of model choice based on location
## (all covariates together? or just one?)
##======================================

# todo?

##==============================================================================



##==============================================================================
## End
##==============================================================================
