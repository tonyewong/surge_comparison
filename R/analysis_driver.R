##==============================================================================
## analysis_driver.R
##
## To be run after calibration_driver.R. This routine will make the surge tide
## projections, compute the model selection metrics and generate plots and
## numbers to report in the manuscript.
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


rm(list=ls())

library(extRemes)
library(Bolstad)
library(plot3D)
library(maps)

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

# yields parameters_maxlike[[model]][site,parameter] and parameters_maxpost...
source("best_models.R")

# yields return level and return period estimates
# ... maximum likelihood:
parameters <- parameters_maxlike
source("make_projections.R")
rp_maxlike <- rp
rl_maxlike <- rl
# ... maximum a posteriori
parameters <- parameters_maxpost
source("make_projections.R")
rp_maxpost <- rp
rl_maxpost <- rl
# Save return levels and periods in RData file
filename.returnlevels <- paste('../output/returnlevels_',today,'.RData', sep='')
print(paste("Saving return levels data to file",filename.returnlevels))
save(list=c('rl_maxlike','rl_maxpost', 'rp_maxlike', 'rp_maxpost', 'site_names', 'covariates', 'return_periods', 'years_proj',
            'nsite', 'nmodel', 'nrp', 'nyear', 'gev_models', 'gpd_models', 'data_calib'), file=filename.returnlevels)

metrics <- vector("list", 4); names(metrics) <- c("NLL","AIC","BIC","NPS")
metrics$NLL <- nll
metrics$AIC <- aic
metrics$BIC <- bic
metrics$NPS <- nps

## number of years of data for each site
num_years <- rep(NA, nsite); names(num_years) <- site_names
for (dd in 1:nsite) {num_years[dd] <- nrow(data_calib$gev[[dd]])}
idx_short <- which(num_years <= median(num_years)) # note: none of the stations (as of this writing) have num_years = median(num_years) (55 years)
idx_long <- which(num_years >= median(num_years))

## the most panel labels you might need - one for each of the 36 sites
panel_labels <- c(expression(bold('a')),expression(bold('b')),expression(bold('c')),expression(bold('d')),
                  expression(bold('e')),expression(bold('f')),expression(bold('g')),expression(bold('h')),
                  expression(bold('i')),expression(bold('j')),expression(bold('k')),expression(bold('l')),
                  expression(bold('m')),expression(bold('n')),expression(bold('o')),expression(bold('p')),
                  expression(bold('q')),expression(bold('r')),expression(bold('s')),expression(bold('t')),
                  expression(bold('u')),expression(bold('v')),expression(bold('w')),expression(bold('x')),
                  expression(bold('y')),expression(bold('z')),expression(bold('aa')),expression(bold('bb')),
                  expression(bold('cc')),expression(bold('dd')),expression(bold('ee')),expression(bold('ff')),
                  expression(bold('gg')),expression(bold('hh')),expression(bold('ii')),expression(bold('jj')))


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
    best_models_datlen[[dl]][[gg]] <- vector("list", length(metrics)); names(best_models_datlen[[dl]][[gg]]) <- names(metrics)
    best_models_all[[gg]] <- vector("list", length(metrics)); names(best_models_all[[gg]]) <- names(metrics) # yes, doing this multiple times
  }
}
model_choices <- c(1,2,3,4,5,6,7,8)
covariate_labels <- c("Time","Temp","Sea level","NAO","Stat") # needs to be in the same order as the second dimension of the `metrics` array: [1] "time" "temp" "sealevel" "nao"
covariate_breaks <- sort(c(seq(from=0.75, to=4.75, by=1), seq(from=1.25, to=5.25, by=1)))
covariate_breaks_short <- sort(c(seq(from=0.75, to=4.75, by=1), seq(from=1.001, to=5.001, by=1)))
covariate_breaks_long <- sort(c(seq(from=0.999, to=4.999, by=1), seq(from=1.25, to=5.25, by=1)))

pdf('../figures/covariate_choice_v.pdf', width=5.5, height=8, pointsize=11, colormodel='cmyk')
par(mfrow=c(4,2), mai=c(.5,.45,.2,.08))
label_cnt <- 1
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

h_offset <- 0.025
pdf('../figures/covariate_choice_h.pdf', width=7.5, height=3.5, pointsize=11, colormodel='cmyk')
par(mfrow=c(2,4))
label_cnt <- 1
for (gg in names_evm) {
  for (metric in names(metrics)) {
    if (label_cnt==1) {mai_panel <- c(.36, .52, .2, 0)  # .45+.08 = .53 total width; .5+.2 = .7 total height
    } else if (label_cnt==5) {mai_panel <- c(.5, .52, .08, 0)
    } else {mai_panel[2] <- mai_panel[2] - h_offset; mai_panel[4] <- mai_panel[4] + h_offset}
    par(mai=mai_panel)
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
    axis(1, at=c(1,3,5), label=rep("",3), tck = -0.03); axis(1, at=c(1,3,5), labels=covariate_labels[c(1,3,5)], line=-0.7, lwd=0, cex.axis=0.9)
    axis(1, at=c(2,4), label=rep("",2), tck = -0.12); axis(1, at=c(2,4), labels=covariate_labels[c(2,4)], line=0.1, lwd=0, cex.axis=0.9)
    axis(2, at=seq(from=0,to=30,by=5), labels=rep("",length(seq(from=0,to=30,by=5))), tck=-0.03)
    axis(2, at=seq(from=0,to=30,by=5), labels=seq(from=0,to=30,by=5), line=-0.5, las=1, lwd=0)
    if (label_cnt>=5) mtext(side=1, text="Covariate", line=2.5, cex=0.85)
    if (label_cnt==1 | label_cnt==5) mtext(side=2, text="# stations", line=1.7, cex=0.85)
    if (label_cnt==1) mtext(side=2, text="GEV", line=3.1, cex=0.85)
    if (label_cnt==5) mtext(side=2, text="GPD", line=3.1, cex=0.85)
    if (label_cnt<=4) mtext(side=3, text=metric, cex=0.85)
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
    best_models_datlen[[dl]][[gg]] <- vector("list", length(metrics)); names(best_models_datlen[[dl]][[gg]]) <- names(metrics)
    best_models_all[[gg]] <- vector("list", length(metrics)); names(best_models_all[[gg]]) <- names(metrics)
  }
}
model_choices <- c(1,2,3,4,5,6,7,8)
model_labels <- vector('list', 2); names(model_labels) <- names_evm
model_labels$gev <- c("none",expression(mu),expression(sigma),expression(xi),expression(mu * ', ' * sigma),expression(mu * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * mu * ', ' * sigma * ', ' * xi))
model_labels$gpd <- c("none",expression(lambda),expression(sigma),expression(xi),expression(lambda * ', ' * sigma),expression(lambda * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * lambda * ', ' * sigma * ', ' * xi))
model_breaks_short <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.001, to=8.001, by=1)))
model_breaks_long <- sort(c(seq(from=0.999, to=7.999, by=1), seq(from=1.4, to=8.4, by=1)))
model_breaks <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.4, to=8.4, by=1)))

pdf('../figures/model_choice.pdf', width=6, height=8, pointsize=11, colormodel='cmyk')
par(mfrow=c(4,2), mai=c(.5,.45,.2,.08))
label_cnt <- 1
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
    hist(best_models_datlen$short[[gg]][[metric]], breaks=model_breaks_short, freq=TRUE, col="gray80", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    hist(best_models_datlen$long[[gg]][[metric]], breaks=model_breaks_long, freq=TRUE, col="gray30", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
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

## individual versions with NLL, AIC, BIC and negative posterior score (NPS)
for (cc in names_covariates) {
  ##  comparison of model choice based on penalization of over-fitting
  best_models_datlen <- vector("list", 2); names(best_models_datlen) <- c("short", "long")
  best_models_all <- vector("list", 2); names(best_models_all) <- names_evm
  for (dl in 1:2) {
    best_models_datlen[[dl]] <- vector("list", 2); names(best_models_datlen[[dl]]) <- names_evm
    for (gg in names_evm) {
      best_models_datlen[[dl]][[gg]] <- vector("list", length(metrics)); names(best_models_datlen[[dl]][[gg]]) <- names(metrics)
      best_models_all[[gg]] <- vector("list", length(metrics)); names(best_models_all[[gg]]) <- names(metrics)
    }
  }
  model_choices <- c(1,2,3,4,5,6,7,8)
  model_labels <- vector('list', 2); names(model_labels) <- names_evm
  model_labels$gev <- c("none",expression(mu),expression(sigma),expression(xi),expression(mu * ', ' * sigma),expression(mu * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * mu * ', ' * sigma * ', ' * xi))
  model_labels$gpd <- c("none",expression(lambda),expression(sigma),expression(xi),expression(lambda * ', ' * sigma),expression(lambda * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * lambda * ', ' * sigma * ', ' * xi))
  model_breaks_short <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.001, to=8.001, by=1)))
  model_breaks_long <- sort(c(seq(from=0.999, to=7.999, by=1), seq(from=1.4, to=8.4, by=1)))
  model_breaks <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.4, to=8.4, by=1)))

  pdf(paste('../figures/model_choice_',cc,'.pdf',sep=''), width=6, height=8, pointsize=11, colormodel='cmyk')
  par(mfrow=c(4,2), mai=c(.5,.45,.2,.08))
  label_cnt <- 1
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
      hist(best_models_all[[gg]][[metric]], breaks=model_breaks, freq=TRUE, xlim=c(0.5,8.5), lty=5, ylim=c(0,35), main="", xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i')
      grid()
      hist(best_models_datlen$short[[gg]][[metric]], breaks=model_breaks_short, freq=TRUE, col="gray80", xlim=c(0.5,8.5), ylim=c(0,17), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
      hist(best_models_datlen$long[[gg]][[metric]], breaks=model_breaks_long, freq=TRUE, col="gray30", xlim=c(0.5,8.5), ylim=c(0,17), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
      if (label_cnt==1) {legend(0.6,34,c("< 55 years","> 55 years","combined"),pch=c(15,15,NA),lty=c(NA,NA,5),col=c("gray80","gray30","black"),pt.cex=2,bty='n')}
      axis(1, at=1:8, labels=model_labels[[gg]])
      axis(2, at=seq(from=0,to=35,by=5), labels=seq(from=0,to=35,by=5), las=1)
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
## NB: this version uses one covariate at a time - loops over them.
##============================

ymaxs <- c(30,20,33,35); names(ymaxs) <- names(metrics)

for (metric in names(metrics)) {
  model_choices <- c(1,2,3,4,5,6,7,8)
  model_labels <- vector('list', 2); names(model_labels) <- names_evm
  model_labels$gev <- c("none",expression(mu),expression(sigma),expression(xi),expression(mu * ', ' * sigma),expression(mu * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * mu * ', ' * sigma * ', ' * xi))
  model_labels$gpd <- c("none",expression(lambda),expression(sigma),expression(xi),expression(lambda * ', ' * sigma),expression(lambda * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * lambda * ', ' * sigma * ', ' * xi))
  model_breaks_short <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.001, to=8.001, by=1)))
  model_breaks_long <- sort(c(seq(from=0.999, to=7.999, by=1), seq(from=1.4, to=8.4, by=1)))
  model_breaks <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.4, to=8.4, by=1)))

  pdf(paste('../figures/model_choice_',metric,'.pdf', sep=''), width=6.5, height=9, pointsize=11, colormodel='cmyk')
  par(mfrow=c(5,2), mai=c(.43,.65,.22,.1))
  label_cnt <- 1

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
    hist(best_models_datlen$all[[gg]], breaks=model_breaks, freq=TRUE, lty=5, xlim=c(0.5,8.5), ylim=c(0,ymaxs[metric]), main="", xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i')
    grid()
    hist(best_models_datlen$short[[gg]], breaks=model_breaks_short, freq=TRUE, col="gray80", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    hist(best_models_datlen$long[[gg]], breaks=model_breaks_long, freq=TRUE, col="gray30", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    if (label_cnt==1) {legend(0.6,ymaxs[metric],c("< 55 years","> 55 years","combined"),pch=c(15,15,NA),lty=c(NA,NA,5),col=c("gray80","gray30","black"),pt.cex=2,bty='n')}
    axis(1, at=1:8, labels=model_labels[[gg]])
    axis(2, at=seq(from=0,to=ymaxs[metric],by=5), labels=seq(from=0,to=ymaxs[metric],by=5), las=1)
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
      hist(best_models_datlen$all[[gg]], breaks=model_breaks, freq=TRUE, lty=5, xlim=c(0.5,8.5), ylim=c(0,ymaxs[metric]), main="", xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i')
      grid()
      hist(best_models_datlen$short[[gg]], breaks=model_breaks_short, freq=TRUE, col="gray80", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
      hist(best_models_datlen$long[[gg]], breaks=model_breaks_long, freq=TRUE, col="gray30", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
      if (label_cnt==1) {legend(0.6,ymaxs[metric],c("< 55 years","> 55 years","combined"),pch=c(15,15,NA),lty=c(NA,NA,5),col=c("gray80","gray30","black"),pt.cex=2,bty='n')}
      axis(1, at=1:8, labels=model_labels[[gg]])
      axis(2, at=seq(from=0,to=ymaxs[metric],by=5), labels=seq(from=0,to=ymaxs[metric],by=5), las=1)
      mtext(side=1, text="Nonstationary parameters", line=2.3, cex=0.85)
      mtext(side=2, text="# stations", line=2.3, cex=0.85)
      if (gg=="gev") {mtext(side=2, text=covariate_labels[match(cc,names_covariates)], line=4, cex=0.85)}
      mtext(side=3, text=panel_labels[label_cnt], cex=0.85, adj=0); label_cnt <- label_cnt + 1
    }
  }
  dev.off()
}

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

for (metric in names(metrics)) {

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

  png(paste("../figures/covariate_choice_map_",metric,"_v.png", sep=""), width=500, height=760)
  par(mfrow=c(2,1))
  label_cnt <- 1
  # GEV model structures
  gg <- "gev"
  map("world", fill=TRUE, col="gray85", bg="white", xlim=c(-105, -66.5), ylim=c(23.5, 47), mar=c(4,10,0,0))
  map("state", fill=TRUE, col="gray85", bg="white", xlim=c(-105, -66.5), ylim=c(23.5, 47), mar=c(4,10,0,0), add=TRUE)
  for (cc in 1:5) {points(lons[covariate_choices_site[[gg]][[cc]]], lats[covariate_choices_site[[gg]][[cc]]], col=covariate_colors[[cc]], cex=2, pch=16)}
  mtext(side=1, text='Longitude', line=2.0, cex=1)
  mtext(side=2, text='Latitude', line=3.8, cex=1)
  mtext(side=3, text=panel_labels[label_cnt], line=0.1, cex=1.3, adj=0); label_cnt <- label_cnt + 1
  mtext(side=3, text='GEV', line=0.1, cex=1.3)
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
  mtext(side=3, text=panel_labels[label_cnt], line=0.1, cex=1.3, adj=0); label_cnt <- label_cnt + 1
  mtext(side=3, text='GPD', line=0.1, cex=1.3)
  axis(side=1, at=seq(from=-110, to=-64, by=5), labels=c("110 °W", "", "100 °W", "", "90 °W", "", "80 °W", "", "70 °W", ""), cex.axis=1)
  axis(side=2, at=seq(from=20, to=50, by=5), labels=c("20 °N", "25 °N", "30 °N", "35 °N", "40 °N", "45 °N", "50 °N"), las=1, cex.axis=1)
  dev.off()
  png(paste("../figures/covariate_choice_map_",metric,"_h.png", sep=""), width=10, height=6, units="in", res=300, pointsize=11)
  par(mfrow=c(1,2))
  label_cnt <- 1
  # GEV model structures
  gg <- "gev"
  map("world", fill=TRUE, col="gray85", bg="white", xlim=c(-103.5, -65.5), ylim=c(22.5, 47.5), mar=c(4,10,0,0))
  map("state", fill=TRUE, col="gray85", bg="white", xlim=c(-103.5, -65.5), ylim=c(22.5, 47.5), mar=c(4,10,0,0), add=TRUE)
  for (cc in 1:5) {points(lons[covariate_choices_site[[gg]][[cc]]], lats[covariate_choices_site[[gg]][[cc]]], col=covariate_colors[[cc]], cex=1.5, pch=16)}
  mtext(side=1, text='Longitude', line=2.0, cex=1)
  mtext(side=2, text='Latitude', line=3.8, cex=1)
  mtext(side=3, text=panel_labels[label_cnt], line=0.1, cex=1, adj=0); label_cnt <- label_cnt + 1
  mtext(side=3, text='GEV', line=0.1, cex=1)
  axis(side=1, at=seq(from=-110, to=-64, by=5), labels=c("110 °W", "", "100 °W", "", "90 °W", "", "80 °W", "", "70 °W", ""), cex.axis=1)
  axis(side=2, at=seq(from=20, to=50, by=5), labels=c("20 °N", "25 °N", "30 °N", "35 °N", "40 °N", "45 °N", "50 °N"), las=1, cex.axis=1)
  legend(-76, 35.5, c("Stat","Time","NAO","Temp","Sea level"), pch=16,
         col=c(covariate_colors$Stat,covariate_colors$Time,covariate_colors$NAO,covariate_colors$Temp,covariate_colors$`Sea level`),
         pt.cex=1.5, cex=1, bty='n')
  # GPD model structures
  gg <- "gpd"
  map("world", fill=TRUE, col="gray85", bg="white", xlim=c(-103.5, -65.5), ylim=c(22.5, 47.5), mar=c(4,10,0,0))
  map("state", fill=TRUE, col="gray85", bg="white", xlim=c(-103.5, -65.5), ylim=c(22.5, 47.5), mar=c(4,10,0,0), add=TRUE)
  for (cc in 1:5) {points(lons[covariate_choices_site[[gg]][[cc]]], lats[covariate_choices_site[[gg]][[cc]]], col=covariate_colors[[cc]], cex=1.5, pch=16)}
  mtext(side=1, text='Longitude', line=2.1, cex=1)
  mtext(side=2, text='Latitude', line=3.6, cex=1)
  mtext(side=3, text=panel_labels[label_cnt], line=0.1, cex=1, adj=0); label_cnt <- label_cnt + 1
  mtext(side=3, text='GPD', line=0.1, cex=1)
  axis(side=1, at=seq(from=-110, to=-64, by=5), labels=c("110 °W", "", "100 °W", "", "90 °W", "", "80 °W", "", "70 °W", ""), cex.axis=1)
  axis(side=2, at=seq(from=20, to=50, by=5), labels=c("20 °N", "25 °N", "30 °N", "35 °N", "40 °N", "45 °N", "50 °N"), las=1, cex.axis=1)
  dev.off()
}


## for the text, how many of each covariate choice are there?
for (gg in names_evm) {for (cc in 1:5) print(paste(gg,covariate_labels[cc],length(covariate_choices_site[[gg]][[cc]])))}
##==============================================================================



##==============================================================================
## Maps and plots of return periods and return levels
##==============================================================================

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

## calculated return periods and levels are the mean of a 21-year moving window
## average, centered on the year in question
metric <- "NPS"
rp_old <- "50"
present <- 2020
year <- 2050
either_side <- 10 # window width is 2*either_side + 1
yi <- match(year, years_proj)
rp_new <- rl_old <- rl_new <- rl_smoothed <- vector("list",length(names_evm))
names(rp_new) <- names(rl_old) <- names(rl_new) <- names(rl_smoothed) <- names_evm
for (gg in names_evm) {
  rp_new[[gg]] <- rl_old[[gg]] <- rl_new[[gg]] <- rep(NA, length(site_names)) # needs to be the return periods in the same order as lats and lons
  rl_smoothed[[gg]] <- vector('list', length(site_names)); names(rl_smoothed[[gg]]) <- site_names
  for (dd in 1:length(site_names)) {
    cc <- best_models_map[[gg]][[metric]]$covar[dd]
    mm <- best_models_map[[gg]][[metric]]$model[dd]
    rl_smoothed[[gg]][[dd]] <- vector("list", length(names_covariates)); names(rl_smoothed[[gg]][[dd]]) <- names_covariates
    if (cc==5) {
      # stationary model, no moving average needed here
      rl_old[[gg]][dd] <- rl[[gg]][[dd]][[1]][rp_old,mm,match(present,years_proj)] # doesn't matter which covariate you take since it's stationary
      rp_new[[gg]][dd] <- as.numeric(rp_old)
      rl_new[[gg]][dd] <- rl_old[[gg]][dd]
    } else {
      # nonstationary model, use moving average
      rl_old[[gg]][dd] <- rl[[gg]][[dd]][[cc]][rp_old,mm,match(present,years_proj)]
      rp_new[[gg]][dd] <- mean(rp[[gg]][[dd]][[cc]][rp_old,mm,(yi-either_side):(yi+either_side)])
      rl_new[[gg]][dd] <- mean(rl[[gg]][[dd]][[cc]][rp_old,mm,(yi-either_side):(yi+either_side)])
    }
    for (c in names_covariates) {
      rl_smoothed[[gg]][[dd]][[c]] <- mat.or.vec(nr=(year-present+1), nc=nmodel)
      for (yy in present:year) {
        yi2 <- match(yy, years_proj)
        for (m in 1:nmodel) {
          rl_smoothed[[gg]][[dd]][[c]][yy-present+1, m] <- mean(rl[[gg]][[dd]][[c]][rp_old,m,(yi2-either_side):(yi2+either_side)])
        }
      }
    }
  }
}


##==============================================================================
## Line plot of return level projections for all covariates and all models, for
## a given return period, and highlight the projection based on the "best" model
## choice (covariate and structure). A panel for each site, and a version for
## each of GEV and GPD
##======================================

year <- 2050
yi <- match(year, years_proj)
yi_present <- match(present, years_proj)

for (metric in names(metrics)) {
  for (gg in names_evm) {
    label_cnt <- 1
    pdf(paste("../figures/projections_long_rl",rp_old,"_",gg,"_",metric,".pdf", sep=""), width=7.5, height=9, pointsize=11, colormodel='cmyk')
    par(mfrow=c(5,4), mai=c(.43,.45,.22,.15))
    for (dd in idx_long) {
      ymin <- min(rl_smoothed[[gg]][[dd]][[1]], na.rm=TRUE); for (c in names_covariates) {if (min(rl_smoothed[[gg]][[dd]][[c]], na.rm=TRUE) < ymin) ymin <- min(rl_smoothed[[gg]][[dd]][[c]], na.rm=TRUE)}
      ymax <- max(rl_smoothed[[gg]][[dd]][[1]], na.rm=TRUE); for (c in names_covariates) {if (max(rl_smoothed[[gg]][[dd]][[c]], na.rm=TRUE) > ymax) ymax <- max(rl_smoothed[[gg]][[dd]][[c]], na.rm=TRUE)}
      cc <- best_models_map[[gg]][[metric]]$covar[dd]
      mm <- best_models_map[[gg]][[metric]]$model[dd]
      if (cc==5) {cc <- 1; mm <- 1}
      # actual plots
      plot(years_proj[yi_present:yi], rl_smoothed[[gg]][[dd]][[cc]][,mm], lwd=2, type='l', col="black", xaxs='i', xlim=c(2020,year), ylim=c(ymin,ymax), main="", xaxt='n', xlab='', ylab='')
      grid()
      for (m in 1:nmodel) {for (c in names_covariates) {lines(years_proj[yi_present:yi], rl_smoothed[[gg]][[dd]][[c]][,m], lwd=0.5, col=covariate_colors[[match(c,names_covariates)]])}}
      lines(years_proj[yi_present:yi], rl_smoothed[[gg]][[dd]][[cc]][,mm], lwd=2.5, col="black")
      axis(1, at=seq(from=2020, to=year, by=5), labels=rep("",length(seq(from=2020, to=year, by=5))))
      axis(1, at=seq(from=2020, to=year, by=10), labels=seq(from=2020, to=year, by=10))
      #axis(2, at=seq(from=ymin,to=ymax,by=500), labels=seq(from=ymin,to=ymax,by=500), las=1)
      mtext(side=1, text="Year", line=2.3, cex=0.85)
      mtext(side=2, text="Height [mm]", line=2.3, cex=0.85)
      mtext(side=3, text=panel_labels[label_cnt], cex=0.85, adj=-0.075); label_cnt <- label_cnt + 1
      mtext(side=3, text=site_dat[dd, "formal"], cex=0.85)
    }
    plot(-1, -1, xlim=c(0,2), ylim=c(0,2), col="white", xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
    legend(0, 2, c("best model", covariate_labels[1:4]), lwd=c(2,1.3,1.3,1.3,1.3), col=c("black", unlist(covariate_colors)), cex=1.3, bty='n')
    dev.off()
  }
}


# same, but for all the sites

for (metric in names(metrics)) {
  for (gg in names_evm) {
    label_cnt <- 1
    pdf(paste("../figures/projections_rl",rp_old,"_",gg,"_",metric,".pdf", sep=""), width=12, height=13, pointsize=11, colormodel='cmyk')
    par(mfrow=c(7,6), mai=c(.43,.45,.22,.15))
    for (dd in 1:length(site_names)) {
      ymin <- min(rl_smoothed[[gg]][[dd]][[1]], na.rm=TRUE); for (c in names_covariates) {if (min(rl_smoothed[[gg]][[dd]][[c]], na.rm=TRUE) < ymin) ymin <- min(rl_smoothed[[gg]][[dd]][[c]], na.rm=TRUE)}
      ymax <- max(rl_smoothed[[gg]][[dd]][[1]], na.rm=TRUE); for (c in names_covariates) {if (max(rl_smoothed[[gg]][[dd]][[c]], na.rm=TRUE) > ymax) ymax <- max(rl_smoothed[[gg]][[dd]][[c]], na.rm=TRUE)}
      cc <- best_models_map[[gg]][[metric]]$covar[dd]
      mm <- best_models_map[[gg]][[metric]]$model[dd]
      if (cc==5) {cc <- 1; mm <- 1}
      # actual plots
      plot(years_proj[yi_present:yi], rl_smoothed[[gg]][[dd]][[cc]][,mm], lwd=2, type='l', col="black", xaxs='i', xlim=c(2020,year), ylim=c(ymin,ymax), main="", xaxt='n', xlab='', ylab='')
      grid()
      for (m in 1:nmodel) {for (c in names_covariates) {lines(years_proj[yi_present:yi], rl_smoothed[[gg]][[dd]][[c]][,m], lwd=0.5, col=covariate_colors[[match(c,names_covariates)]])}}
      lines(years_proj[yi_present:yi], rl_smoothed[[gg]][[dd]][[cc]][,mm], lwd=2.5, col="black")
      axis(1, at=seq(from=2020, to=year, by=5), labels=rep("",length(seq(from=2020, to=year, by=5))))
      axis(1, at=seq(from=2020, to=year, by=10), labels=seq(from=2020, to=year, by=10))
      #axis(2, at=seq(from=ymin,to=ymax,by=500), labels=seq(from=ymin,to=ymax,by=500), las=1)
      mtext(side=1, text="Year", line=2.3, cex=0.85)
      mtext(side=2, text="Height [mm]", line=2.3, cex=0.85)
      mtext(side=3, text=panel_labels[label_cnt], cex=0.85, adj=-0.075); label_cnt <- label_cnt + 1
      mtext(side=3, text=site_dat[dd, "formal"], cex=0.85)
    }
    plot(-1, -1, xlim=c(0,2), ylim=c(0,2), col="white", xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
    legend(0, 2, c("best model", covariate_labels[1:4]), lwd=c(2,1.3,1.3,1.3,1.3), col=c("black", unlist(covariate_colors)), cex=1.3, bty='n')
    dev.off()
  }
}


## Numbers for paper:
# Key West
metric <- "NPS"; gg <- "gpd"; dd <- "KeyWestFL"
cc <- best_models_map[[gg]][[metric]]$covar[dd]; mm <- best_models_map[[gg]][[metric]]$model[dd]
print(paste(cc, mm)); print(rl_smoothed[[gg]][[dd]][[cc]][year-present+1,mm]); print(diff(range(rl_smoothed[[gg]][[dd]][[cc]][year-present+1,])))
# Boston
metric <- "NPS"; gg <- "gpd"; dd <- "BostonMA"
cc <- best_models_map[[gg]][[metric]]$covar[dd]; mm <- best_models_map[[gg]][[metric]]$model[dd]
print(paste(cc, mm)); print(rl_smoothed[[gg]][[dd]][[cc]][year-present+1,mm]); print(diff(range(rl_smoothed[[gg]][[dd]][[cc]][year-present+1,])))


## calculated return periods and levels are the mean of a 21-year moving window
## average, centered on the year in question
metric <- "NPS"
rp_old <- "50"
present <- 2020
year <- 2050
either_side <- 10 # window width is 2*either_side + 1
yi <- match(year, years_proj)
rp_new <- rl_old <- rl_new <- rl_smoothed <- vector("list",length(names_evm))
names(rp_new) <- names(rl_old) <- names(rl_new) <- names(rl_smoothed) <- names_evm
for (gg in names_evm) {
  rp_new[[gg]] <- rl_old[[gg]] <- rl_new[[gg]] <- rep(NA, length(site_names)) # needs to be the return periods in the same order as lats and lons
  rl_smoothed[[gg]] <- vector('list', length(site_names)); names(rl_smoothed[[gg]]) <- site_names
  for (dd in 1:length(site_names)) {
    cc <- best_models_map[[gg]][[metric]]$covar[dd]
    mm <- best_models_map[[gg]][[metric]]$model[dd]
    rl_smoothed[[gg]][[dd]] <- vector("list", length(names_covariates)); names(rl_smoothed[[gg]][[dd]]) <- names_covariates
    if (cc==5) {
      # stationary model, no moving average needed here
      rl_old[[gg]][dd] <- rl[[gg]][[dd]][[1]][rp_old,mm,match(present,years_proj)] # doesn't matter which covariate you take since it's stationary
      rp_new[[gg]][dd] <- as.numeric(rp_old)
      rl_new[[gg]][dd] <- rl_old[[gg]][dd]
    } else {
      # nonstationary model, use moving average
      rl_old[[gg]][dd] <- rl[[gg]][[dd]][[cc]][rp_old,mm,match(present,years_proj)]
      rp_new[[gg]][dd] <- mean(rp[[gg]][[dd]][[cc]][rp_old,mm,(yi-either_side):(yi+either_side)])
      rl_new[[gg]][dd] <- mean(rl[[gg]][[dd]][[cc]][rp_old,mm,(yi-either_side):(yi+either_side)])
    }
    for (c in names_covariates) {
      rl_smoothed[[gg]][[dd]][[c]] <- mat.or.vec(nr=(year-present+1), nc=nmodel)
      for (yy in present:year) {
        yi2 <- match(yy, years_proj)
        for (m in 1:nmodel) {
          rl_smoothed[[gg]][[dd]][[c]][yy-present+1, m] <- mean(rl[[gg]][[dd]][[c]][rp_old,m,(yi2-either_side):(yi2+either_side)])
        }
      }
    }
  }
}


##======================================
## Same, but for 100-year return period

## calculated return periods and levels are the mean of a 21-year moving window
## average, centered on the year in question
metric <- "NPS"
rp_old <- "100"
present <- 2020
year <- 2050
either_side <- 10 # window width is 2*either_side + 1
yi <- match(year, years_proj)
rp_new <- rl_old <- rl_new <- rl_smoothed <- vector("list",length(names_evm))
names(rp_new) <- names(rl_old) <- names(rl_new) <- names(rl_smoothed) <- names_evm
for (gg in names_evm) {
  rp_new[[gg]] <- rl_old[[gg]] <- rl_new[[gg]] <- rep(NA, length(site_names)) # needs to be the return periods in the same order as lats and lons
  rl_smoothed[[gg]] <- vector('list', length(site_names)); names(rl_smoothed[[gg]]) <- site_names
  for (dd in 1:length(site_names)) {
    cc <- best_models_map[[gg]][[metric]]$covar[dd]
    mm <- best_models_map[[gg]][[metric]]$model[dd]
    rl_smoothed[[gg]][[dd]] <- vector("list", length(names_covariates)); names(rl_smoothed[[gg]][[dd]]) <- names_covariates
    if (cc==5) {
      # stationary model, no moving average needed here
      rl_old[[gg]][dd] <- rl[[gg]][[dd]][[1]][rp_old,mm,match(present,years_proj)] # doesn't matter which covariate you take since it's stationary
      rp_new[[gg]][dd] <- as.numeric(rp_old)
      rl_new[[gg]][dd] <- rl_old[[gg]][dd]
    } else {
      # nonstationary model, use moving average
      rl_old[[gg]][dd] <- rl[[gg]][[dd]][[cc]][rp_old,mm,match(present,years_proj)]
      rp_new[[gg]][dd] <- mean(rp[[gg]][[dd]][[cc]][rp_old,mm,(yi-either_side):(yi+either_side)])
      rl_new[[gg]][dd] <- mean(rl[[gg]][[dd]][[cc]][rp_old,mm,(yi-either_side):(yi+either_side)])
    }
    for (c in names_covariates) {
      rl_smoothed[[gg]][[dd]][[c]] <- mat.or.vec(nr=(year-present+1), nc=nmodel)
      for (yy in present:year) {
        yi2 <- match(yy, years_proj)
        for (m in 1:nmodel) {
          rl_smoothed[[gg]][[dd]][[c]][yy-present+1, m] <- mean(rl[[gg]][[dd]][[c]][rp_old,m,(yi2-either_side):(yi2+either_side)])
        }
      }
    }
  }
}

year <- 2050
yi <- match(year, years_proj)
yi_present <- match(present, years_proj)

for (metric in names(metrics)) {
  for (gg in names_evm) {
    label_cnt <- 1
    pdf(paste("../figures/projections_long_rl",rp_old,"_",gg,"_",metric,".pdf", sep=""), width=7.5, height=9, pointsize=11, colormodel='cmyk')
    par(mfrow=c(5,4), mai=c(.43,.45,.22,.15))
    for (dd in idx_long) {
      ymin <- min(rl_smoothed[[gg]][[dd]][[1]], na.rm=TRUE); for (c in names_covariates) {if (min(rl_smoothed[[gg]][[dd]][[c]], na.rm=TRUE) < ymin) ymin <- min(rl_smoothed[[gg]][[dd]][[c]], na.rm=TRUE)}
      ymax <- max(rl_smoothed[[gg]][[dd]][[1]], na.rm=TRUE); for (c in names_covariates) {if (max(rl_smoothed[[gg]][[dd]][[c]], na.rm=TRUE) > ymax) ymax <- max(rl_smoothed[[gg]][[dd]][[c]], na.rm=TRUE)}
      cc <- best_models_map[[gg]][[metric]]$covar[dd]
      mm <- best_models_map[[gg]][[metric]]$model[dd]
      if (cc==5) {cc <- 1; mm <- 1}
      # actual plots
      plot(years_proj[yi_present:yi], rl_smoothed[[gg]][[dd]][[cc]][,mm], lwd=2, type='l', col="black", xaxs='i', xlim=c(2020,year), ylim=c(ymin,ymax), main="", xaxt='n', xlab='', ylab='')
      grid()
      for (m in 1:nmodel) {for (c in names_covariates) {lines(years_proj[yi_present:yi], rl_smoothed[[gg]][[dd]][[c]][,m], lwd=0.5, col=covariate_colors[[match(c,names_covariates)]])}}
      lines(years_proj[yi_present:yi], rl_smoothed[[gg]][[dd]][[cc]][,mm], lwd=2.5, col="black")
      axis(1, at=seq(from=2020, to=year, by=5), labels=rep("",length(seq(from=2020, to=year, by=5))))
      axis(1, at=seq(from=2020, to=year, by=10), labels=seq(from=2020, to=year, by=10))
      #axis(2, at=seq(from=ymin,to=ymax,by=500), labels=seq(from=ymin,to=ymax,by=500), las=1)
      mtext(side=1, text="Year", line=2.3, cex=0.85)
      mtext(side=2, text="Height [mm]", line=2.3, cex=0.85)
      mtext(side=3, text=panel_labels[label_cnt], cex=0.85, adj=-0.075); label_cnt <- label_cnt + 1
      mtext(side=3, text=site_dat[dd, "formal"], cex=0.85)
    }
    plot(-1, -1, xlim=c(0,2), ylim=c(0,2), col="white", xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
    legend(0, 2, c("best model", covariate_labels[1:4]), lwd=c(2,1.3,1.3,1.3,1.3), col=c("black", unlist(covariate_colors)), cex=1.3, bty='n')
    dev.off()
  }
}


##======================================
## Same, but for 10-year return period

## calculated return periods and levels are the mean of a 21-year moving window
## average, centered on the year in question
metric <- "NPS"
rp_old <- "10"
present <- 2020
year <- 2050
either_side <- 10 # window width is 2*either_side + 1
yi <- match(year, years_proj)
rp_new <- rl_old <- rl_new <- rl_smoothed <- vector("list",length(names_evm))
names(rp_new) <- names(rl_old) <- names(rl_new) <- names(rl_smoothed) <- names_evm
for (gg in names_evm) {
  rp_new[[gg]] <- rl_old[[gg]] <- rl_new[[gg]] <- rep(NA, length(site_names)) # needs to be the return periods in the same order as lats and lons
  rl_smoothed[[gg]] <- vector('list', length(site_names)); names(rl_smoothed[[gg]]) <- site_names
  for (dd in 1:length(site_names)) {
    cc <- best_models_map[[gg]][[metric]]$covar[dd]
    mm <- best_models_map[[gg]][[metric]]$model[dd]
    rl_smoothed[[gg]][[dd]] <- vector("list", length(names_covariates)); names(rl_smoothed[[gg]][[dd]]) <- names_covariates
    if (cc==5) {
      # stationary model, no moving average needed here
      rl_old[[gg]][dd] <- rl[[gg]][[dd]][[1]][rp_old,mm,match(present,years_proj)] # doesn't matter which covariate you take since it's stationary
      rp_new[[gg]][dd] <- as.numeric(rp_old)
      rl_new[[gg]][dd] <- rl_old[[gg]][dd]
    } else {
      # nonstationary model, use moving average
      rl_old[[gg]][dd] <- rl[[gg]][[dd]][[cc]][rp_old,mm,match(present,years_proj)]
      rp_new[[gg]][dd] <- mean(rp[[gg]][[dd]][[cc]][rp_old,mm,(yi-either_side):(yi+either_side)])
      rl_new[[gg]][dd] <- mean(rl[[gg]][[dd]][[cc]][rp_old,mm,(yi-either_side):(yi+either_side)])
    }
    for (c in names_covariates) {
      rl_smoothed[[gg]][[dd]][[c]] <- mat.or.vec(nr=(year-present+1), nc=nmodel)
      for (yy in present:year) {
        yi2 <- match(yy, years_proj)
        for (m in 1:nmodel) {
          rl_smoothed[[gg]][[dd]][[c]][yy-present+1, m] <- mean(rl[[gg]][[dd]][[c]][rp_old,m,(yi2-either_side):(yi2+either_side)])
        }
      }
    }
  }
}

year <- 2050
yi <- match(year, years_proj)
yi_present <- match(present, years_proj)

for (metric in names(metrics)) {
  for (gg in names_evm) {
    label_cnt <- 1
    pdf(paste("../figures/projections_long_rl",rp_old,"_",gg,"_",metric,".pdf", sep=""), width=7.5, height=9, pointsize=11, colormodel='cmyk')
    par(mfrow=c(5,4), mai=c(.43,.45,.22,.15))
    for (dd in idx_long) {
      ymin <- min(rl_smoothed[[gg]][[dd]][[1]], na.rm=TRUE); for (c in names_covariates) {if (min(rl_smoothed[[gg]][[dd]][[c]], na.rm=TRUE) < ymin) ymin <- min(rl_smoothed[[gg]][[dd]][[c]], na.rm=TRUE)}
      ymax <- max(rl_smoothed[[gg]][[dd]][[1]], na.rm=TRUE); for (c in names_covariates) {if (max(rl_smoothed[[gg]][[dd]][[c]], na.rm=TRUE) > ymax) ymax <- max(rl_smoothed[[gg]][[dd]][[c]], na.rm=TRUE)}
      cc <- best_models_map[[gg]][[metric]]$covar[dd]
      mm <- best_models_map[[gg]][[metric]]$model[dd]
      if (cc==5) {cc <- 1; mm <- 1}
      # actual plots
      plot(years_proj[yi_present:yi], rl_smoothed[[gg]][[dd]][[cc]][,mm], lwd=2, type='l', col="black", xaxs='i', xlim=c(2020,year), ylim=c(ymin,ymax), main="", xaxt='n', xlab='', ylab='')
      grid()
      for (m in 1:nmodel) {for (c in names_covariates) {lines(years_proj[yi_present:yi], rl_smoothed[[gg]][[dd]][[c]][,m], lwd=0.5, col=covariate_colors[[match(c,names_covariates)]])}}
      lines(years_proj[yi_present:yi], rl_smoothed[[gg]][[dd]][[cc]][,mm], lwd=2.5, col="black")
      axis(1, at=seq(from=2020, to=year, by=5), labels=rep("",length(seq(from=2020, to=year, by=5))))
      axis(1, at=seq(from=2020, to=year, by=10), labels=seq(from=2020, to=year, by=10))
      #axis(2, at=seq(from=ymin,to=ymax,by=500), labels=seq(from=ymin,to=ymax,by=500), las=1)
      mtext(side=1, text="Year", line=2.3, cex=0.85)
      mtext(side=2, text="Height [mm]", line=2.3, cex=0.85)
      mtext(side=3, text=panel_labels[label_cnt], cex=0.85, adj=-0.075); label_cnt <- label_cnt + 1
      mtext(side=3, text=site_dat[dd, "formal"], cex=0.85)
    }
    plot(-1, -1, xlim=c(0,2), ylim=c(0,2), col="white", xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
    legend(0, 2, c("best model", covariate_labels[1:4]), lwd=c(2,1.3,1.3,1.3,1.3), col=c("black", unlist(covariate_colors)), cex=1.3, bty='n')
    dev.off()
  }
}

##==============================================================================




##==============================================================================
## Check if any of the GPD rate-nonstationary
## models fits worse than a stationary model
##===========================================

print(any(metrics$NPS$gpd[,,1]-metrics$NPS$gpd[,,2] < 0))
print(apply(metrics$NPS$gev[,,3]-metrics$NPS$gev[,,1], 2, mean))
print(apply(metrics$NPS$gpd[,,5]-metrics$NPS$gpd[,,1], 2, mean))
##==============================================================================



##==============================================================================
## Plots and maps for only the 1-parameter models
## WARNING: at this point, the script below is going to overwrite all of the
## main analysis arrays and lists. For example, the `best_models_datlen` list.

source("plots_one_parameter_models.R")
##==============================================================================



##==============================================================================
## End
##==============================================================================
