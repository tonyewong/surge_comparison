##==============================================================================
## plots_one_parameter_models.R
##
## To be run at the very end of `analysis_driver.R`. This routine will perform
## the same comparisons and generate the same plots as the main analyses, but
## will only consider 1-parameter models with either sigma, mu (for GEV) or
## lambda (for GPD) nonstationary. xi-nonstationary is not considered in view
## of the arguments against it in (for example) Ceres et al. (2017; doi:
## 10.1007/s10584-017-2075-0).
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
model_choices <- c(1,2,3)
covariate_labels <- c("Time","Temp","Sea level","NAO","Stat") # needs to be in the same order as the second dimension of the `metrics` array: [1] "time" "temp" "sealevel" "nao"
covariate_breaks <- sort(c(seq(from=0.75, to=4.75, by=1), seq(from=1.25, to=5.25, by=1)))
covariate_breaks_short <- sort(c(seq(from=0.75, to=4.75, by=1), seq(from=1.001, to=5.001, by=1)))
covariate_breaks_long <- sort(c(seq(from=0.999, to=4.999, by=1), seq(from=1.25, to=5.25, by=1)))

h_offset <- 0.025
pdf('../figures/covariate_choice_1param.pdf', width=7.5, height=3.5, pointsize=11, colormodel='cmyk')
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
model_choices <- c(1,2,3)
model_labels <- vector('list', 2); names(model_labels) <- names_evm
model_labels$gev <- c("none",expression(mu),expression(sigma),expression(xi),expression(mu * ', ' * sigma),expression(mu * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * mu * ', ' * sigma * ', ' * xi))
model_labels$gpd <- c("none",expression(lambda),expression(sigma),expression(xi),expression(lambda * ', ' * sigma),expression(lambda * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * lambda * ', ' * sigma * ', ' * xi))
model_breaks_short <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.001, to=8.001, by=1)))
model_breaks_long <- sort(c(seq(from=0.999, to=7.999, by=1), seq(from=1.4, to=8.4, by=1)))
model_breaks <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.4, to=8.4, by=1)))

pdf('../figures/model_choice_1param.pdf', width=6, height=8, pointsize=11, colormodel='cmyk')
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
    hist(best_models_all[[gg]][[metric]], breaks=model_breaks, freq=TRUE, lty=5, xlim=c(0.5,8.5), ylim=c(0,36), main="", xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i')
    grid()
    hist(best_models_datlen$short[[gg]][[metric]], breaks=model_breaks_short, freq=TRUE, col="gray80", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    hist(best_models_datlen$long[[gg]][[metric]], breaks=model_breaks_long, freq=TRUE, col="gray30", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    if (label_cnt==1) {legend(0.6,37,c("< 55 years","> 55 years","combined"),pch=c(15,15,NA),lty=c(NA,NA,5),col=c("gray80","gray30","black"),pt.cex=2,bty='n')}
    axis(1, at=1:8, labels=model_labels[[gg]])
    axis(2, at=seq(from=0,to=36,by=5), labels=seq(from=0,to=36,by=5), las=1)
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
  model_choices <- c(1,2,3)
  model_labels <- vector('list', 2); names(model_labels) <- names_evm
  model_labels$gev <- c("none",expression(mu),expression(sigma),expression(xi),expression(mu * ', ' * sigma),expression(mu * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * mu * ', ' * sigma * ', ' * xi))
  model_labels$gpd <- c("none",expression(lambda),expression(sigma),expression(xi),expression(lambda * ', ' * sigma),expression(lambda * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * lambda * ', ' * sigma * ', ' * xi))
  model_breaks_short <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.001, to=8.001, by=1)))
  model_breaks_long <- sort(c(seq(from=0.999, to=7.999, by=1), seq(from=1.4, to=8.4, by=1)))
  model_breaks <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.4, to=8.4, by=1)))

  pdf(paste('../figures/model_choice_',cc,'_1param.pdf',sep=''), width=6, height=8, pointsize=11, colormodel='cmyk')
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
  model_choices <- c(1,2,3)
  model_labels <- vector('list', 2); names(model_labels) <- names_evm
  model_labels$gev <- c("none",expression(mu),expression(sigma),expression(xi),expression(mu * ', ' * sigma),expression(mu * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * mu * ', ' * sigma * ', ' * xi))
  model_labels$gpd <- c("none",expression(lambda),expression(sigma),expression(xi),expression(lambda * ', ' * sigma),expression(lambda * ', ' * xi),expression(sigma * ', ' * xi),expression(' ' * lambda * ', ' * sigma * ', ' * xi))
  model_breaks_short <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.001, to=8.001, by=1)))
  model_breaks_long <- sort(c(seq(from=0.999, to=7.999, by=1), seq(from=1.4, to=8.4, by=1)))
  model_breaks <- sort(c(seq(from=0.6, to=7.6, by=1), seq(from=1.4, to=8.4, by=1)))

  pdf(paste('../figures/model_choice_',metric,'_1param.pdf', sep=''), width=6.5, height=9, pointsize=11, colormodel='cmyk')
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
    # hist(best_models_datlen$all[[gg]], breaks=model_breaks, freq=TRUE, lty=5, xlim=c(0.5,8.5), ylim=c(0,ymaxs[metric]), main="", xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i')
    # grid()
    # hist(best_models_datlen$short[[gg]], breaks=model_breaks_short, freq=TRUE, col="gray80", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    # hist(best_models_datlen$long[[gg]], breaks=model_breaks_long, freq=TRUE, col="gray30", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    # if (label_cnt==1) {legend(5,ymaxs[metric],c("< 55 years","> 55 years","combined"),pch=c(15,15,NA),lty=c(NA,NA,5),col=c("gray80","gray30","black"),pt.cex=2,bty='n')}
    # axis(1, at=1:8, labels=model_labels[[gg]])
    # axis(2, at=seq(from=0,to=ymaxs[metric],by=5), labels=seq(from=0,to=ymaxs[metric],by=5), las=1)
    # mtext(side=1, text="Nonstationary parameters", line=2.3, cex=0.85)
    # mtext(side=2, text="# stations", line=2.3, cex=0.85)
    # if (gg=="gev") {mtext(side=2, text="All covariates", line=4, cex=0.85)}
    # if(gg=="gev") {mtext(side=3, text="GEV", line=.6, cex=0.85)} else {mtext(side=3, text="GPD", line=.6, cex=0.85)}
    # mtext(side=3, text=panel_labels[label_cnt], cex=0.85, adj=0); label_cnt <- label_cnt + 1
    hist(best_models_datlen$all[[gg]], breaks=model_breaks, freq=TRUE, lty=5, xlim=c(0.5,length(model_choices)+.5), ylim=c(0,ymaxs[metric]), main="", xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i', col="white")
    grid()
    hist(best_models_datlen$short[[gg]], breaks=model_breaks_short, freq=TRUE, col="gray80", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    hist(best_models_datlen$long[[gg]], breaks=model_breaks_long, freq=TRUE, col="gray30", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
    if (label_cnt==1) {legend(0.6,ymaxs[metric]+1,c("< 55 years","> 55 years","combined"),pch=c(15,15,NA),lty=c(NA,NA,5),col=c("gray80","gray30","black"),pt.cex=2,bty='n')}
    axis(1, at=1:length(model_choices), labels=model_labels[[gg]][model_choices])
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
      # hist(best_models_datlen$all[[gg]], breaks=model_breaks, freq=TRUE, lty=5, xlim=c(0.5,8.5), ylim=c(0,ymaxs[metric]), main="", xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i')
      # grid()
      # hist(best_models_datlen$short[[gg]], breaks=model_breaks_short, freq=TRUE, col="gray80", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
      # hist(best_models_datlen$long[[gg]], breaks=model_breaks_long, freq=TRUE, col="gray30", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
      # if (label_cnt==1) {legend(5,ymaxs[metric],c("< 55 years","> 55 years","combined"),pch=c(15,15,NA),lty=c(NA,NA,5),col=c("gray80","gray30","black"),pt.cex=2,bty='n')}
      # axis(1, at=1:8, labels=model_labels[[gg]])
      # axis(2, at=seq(from=0,to=ymaxs[metric],by=5), labels=seq(from=0,to=ymaxs[metric],by=5), las=1)
      # mtext(side=1, text="Nonstationary parameters", line=2.3, cex=0.85)
      # mtext(side=2, text="# stations", line=2.3, cex=0.85)
      # if (gg=="gev") {mtext(side=2, text=covariate_labels[match(cc,names_covariates)], line=4, cex=0.85)}
      # mtext(side=3, text=panel_labels[label_cnt], cex=0.85, adj=0); label_cnt <- label_cnt + 1
      hist(best_models_datlen$all[[gg]], breaks=model_breaks, freq=TRUE, lty=5, xlim=c(0.5,length(model_choices)+.5), ylim=c(0,ymaxs[metric]), main="", xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i', col="white")
      grid()
      hist(best_models_datlen$short[[gg]], breaks=model_breaks_short, freq=TRUE, col="gray80", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
      hist(best_models_datlen$long[[gg]], breaks=model_breaks_long, freq=TRUE, col="gray30", xlim=c(0.5,8.5), ylim=c(0,30), main="", xaxt='n', yaxt='n', xlab='', ylab='', add=TRUE)
      if (label_cnt==1) {legend(0.6,ymaxs[metric],c("< 55 years","> 55 years","combined"),pch=c(15,15,NA),lty=c(NA,NA,5),col=c("gray80","gray30","black"),pt.cex=2,bty='n')}
      axis(1, at=1:length(model_choices), labels=model_labels[[gg]][model_choices])
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

best_models_map <- vector("list", 2); names(best_models_map) <- names_evm
for (gg in names_evm) {best_models_map[[gg]] <- vector("list", length(metrics)); names(best_models_map[[gg]]) <- names(metrics)}
model_choices <- c(1,2,3)
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

  png(paste("../figures/covariate_choice_map_",metric,"_1param.png", sep=""), width=10, height=6, units="in", res=300, pointsize=11)
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

##==============================================================================



##==============================================================================
## End
##==============================================================================
