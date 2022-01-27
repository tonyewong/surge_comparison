##==============================================================================
## install_packages.R
##
## Install relevant R packages and data that will be used in this program.
##
## Tony Wong (aewsma@rit.edu)
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

# packages
install.packages("abind")
install.packages("coda")
install.packages("extRemes")
install.packages("ncdf4")
install.packages("DEoptim")
install.packages("Hmisc")
install.packages("date")
install.packages("Bolstad")
install.packages("xlsx")
install.packages("plot3D")
install.packages("maps")

# data - temperature projections
local <- "../input_data/global.tas.aann.CNRM-CM5.historical+rcp85.r1i1p1.18500101-21001231.nc"
remote <- "https://drive.google.com/u/1/uc?id=18b91t463kHr-bpuWR0CTWcHVuQt_DhMp&export=download"
download.file(remote, local)

# WARNING! These two might take a while depending on your internet speed.

# data - sea-level projections
local <- "../input_data/BRICK_GMSL_WongKeller2017.nc"
remote <- "https://drive.google.com/u/1/uc?id=1iU32zKuz_UV0QEVOf2ZbX8sdXMNOc3LH&export=download"
download.file(remote, local)

# data - sea-level pressure projections
local <- "../input_data/DMIEH5_SRA1B_4_MM_psl.1-1200.nc"
remote <- "https://drive.google.com/u/1/uc?id=1M74VXYA2j-BLMXOU9npPs_bBRTe4p6sE&export=download"
download.file(remote, local)

##==============================================================================
## End
##==============================================================================
