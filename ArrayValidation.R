source("./config.R")

library(tidyverse)
library(glmnet)
library(zeallot)
library(readxl)
library(plyr)
library(gplots)

select = dplyr::select

source("./LASSOFunctions.R")

# Suppress warnings
defaultW <- getOption("warn")
options(warn=-1)

# Parse command line argument for index
i <- as.integer(commandArgs()[6])

xy <- readRDS(xy_rds)

alphas <- seq(0, 1, length.out=20)

# Format for reading data in
demark <- formatC(i, width=3, flag="0")

# Retrieve alpha grid
rds_out <- readRDS(paste(grid_path, demark, ".rds", sep="")); alpha_df <- rds_out[[3]]

# Retrieve test set
test_indices <- read.table(paste(grid_path, demark, "_indices.txt", sep=""))$x

dat_df <- LOO_validation(xy, alpha_df, test_indices)

dat_df$trial <- i

write.table(dat_df, val_csv, row.names=F, append=T, sep=",", col.names=!file.exists(val_csv))

# Return warning status
options(warn=defaultW)
