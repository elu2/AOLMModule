source("./config.R")

library(tidyverse)
library(glmnet)
library(zeallot)
library(readxl)
library(plyr)
library(gplots)

select = dplyr::select

source("./LASSOFunctions.R")

defaultW <- getOption("warn")
options(warn=-1)

xy <- readRDS(xy_rds)

i <- as.integer(commandArgs()[6])

# Reformat string for writing out to RDS
demark <- formatC(i, width=3, flag="0")

# Split data into training and testing sets
c(train, test, test_i) %<-% train_test(xy, .9)
test_y <- test %>% pull(y)

# Test a grid of alpha values
out <- LOOCV_alpha_testing(train, 20, family=family, type.measure=type.measure)
alpha_df <- out[[3]]

# Save output as rds
saveRDS(out, paste(grid_path, demark, ".rds", sep=""))
# Write out test indices
write.table(test_i, paste(grid_path, demark, "_indices.txt", sep=""))

options(warn=1)
