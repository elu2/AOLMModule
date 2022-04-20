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
i <- as.integer(commandArgs()[6]); demark <- formatC(i, width=3, flag="0")

# Read in data
xy <- readRDS(xy_rds)

c(res_xy, test_xy, test_indices) %<-% train_test(xy, 0.9)

# Get frequency of predictors
freq_df <- frequency_test(res_xy, n_iter=n_iter)
top_ten_vars <- freq_df %>% filter(freq >= n_iter * .1) %>% pull(ID)

# Construct new test and train set for alpha optimization
ss_xy <- res_xy %>% select(top_ten_vars, y)
test_x <- as.matrix(test_xy %>% select(top_ten_vars))
test_y <- test_xy %>% pull(y)

# Get best model using frequent predictors and make a prediction
mat_out <- mod_alpha_testing(ss_xy)

errs <- mat_out$Errors
best_model <- mat_out$Fits[which(errs == min(errs))][[1]]

pred <- predict(best_model, newx=test_x, type="response", s=best_model$lambda.min)

# Write predictions into dataframe and save in csv
dat_df <- data.frame(pred=pred[,1], actual=test_y, i=i)

write.table(dat_df, freq_val_csv, row.names=F, append=T, sep=",", col.names=!file.exists(freq_val_csv))
