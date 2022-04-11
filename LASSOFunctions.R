library(glmnet)
library(tidyverse)


rds_convert <- function(df, rv){
  df <- rename_rv(df, rv)
}


run <- function(x, y, alpha=1, foldid=NULL){
  fit <- glmnet(x, y, family="binomial", alpha=alpha)
  cv_fit <- cv.glmnet(as.matrix(x), as.matrix(y), family="binomial", type.measure="class", alpha=alpha, foldid=foldid, nfolds=10)
  return(list(fit, cv_fit))
}


rename_rv <- function(df, rv){
  names(df)[names(df) == rv] <- "y"
  return(df)
}


train_test <- function(xy, train_prop){
  train_i <- sample(1:nrow(xy), floor(nrow(xy) * train_prop))
  test_i <- (1:nrow(xy))[!(1:nrow(xy) %in% train_i)]
  
  train_ss <- xy[train_i,][sample(1:length(train_i), length(train_i), replace=F),]
  test_ss <- xy[test_i,][sample(1:length(test_i), length(test_i), replace=F),]
  
  out <- setNames(list(train_ss, test_ss, test_i), c("train", "test", "testIndices"))
  return(out)
}


# From an x and y, extract the cross-validated error
cv_error <- function(run_x, run_y, weights){
  # Ridge-to-LASSO parameter vector
  alpha_values <- seq(0, 1, length.out=100)
  fold_ids <- sample(1:10, nrow(run_x), replace=T)
  # Vector to store errors
  min_error <- vector()
  
  i <- 1
  for (alpha in alpha_values){
    # Verbose runtime confirmation
    print(paste("Running alpha=", alpha))
    
    run_fits <- run(run_x, run_y, alpha=alpha, foldid=fold_ids)
    
    # Cross validated LASSO at specified alpha value
    cv.fit <- run_fits[[2]]
    
    # Index of a lambda value (lambda min in this case)
    index <- which(cv.fit$lambda == cv.fit$lambda.min)
    
    # ...to index the error value
    mse.min <- cv.fit$cvm[index]
    min_error[i] <- mse.min
    
    i <- i + 1
  }
  print(plot(alpha_values, min_error))
  return(min_error)
}


# fit: cross-validated fit
# Can/will need to manually change from lambda.min in the lm_i object
AIC_BIC <- function(fit){
  lm_i <- fit$index[1]
  lm_dr <- fit$glmnet.fit$dev.ratio[lm_i]
  nd <- fit$glmnet.fit$nulldev
  
  ll <- nd - nd*(1-lm_dr)
  k <- fit$glmnet.fit$df[lm_i]
  n <- fit$glmnet.fit$nobs
  
  AIC <- -ll + 2*k*(k+1)/n-k-1
  BIC <- -log(n)*k - ll
  
  AB_list <- list(AIC, BIC)
  names(AB_list) <- c("AIC", "BIC")
  
  return(AB_list)
}


# neit - dataframe with neither bistable or rebi
# secondary - dataframe with either bistable or rebi
# n - number of total rows
# prop - proportion of neit to secondary data rows

enrich_combine <- function(neit, secondary, n_rows, prop){
  n1 <- ceiling(n_rows * prop)
  n2 <- ceiling(n_rows * (1 - prop))
  
  rNeitIndices <- sample(nrow(neit), n1)
  rBiSubsetIndices <- sample(nrow(secondary), n2)
  
  df_neit <- neit[rNeitIndices, ]
  df_bi <- secondary[rBiSubsetIndices, ]
  
  combined <- rbind(df_neit, df_bi)
  
  return(combined)
}


# Function that takes data and partitions out x, y, and the data unsampled from

train_set <- function(data, n) {
  last_col <- ncol(data)
  train_indices <- sample(nrow(data), n)
  
  tr_x <- as.matrix(data[train_indices, 1:(last_col-1)])
  tr_y <- as.matrix(rev(data)[train_indices, 1])
  
  remainder <- data[-train_indices,]
  
  return(list(tr_x, tr_y, remainder))
}


# fit: fit object
# data: data to sample (randomly) from for newx
# s: specific lambda value

prediction_fun <- function(fit, data, s){
  test_indices <- sample(nrow(data), 10000)
  last_col <- ncol(data)
  te_x <- as.matrix(data[test_indices, 1:(last_col)])
  return(predict(fit, newx=te_x, s=s, type="class"))
}


# data: dataframe of data to sample from

fitter <- function(data, n){
  train_data <- train_set(data, n)
  tr_x <- train_data[[1]]
  tr_y <- train_data[[2]]
  
  fit <- glmnet(tr_x, tr_y, family="binomial")
  cv_fit <- cv.glmnet(tr_x, tr_y, family="binomial", type.measure="class")
  
  return(list(fit, cv_fit))
}


# Returns weights for samples corresponding to class imbalances. Same order as xy
y_weights <- function(xy){
  cts <- count(xy$y)
  tot <- sum(cts$freq)
  cts$Weight <- 1 - cts$freq / tot; cts$freq <- NULL; names(cts) <- c("y", "Weight")
  
  wdf <- join(data.frame(y=xy$y), cts, by="y")
  weights <- wdf$Weight
  
  return(weights)
}


nzcs <- function(coefs){
  nzc <- coefs@Dimnames[[1]][coefs@i + 1]
  nzc <- nzc[2:length(nzc)]
  
  return(nzc)
}

# Extracts the non-zero coefficients from a cross-validated model and a lambda value
extract_nzc <- function(fit, lm, family){
  coefs <- coef(fit, s=lm)
  if (family == "multinomial"){nzc <- unique(unlist(lapply(coefs, nzcs)))}
  else(nzc <- nzcs(coefs))

  return(nzc)
}


mce_fnr_fpr <- function(actual, pred){
  # Record error rates
  incorrect <- actual[!actual == pred]
  correct <- actual[actual == pred]
  val_mce <- length(incorrect) / (length(incorrect) + length(correct))
  
  fnr <- sum(!incorrect)/(sum(!incorrect)+sum(correct))
  if (is.nan(fnr)){fnr <- 0}
  fpr <- sum(incorrect)/(sum(sum(incorrect)+sum(!correct)))
  if (is.nan(fpr)){fpr <- 0}
  
  out <- setNames(list(val_mce, fnr, fpr), c("MCE", "FNR", "FPR"))
  
  return(out)
}


# Create a confusion matrix from
# x, y: data and labels for running a cross validation. (At time of creation, these are ordered)
# zflm: the zoom-fit lambda min. Previous cross-validation run's lambda.min
# iterations: number of times to run cross-validation
# folds: i.e. the proportion of data to hold out.
confusion_f <- function(x, y, zflm, iterations, weights, alpha, family="binomial", type.measure="class", folds=5){
  all_fnr <- c()
  all_fpr <- c()
  for (i in 1:iterations){
    cutoff <- floor(nrow(x)/folds)
    
    shuffle_i <- sample(nrow(x))
    while (mean(y[shuffle_i,][1:cutoff]) == 1){
      shuffle_i <- sample(nrow(x))
    }
    
    # Holdout/Validation set
    held_szx <- x[shuffle_i,][1:cutoff,]
    held_szy <- y[shuffle_i,][1:cutoff]
    
    # Shuffled training set
    szx <- x[shuffle_i,][cutoff:nrow(x),]
    szy <- y[shuffle_i,][cutoff:nrow(y)]
    
    cweights <- weights[shuffle_i][cutoff:nrow(y)]
    
    fit <- glmnet(szx, szy, family=family, type.measure=type.measure, alpha=alpha, lambda=zflm, weights=cweights)
    
    type <- type.measure
    if (type.measure == "mse"){type <- "response"}
    pred <- as.numeric(predict(fit, newx=held_szx, type=type, s=fit$lambda.min))
    
    incorrect <- held_szy[!held_szy == pred]
    correct <- held_szy[held_szy == pred]
    
    # Calculate, record false negative and false positive rates
    fnr <- sum(!incorrect)/(sum(!incorrect)+sum(correct))
    all_fnr[i] <- fnr
    fpr <- sum(incorrect)/(sum(sum(incorrect)+sum(!correct)))
    all_fpr[i] <- fpr
  }
  
  # Remove NaNs and average the false rates
  all_fnr <- na.omit(all_fnr)
  mean_fnr <- sum(all_fnr)/length(all_fnr)
  all_fpr <- na.omit(all_fpr)
  mean_fpr <- sum(all_fpr)/length(all_fpr)

  return(c(mean_fpr, mean_fnr))
}


# layers: How many iterations to run the model for. Minimum 1.
# min_params: keep resampling until at least min_params features are left after first layer
# alpha1: Elastic Net alpha value for first layer
# alphai: Elastic Net alpha value for layers 2+
# weighted: add weights to observations. Suggested only for binomial classification
# save: write out plots and errors
layered_model <- function(xy, layers, nfolds=5, alpha1=1, alphai=1, min_params=1, family="binomial", type.measure="class", weighted=T, save=F, parallel=F, verbose=F){
  # Store false rates
  fnrs <- rep(0, layers)
  fprs <- rep(0, layers)
  oErrs <- rep(0, layers)
  lms <- rep(0, layers)
  
  x <- as.matrix(xy %>% select(-y))
  y <- as.matrix(xy %>% select(y))
  
  if (weighted){weights <- y_weights(xy)}
  else{weights <- NULL}
  
  # ----- LAYER 1 -----
  l1_fit <- cv.glmnet(x, y, nfolds=nfolds, family=family, type.measure=type.measure, alpha=alpha1, weights=weights, parallel=parallel)
  l1_lm <- l1_fit$lambda.min
  l1_nzc <- extract_nzc(l1_fit, l1_lm, family=family)

  # Do not use models that have fit too few parameters
  stall_count <- 0
  while (length(l1_nzc) < min_params){
    l1_fit <- cv.glmnet(x, y, nfolds=nfolds, family=family, type.measure=type.measure, alpha=alpha1, weights=weights, parallel=parallel)
    l1_lm <- l1_fit$lambda.min
    l1_nzc <- extract_nzc(l1_fit, l1_lm, family=family)
    stall_count <- stall_count + 1
    if (stall_count > 10){print(paste("Over 10 fits attempted with min_params=", min_params, ". Consider lowering value.", sep=""))}
  }
  
  if (save){
  png("./layer1.png")
  plot(l1_fit)
  dev.off()
  }
  
  if (verbose){plot(l1_fit); print(paste(length(l1_nzc), "parameters left after first layer."))}
  
  frs <- confusion_f(x, y, l1_lm, 100, alpha=alpha1, weights=weights, folds=nfolds, family=family, type.measure=type.measure)
  oErrs[[1]] <- min(l1_fit$cvm)
  fprs[[1]] <- frs[[1]]
  fnrs[[1]] <- frs[[2]]
  lms[[1]] <- l1_lm
  
  prev_nzc <- l1_nzc
  li_lm <- l1_lm
  
  final_fit <- NULL
  # ----- LAYERS 2+ ------
  optimized <- F

  for (i in seq(2, layers)){
    if (layers == 1){li_fit <- l1_fit; li_lm <- l1_lm; final_fit <- li_fit; break}

    li_xy <- xy[append(prev_nzc, c("y"))]
    li_x <- as.matrix(li_xy %>% select(-y))
    li_y <- as.matrix(li_xy %>% select(y))
    
    li_fit <- cv.glmnet(li_x, li_y, nfolds=nfolds, family=family, type.measure=type.measure, alpha=alphai, weights=weights)
    
    # If lambda.min is no longer changing, we know it's optimized.
    if (li_fit$lambda.min == li_lm & optimized == F){
      #print(paste("Optimized at layer", i))
      optimized <- T
    }
    li_lm <- li_fit$lambda.min
    li_nzc <- extract_nzc(li_fit, li_lm, family=family)
    prev_nzc <- li_nzc
    
    #print(paste("Layer", i+1, "MCE:", min(li_fit$cvm)))
    if (save){
    png(paste("./layer", i+1, ".png", sep=""))
    plot(li_fit)
    dev.off()
    }
    
    plot(li_fit)
    
    frs <- confusion_f(li_x, li_y, li_lm, 100, alpha=alphai, weights=weights, folds=nfolds, family=family, type.measure=type.measure)
    oErrs[[i]] <- min(li_fit$cvm)
    fprs[[i]] <- frs[1]
    fnrs[[i]] <- frs[2]
    lms[[i]] <- li_lm
    final_fit <- li_fit
  }
  
  err_df <- data.frame(layer=seq(1, length(oErrs)), MCE=oErrs, FNR=fnrs, FPR=fprs)
  if (save) {write.csv(err_df, "Errors.csv", row.names=F)}
  
  out <- setNames(list(prev_nzc, lms, err_df, oErrs, fnrs, fprs, coef(li_fit, li_lm), final_fit), c("FinalCoefficients", "FinalLambdaMin", "ErrorTable", "MCE", "FNR", "FPR", "FinalCoefficientValues", "final_fit"))
  
  return(out)
}


LOOCV <- function(xy, alpha, family="binomial", type.measure="class", weighted=T){
  x <- as.matrix(xy %>% select(-y))
  y <- as.matrix(xy %>% select(y))
  
  if (weighted){weights <- y_weights(xy)}
  else{weights <- NULL}

  fit <- cv.glmnet(x, y, nfolds=nrow(xy), family=family, type.measure=type.measure, weights=weights, alpha=alpha)
  
  return(fit)
}


layered_LOOCV <- function(xy, alpha, family="binomial"){
  alpha <- alpha
  l1.fit <- LOOCV(xy, alpha)
  l1.nzcs <- extract_nzc(l1.fit, l1.fit$lambda.min, family=family)
  l2.xy <- xy[,c(l1.nzcs, "y")]
  
  while (dim(l2.xy)[2] < 2){
    alpha <- alpha - 0.05
    l1.fit <- LOOCV(xy, alpha)
    l1.nzcs <- extract_nzc(l1.fit, l1.fit$lambda.min, family=family)
    l2.xy <- xy[,c(l1.nzcs, "y")]
  }
  
  if (dim(l2.xy)[2] == 2) {
    l2.fit <- l1.fit
    l2.nzcs <- extract_nzc(l2.fit, l2.fit$lambda.min, family=family)
  }
  
  else{
    l2.fit <- LOOCV(l2.xy, alpha)
    l2.nzcs <- extract_nzc(l2.fit, l2.fit$lambda.min, family=family)
  }
  
  return(list(l2.fit, l2.nzcs))
}


# Generate error rates for all alpha values between 0 and 1 for data. 
# n_alphas: number of alphas to test. Granularity.
LOOCV_alpha_testing <- function(xy, n_alphas, family="binomial", type.measure="class", weighted=T){
  alphas <- seq(0, 1, length.out=n_alphas)
  
  alpha1s <- c()
  alpha2s <- c()
  err_col <- c()
  coefs <- c()
  fits <- list()
  
  l1_errors <- vector(length=length(alphas))
  l2_errors <- matrix(nrow=length(alphas),ncol=length(alphas))
  k <- 1
  i <- 1
  for (alpha1 in alphas){
    writeLines(paste("Alpha 1: ", alpha1, sep=""))
    loocv1 <- LOOCV(xy, alpha1, family=family, type.measure=type.measure, weighted=weighted)
    nzcs <- extract_nzc(loocv1, loocv1$lambda.min, family=family)
    l1_errors[i] <- min(loocv1$cvm)
    subset <- xy[append(nzcs, c("y"))]
    
    j <- 1
    for (alpha2 in alphas){
      alpha1s <- c(alpha1s, alpha1)
      alpha2s <- c(alpha2s, alpha2)
      
      loocv2 <- LOOCV(subset, alpha2, family=family, type.measure=type.measure, weighted=weighted)
      nzcs <- extract_nzc(loocv2, loocv2$lambda.min, family=family)
      l2_errors[i, j] <- min(loocv2$cvm)
      
      err_col <- c(err_col, min(loocv2$cvm))
      coefs <- c(coefs, list(nzcs))
      fits[[k]] <- loocv2
      
      j <- j + 1
      k <- k + 1
    }
    i <- i + 1
  }
  
  gradient <- colorRampPalette(c("black", "yellow", "red"))(20)
  heatmap.2(l2_errors, Rowv=NA, Colv=NA, labRow=round(alphas, 2), labCol=round(alphas, 2), xlab="Alpha 2", ylab="Alpha 1", tracecol=NA, col=gradient)
  
  data_df <- data.frame(Alpha1=alpha1s, Alpha2=alpha2s, Errors=err_col)
  data_df$Coefs <- coefs
  data_df$Fits <- fits
  
  out <- list(l1_errors, l2_errors, data_df)
  
  return(out)
}


# n_iter: Number of validation sets to run
# train_prop: Proportion of data to partition into training set
# write_path: path to write out data to
# warnings: Show glmnet warnings
LOOAT_val <- function(xy, n_iter, train_prop=0.9, write_path="./", warnings=F){
  if (!warnings) {defaultW <- getOption("warn"); options(warn=-1)}
  
  all_fnr <- c()
  all_fpr <- c()
  val_mce <- c()
  
  for(i in 1:n_iter){
    print(paste("Starting run", i))
    # Reformat string for writing out to RDS
    demark <- formatC(i, width=3, flag="0")
    
    # Split data into training and testing sets
    c(train, test) %<-% train_test(xy, .9)
    test_y <- test %>% pull(y)
    
    # Test a grid of alpha values
    out <- LOOCV_alpha_testing(train, 20)
    alpha_df <- out[[3]]
    
    # Save output as rds
    saveRDS(out, paste(write_path, demark, ".rds", sep=""))
    # Write out test indices
    write.table(sample(1:nrow(xy), 8, replace=F), paste(write_path, demark, "_indices.txt", sep=""))
  }
  
  if (!warnings) {options(warn=defaultW)}
}


# alpha_df: dataframe of alpha grid. Output from alpha_testing.
# test_indices: indices of xy which were used to test
LOO_validation <- function(xy, alpha_df, test_indices){
  test <- xy[test_indices, ]
  test_y <- test %>% pull(y)
  train <- xy[-test_indices, ]
  
  
  # Extract species from minimum error rates' alphas
  errs <- alpha_df$Errors
  min_errs <- alpha_df[which(errs == min(errs)), ]
  min_prots <- unique(c(unlist(min_errs$Coefs)))
  agg_subset <- train %>% select(min_prots, y)
  test_x <- test %>% select(min_prots) %>% as.matrix()
  
  cv_mces <- c()
  # 1-Dimensional optimization of alphas
  for (j in 1:length(alphas)){
    mod_out <- layered_model(agg_subset, 1, nfolds=nrow(agg_subset), alpha1=alphas[j], min_params=2, family=family, type.measure=type.measure)
    cv_mces[j] <- mod_out$MCE[1]
  }
  
  # Select alpha with lowest cv MCE and re-fit the LOOCV model
  min_a1 <- alphas[which(cv_mces == min(cv_mces))]
  mod_out <- layered_model(agg_subset, 1, nfolds=nrow(agg_subset), alpha1=min_a1, min_params=2, family=family, type.measure=type.measure)
  s <- mod_out$FinalLambdaMin
  fit <- mod_out$final_fit
  
  pla2_pos <- as.numeric("P14555" %in% mod_out$FinalCoefficients)
  write.table(pla2_pos, "./pla2Pos.csv", row.names=F, append=T, sep=",", col.names=!file.exists("./pla2Pos.csv"))
  
  
  # Make prediction
  val_pred <- c(predict(fit, newx=test_x, type="response", s=fit$lambda.min))
  
  # Reformat and write
  dat_df <- data.frame(pred=val_pred, actual=test_y)
  
  return(dat_df)
}


# alpha_df: dataframe of alpha grid. Output from alpha_testing.
# test_indices: indices of xy which were used to test
LOO_mean_validation <- function(xy, alpha_df, test_indices){
  test <- xy[test_indices, ]
  test_y <- test %>% pull(y)
  train <- xy[-test_indices, ]
  
  # Extract species from minimum error rates' alphas
  errs <- alpha_df$Errors
  min_errs <- alpha_df[which(errs == min(errs)), ]
  
  preds <- data.frame()
  for (i in 1:nrow(min_errs)){
    skip <- F
    
    fit <- min_errs$Fits[[i]]
    
    coef_fit <- coef(fit)
    subset_coefs <- rownames(coef_fit)[2:length(coef_fit)]
    test_x <- as.matrix(test %>% select(subset_coefs))
    
    fit$lambda.min
    
    #pred <- tryCatch({predict(fit, newx=test_x, type="response", s=fit$lambda.min)}, error=function(e){skip <<- T; print("skipped 1")})
    #if(skip){next}
    
    pred <- predict(fit, newx=test_x, type="response", s=fit$lambda.min)
    
    pred <- t(data.frame(pred))
    preds <- rbind(preds, pred)
  }
  
  mean_preds <- colMeans(preds)
  
  dat_df <- data.frame(meanPred=mean_preds, actual=test_y, row.names=NULL)
  
  return(dat_df)
}


# actual: Actual response values
# pred: raw prediction probabilities
roc_validation <- function(actual, pred){
  # Calculate roc and auroc
  roc_data <- roc(actual, pred)
  AUC <- round(calculateAUC(roc_data$sensitivities, roc_data$specificities), 3)
  
  ggroc(roc_data, color="red") + 
    geom_label(aes(.1, .1, label=paste("AUC:", AUC))) + 
    labs(x="Specificity", y="Sensitivity", title="") + 
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size=.5)) + 
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed") + 
    xlim(1, 0) + 
    ylim(0, 1)
}
