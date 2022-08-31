#' Grouped k-fold cross-validation function used for leave-one-patient-out CV
#'
#' @param cvmodel The model to be cross-validated
#'
#' @param n_folds The number of folds to use (should be equal to the number of patients for LOO)
#'
#' @param pca_features The variables/features used in the PCA
#'
#' @param data_training The data that was used to train the model
#'
#' @param data_test The data on which the cross-validation should be performed
#'
#' @param sampling_mode The strategy for sampling from random effects for unobserved groups/patients.
#' See argument sample_new_levels in posterior_epred from brms.
#'
#' @return A data frame with the model predictions for the observations in data_test, obtained through k-fold cross-validation
kfold_manual <- function(cvmodel, n_folds, pca_features, data_training = NULL, data_test, sampling_mode = "uncertainty") {
  folds <- kfold_split_grouped(K = n_folds, x = cv_data_update_test$ID)

  if (is.null(data_training)) data <- cvmodel$data

  if (nrow(data_training) != nrow(cvmodel$data)) {
    print("CV training data and data from model are not identical!")
  }

  predictions <- list()

  for (fold in unique(folds)) {
    print(paste("Fitting model for fold", fold, "out of", n_folds))

    left_out <- data_test[folds == fold, ]
    subset <- data_training %>% filter(!(ID %in% unique(left_out$ID)))

    # PCA
    pca_results <- transform_pca(subset, list(subset, left_out), pca_features)

    subset_pca <- pca_results[[1]]
    left_out_pca <- pca_results[[2]]

    # fit model and predict left out observations
    updated_model <- update(cvmodel, newdata = subset_pca, recompile = F, inits = 0)
    predictions[[fold]] <- posterior_epred(updated_model, newdata = left_out_pca, allow_new_levels = T, summary = F, sample_new_levels = sampling_mode)
    # print(predictions[[fold]])
  }

  cvres <- do.call(abind, list(predictions, along = 2))

  kfold_pred <- bind_rows(
    lapply(1:dim(cvres)[1], function(i) {
      data_test %>%
        select(ID, rel_day, event) %>%
        bind_cols(
          cvres[i, , ] %>% as.data.frame() %>% rename(dismiss = 1, stay = 2, icu = 3) %>% mutate(draw = i, .before = 1)
        )
    })
  )
  return(kfold_pred)
}

# Estimation of time-dependent ROC and AUROC
library(ROCR)
source("utils/performance/NNE-estimate.R")
source("utils/performance/NNE-CrossValidation.R")
source("utils/performance/interpolate.R")
source("utils/performance/dynamicTP.R")
source("utils/performance/NNE-estimate_TPR.R")
source("utils/performance/dynamicIntegrateAUC.R")

# Mean rank function
MeanRank <- function(survival.time, survival.status, marker, start = NULL) {
  if (is.null(start)) start <- rep(0, length(survival.time))
  #####
  #
  ##### drop any missing marker values...
  #
  keep <- !is.na(survival.time) &
    !is.na(survival.status) &
    !is.na(marker) &
    !is.na(start)
  survival.time <- survival.time[keep]
  survival.status <- survival.status[keep]
  marker <- marker[keep]
  start <- start[keep]
  #
  utimes <- unique(survival.time[survival.status == 1])
  utimes <- utimes[order(utimes)]
  #
  nonparmAUC <- rep(NA, length(utimes))
  nControls <- rep(NA, length(utimes))
  TheMarker <- marker
  #
  for (j in 1:length(utimes)) {
    dead.guy <- TheMarker[(survival.time == utimes[j]) & (survival.status == 1)]
    is.started <- start < utimes[j]
    control.set <- TheMarker[is.started & (survival.time >= utimes[j]) & (survival.status == 0)] # TheMarker[ is.started & (survival.time > utimes[j]) ]
    set.size <- length(control.set)
    nControls[j] <- set.size
    ndead <- length(dead.guy)
    if (ndead == 1) {
      nonparmAUC[j] <- sum(rep(dead.guy, set.size) > control.set) / set.size
    } else {
      mean.rank <- 0
      for (k in 1:ndead) {
        this.rank <- sum(rep(dead.guy[k], set.size) > control.set) / set.size
        mean.rank <- mean.rank + this.rank / ndead
      }
      nonparmAUC[j] <- mean.rank
    }
  }
  out <- list(time = utimes, mean.rank = nonparmAUC, nControls = nControls)
  out
}

# MeanRank function for competing outcomes
MeanRank_competing <- function(survival.time, survival.status, marker, start = NULL, event_code = 1) {
  if (is.null(start)) start <- rep(0, length(survival.time))
  #####
  #
  ##### drop any missing marker values...
  #
  keep <- !is.na(survival.time) &
    !is.na(survival.status) &
    !is.na(marker) &
    !is.na(start)
  survival.time <- survival.time[keep]
  survival.status <- survival.status[keep]
  marker <- marker[keep]
  start <- start[keep]
  #
  utimes <- unique(survival.time[survival.status == event_code])
  utimes <- utimes[order(utimes)]
  #
  nonparmAUC <- rep(NA, length(utimes))
  nControls <- rep(NA, length(utimes))
  TheMarker <- marker
  #
  for (j in 1:length(utimes)) {
    dead.guy <- TheMarker[(survival.time == utimes[j]) & (survival.status == event_code)]
    is.started <- start < utimes[j]
    # as control set we only take the individuals which remain event-free through time t
    control.set <- TheMarker[is.started & (survival.time >= utimes[j]) & (survival.status == 0)] # TheMarker[ is.started & (survival.time > utimes[j]) ]
    set.size <- length(control.set)
    nControls[j] <- set.size
    ndead <- length(dead.guy)
    if (ndead == 1) {
      nonparmAUC[j] <- sum(rep(dead.guy, set.size) > control.set) / set.size
    } else {
      mean.rank <- 0
      for (k in 1:ndead) {
        this.rank <- sum(rep(dead.guy[k], set.size) > control.set) / set.size
        mean.rank <- mean.rank + this.rank / ndead
      }
      nonparmAUC[j] <- mean.rank
    }
  }
  out <- list(time = utimes, mean.rank = nonparmAUC, nControls = nControls)
  out
}

# Computation and visualization of time-dependent AUROC
get_AUC_statistics <- function(df, max_time = NULL) {
  mmmTV_dismiss <- MeanRank_competing(
    survival.time = df$tstop,
    survival.status = df$event,
    marker = df$dismiss,
    start = df$tstart,
    event_code = 1
  )

  mmmTV_icu <- MeanRank_competing(
    survival.time = df$tstop,
    survival.status = df$event,
    marker = df$icu,
    start = df$tstart,
    event_code = 2
  )

  if (!is.null(max_time)) {
    keep <- mmmTV_dismiss$time <= max_time
    mmmTV_dismiss$time <- mmmTV_dismiss$time[keep]
    mmmTV_dismiss$mean.rank <- mmmTV_dismiss$mean.rank[keep]
    mmmTV_dismiss$nControls <- mmmTV_dismiss$nControls[keep]

    keep <- mmmTV_icu$time <= max_time
    mmmTV_icu$time <- mmmTV_icu$time[keep]
    mmmTV_icu$mean.rank <- mmmTV_icu$mean.rank[keep]
    mmmTV_icu$nControls <- mmmTV_icu$nControls[keep]
  }

  # smooth
  nnn_dismiss <- nne(x = mmmTV_dismiss$time, y = mmmTV_dismiss$mean.rank, lambda = 1, nControls = mmmTV_dismiss$nControls)
  nnn_icu <- nne(x = mmmTV_icu$time, y = mmmTV_icu$mean.rank, lambda = 1, nControls = mmmTV_icu$nControls)

  # With the generalized concordance-score we here estimate the probability that the predictions for a random pair of subjects are concordant with their outcomes,
  # given that the smaller event time occurs in (0, cutoffTime). Concordance means that the subject with the higher risk score will have their event first.
  cindex_dismiss <- dynamicIntegrateAUC(
    survival.time = df$tstop,
    survival.status = df$event,
    start = df$tstart,
    marker = df$dismiss,
    event_code = 1,
    cutoffTime = 12
  )

  cindex_icu <- dynamicIntegrateAUC(
    survival.time = df$tstop,
    survival.status = df$event,
    start = df$tstart,
    marker = df$icu,
    event_code = 2,
    cutoffTime = 12
  )

  return(list(
    dismiss = list(auc = mmmTV_dismiss, auc_smooth = nnn_dismiss, cindex = cindex_dismiss),
    icu = list(auc = mmmTV_icu, auc_smooth = nnn_icu, cindex = cindex_icu)
  ))
}
