#' Select only observations which have a sufficient coverage of measurements during the night
#'
#' @param df data frame with the observations
#'
#' @param predictor for what predictor is the check made?
#'
#' @param observation_coverage data frame with the observations coverage information
#'
#' @param coverage_factor what is the minimal coverage required?
#'
#' @return data frame with observations of sufficient quality
ensure_quality_observations <- function(df, predictor = "HEART_RATE_mean", observation_coverage, coverage_factor) {
  if (str_starts(predictor, "HEART_RATE_")) {
    valid_observations <- observation_coverage %>%
      filter(HEART_RATE_n_obs >= coverage_factor) %>%
      select(ID, rel_day)
  } else if (str_starts(predictor, "RESPIRATION_")) {
    valid_observations <- observation_coverage %>%
      filter(RESPIRATION_n_obs >= coverage_factor) %>%
      select(ID, rel_day)
  } else if (str_starts(predictor, "hrv_")) {
    valid_observations <- observation_coverage %>%
      filter(hrv_n_windows_notna >= coverage_factor) %>%
      select(ID, rel_day)
  } else {
    stop("Predictor could not be classified")
  }

  return(df %>% inner_join(valid_observations, by = c("ID", "rel_day")))
}

#' Use LASSO to select PCs with predictive value
#'
#' @param pca_df data frame with the PCs
#'
#' @param seed a seed for the LASSO regression
#'
#' @return formula-style string with selected PCs
get_PCs_LASSO <- function(pca_df, seed = 123) {
  suppressWarnings({
    set.seed(seed)

    # discharge
    cl_data_uninformative <- pca_df %>%
      mutate(event = ifelse(event == "icu", "stay", as.character(event))) %>%
      select(c(event, age, sex), starts_with("PC")) %>%
      drop_na()

    cv.lasso <- cv.glmnet(data.matrix(cl_data_uninformative %>% select(-c(event))), cl_data_uninformative$event, alpha = 1, family = "binomial")

    lambda_opt <- cv.lasso$lambda[min(which(cv.lasso$cvm <= 0 * cv.lasso$cvsd[which.min(cv.lasso$cvm)] + cv.lasso$cvm[which.min(cv.lasso$cvm)]))]

    PC_select_discharge <- as.data.frame(as.matrix(coef(cv.lasso, lambda_opt))) %>%
      rownames_to_column() %>%
      set_colnames(c("name", "value")) %>%
      mutate(value = round(value, 10)) %>%
      filter(value != 0) %>%
      arrange(desc(abs(value)))

    # ICU
    cl_data_uninformative <- pca_df %>%
      mutate(event = ifelse(event == "dismiss", "stay", as.character(event))) %>%
      select(c(event, age, sex), starts_with("PC")) %>%
      drop_na()

    cv.lasso_icu <- cv.glmnet(data.matrix(cl_data_uninformative %>% select(-c(event))), cl_data_uninformative$event, alpha = 1, family = "binomial")

    lambda_opt <- cv.lasso_icu$lambda[min(which(cv.lasso_icu$cvm <= 0 * cv.lasso_icu$cvsd[which.min(cv.lasso_icu$cvm)] + cv.lasso_icu$cvm[which.min(cv.lasso_icu$cvm)]))]
    PC_select_icu <- as.data.frame(as.matrix(coef(cv.lasso_icu, lambda_opt))) %>%
      rownames_to_column() %>%
      set_colnames(c("name", "value")) %>%
      mutate(value = round(value, 10)) %>%
      filter(value != 0) %>%
      arrange(desc(abs(value)))

    PC_select <- rbind(PC_select_discharge, PC_select_icu) %>%
      filter(grepl("PC", name)) %>%
      arrange(name)
  })
  return(PC_select %>% pull(name) %>% paste(collapse = " + "))
}

#' Apply PCA on "new" observations. This is used during cross-validation to also include the PCA in the
#' leave-one-patient out CV scheme (PCs are only defined based on the training data).
#'
#' @param training Data frame with the training data
#'
#' @param newdatalist List of data frames with new (test) data
#'
#' @param variables Vector with the names of the variables to include in the PCA
#'
#' @param standardize Should the PCs be standardized?
#'
#' @return A list of data frames with the PCs for the new (test) data
transform_pca <- function(training, newdatalist, variables, standardize = T) {

  # run PCA
  pca_result <- prcomp(training %>% select(all_of(variables)) %>% drop_na(), scale. = T)

  if (standardize) {
    pc_predict <- data.frame(predict(pca_result, training))
    pc_mean <- pc_predict %>% summarize(across(everything(), mean, na.rm = T))
    pc_sd <- pc_predict %>% summarize(across(everything(), sd, na.rm = T))
  }

  result_df_list <- lapply(newdatalist, function(newdata) {
    # predict PCs and add to original df
    result_df <- newdata %>%
      bind_cols(data.frame(predict(pca_result, newdata))) %>%
      relocate(starts_with("PC"), .after = sex)

    if (standardize) {
      result_df <- bind_rows(pc_mean, pc_sd, result_df) %>%
        mutate(across(starts_with("PC"), ~ (. - .[1]) / .[2])) %>%
        slice(3:n())
    }

    return(result_df)
  })

  return(result_df_list)
}
