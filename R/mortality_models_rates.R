library(mgcv)
library(data.table)

# 1. The Fitting Function
fit_two_step_model <- function(data, dof = 10) {
  # Step 1: Fit the baseline parametric model
  # This forces the model to explain as much as possible using the standard curve
  baseline_fit <- lm(log_mx ~ 1 + log_mx_std, data = data)
  
  # Step 2: Calculate the residuals
  # These are the localized bumps/dips that the standard curve missed
  resids <- vector(mode='numeric',length = nrow(data))
  resids[] <- 0
  resids[!is.na(data$log_mx)] <- residuals(baseline_fit)
  data$resids <- resids
  
  # Step 3: Fit the GAM on the residuals to smooth the deviations across age
  smooth_fit <- gam(resids ~ s(x, k = dof), data = data)
  
  # Return both models in a list so we can use them together later
  return(list(
    baseline = baseline_fit,
    smooth = smooth_fit
  ))
}

# 2. The Prediction Function
predict_two_step_model <- function(model_list, newdata) {
  # Predict the main trend from the linear model
  base_pred <- predict(model_list$baseline, newdata = newdata)
  
  # Predict the age-specific deviations from the GAM
  smooth_pred <- predict(model_list$smooth, newdata = newdata)
  
  # Combine them for the final predicted log(mx)
  final_pred <- base_pred + smooth_pred
  
  return(final_pred)
}

# 3. Smooth predictions with Kannisto model
fit_kannisto <- function(model_list, newdata,limit_ages=c(80,95)) {
  
  # Get predicted rates
  predicted_rates <- exp(predict_two_step_model(model_list,newdata = copy(newdata)))
  
  # Create model data
  model_data <- data.table(
    x=newdata$x,
    mx=predicted_rates
  )
  
  # Add logit of rates
  model_data[,logit_mx:=log(mx/(1-mx))]
  
  # Fit Kannisto model
  kannisto_fit <- lm(logit_mx ~ 1 + x,data=model_data[between(x,limit_ages[1],limit_ages[2])])
  
  return(kannisto_fit)
}

#' Fit a Flexible Relational Mortality Model (Rates)
#'
#' @description 
#' `singularize` converts abridged or sparse empirical mortality rates into a 
#' smooth, single-year mortality schedule. It achieves this by anchoring the 
#' empirical rates to a standard demographic schedule and using a GAM to smooth 
#' the localized residuals.
#'
#' @details 
#' This function models the log-transformed mortality rates using a sequential, 
#' two-step architecture. To prevent the flexible spline from absorbing the 
#' structural variance that should be explained by the standard curve (concurvity), 
#' the function first fits an OLS model to extract the global trend. A penalized 
#' spline is then fit exclusively to the residuals of the baseline model. 
#' 
#' Note: This approach assumes constant variance across log-rates and will fail 
#' if `mx` contains exact zeros. For data with structural zeros or severe 
#' exposure imbalances, a count-based model is recommended.
#'
#' @param mx_data A `data.table` or `data.frame` containing the empirical data. 
#'   Must include columns `x` (age) and `mx` (observed mortality rate).
#' @param mx_std_data A `data.table` containing the standard mortality schedule. 
#'   Must include columns `x` (age) and `mx_std` (standard mortality rate).
#' @param dof Integer. The dimension of the basis used to represent the smooth 
#'   term in the GAM (the `k` argument in `s()`). Controls maximum flexibility 
#'   of the residuals. Default is 10.
#' @param skip_open Logical. If TRUE, the highest age in `mx_data` is treated 
#'   as NA and the model interpolates/extrapolates over it. Default is FALSE.
#' @param skip_ages A numeric vector of length 2, or a list of such vectors, 
#'   defining lower and upper bounds of ages to exclude (treat as NA) during 
#'   fitting. Useful for bridging over known data anomalies or age-heaping. 
#'   Default is NULL.
#' @param ages A numeric vector defining the granular age grid for continuous 
#'   predictions. Default is `c(0.14, 1.5:109.5, 111.5)`.
#' @param smoooth_old_ages Logical. If TRUE, fits a Kannisto (logistic) model 
#'   to the older ages to ensure asymptotic mortality closure. Default is FALSE.
#' @param kannisto_limit_ages Numeric vector of length 2 defining the age bounds 
#'   used to fit the Kannisto extrapolation. Default is `c(80, 95)`.
#' @param kannisto_extrapolate Numeric. The exact age at which the Kannisto 
#'   extrapolation seamlessly takes over from the GAM predictions. Default is 85.
#'
#' @return A `data.table` containing the original granular `ages`, the merged 
#'   standard rates, the empirical data (if matched), and the newly predicted 
#'   smoothed continuous mortality rates in a column named `mx_s`.
#' 
#' @export

# Putting all together
singularize <- function(
    mx_data,
    mx_std_data,
    dof=10,
    skip_open=F,
    skip_ages=NULL,
    ages=c(0.14,1.5:109.5,111.5),
    smoooth_old_ages=F,
    kannisto_limit_ages=c(80,95),
    kannisto_extrapolate=85
    ) {
  
  data <- merge(mx_data,mx_std_data,by=c('x'),all.x=T)
  
  data[,log_mx:=log(mx)]
  data[,log_mx_std:=log(mx_std)]
  
  if (skip_open==T) {
    data[x==max(x),log_mx:=NA]
  }
  
  # NEW: Handle either a single vector or a list of vectors for skip_ages
  if (!is.null(skip_ages)) {
    # If a single vector is passed, wrap it in a list so we can iterate uniformly
    if (!is.list(skip_ages)) {
      skip_ages <- list(skip_ages)
    }
    
    # Iterate over the list of bounds
    for (age_range in skip_ages) {
      if (is.numeric(age_range) && length(age_range) == 2) {
        data[x %between% age_range, log_mx := NA]
      } else {
        warning("Elements in skip_ages must be numeric vectors of length 2 (e.g., c(min, max)). Ignoring malformed input.")
      }
    }
  }
  
  # Create predicted rates with granular age
  new_data <- expand.grid(
    x=ages
  )
  
  model_list <- fit_two_step_model(data,dof=dof)
  
  setDT(new_data)
  
  new_data <- merge(new_data,mx_std_data,by='x')
  new_data <- merge(new_data,mx_data,by='x',all.x=T)
  
  new_data[,log_mx_std:=log(mx_std)]
  
  predicted_rates <- exp(predict_two_step_model(model_list,newdata = copy(new_data)))
  
  if (smoooth_old_ages==T) {
    kannisto_fit <- fit_kannisto(model_list,new_data,kannisto_limit_ages)
    predicted_rates[new_data$x>kannisto_extrapolate] <- inv_logit(predict(kannisto_fit,newdata = new_data[x>kannisto_extrapolate]))
  }

  new_data[,mx_s:=predicted_rates]
  new_data
  
}

inv_logit <- function(x) {
  exp(x)/(1+exp(x))
}
