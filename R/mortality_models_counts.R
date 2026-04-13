library(mgcv)
library(data.table)

# 1. The Fitting Function
fit_two_step_model <- function(data, dof = 10) {
  # We require Ex (exposure) and Dx (deaths) in the data
  if(any(data$Ex <= 0, na.rm = TRUE)) stop("Exposure must be strictly positive.")
  
  # Step 1: Fit the baseline parametric model using a Negative Binomial likelihood
  # log(Ex) acts as the offset to convert the count process to a rate
  baseline_fit <- gam(Dx ~ 1 + log_mx_std + offset(log(Ex)), 
                      family = nb(), data = data)
  
  # Step 2: Extract the baseline prediction on the log scale
  # predict(..., type="link") includes the offset. We subtract log(Ex) so eta_base
  # represents pure log(mx) from the baseline model.
  eta_base <- predict(baseline_fit, type = "link",newdata = data) - log(data$Ex)
  data$eta_base <- eta_base
  
  # Step 3: Fit the GAM for the deviations
  # We remove the intercept (-1) so the GAM only models the structural deviations.
  # The new offset is the exposure plus the baseline trend.
  smooth_fit <- gam(Dx ~ -1 + s(x, k = dof) + offset(log(Ex) + eta_base), 
                    family = nb(), data = data)
  
  return(list(
    baseline = baseline_fit,
    smooth = smooth_fit
  ))
}

# 2. The Prediction Function
predict_two_step_model <- function(model_list, newdata) {
  # Create a mock dataset to safely predict underlying rates, not counts.
  newdata_mock <- copy(newdata)
  
  # 1. Mock the exposure to 1 for ALL rows (even extrapolated ones)
  # This prevents any log(NA) issues and zeroes out the log(Ex) offset.
  newdata_mock$Ex <- 1 
  
  # 2. Predict baseline log-rate
  # Since log(Ex) = 0, this outputs pure log(mx) from the parametric trend
  base_log_mx <- predict(model_list$baseline, newdata = newdata_mock, type = "link")
  
  # 3. Satisfy the GAM's formula requirements
  # mgcv remembers the offset(log(Ex) + eta_base) formula and expects 
  # eta_base to exist in newdata, even though type="terms" won't mathematically 
  # use it. We supply base_log_mx to prevent "variable missing" errors.
  newdata_mock$eta_base <- base_log_mx
  
  # 4. Predict the age-specific deviations
  # Using type="terms" isolates the s(x) component matrix.
  # We use newdata_mock to guarantee no NAs cause dropped rows.
  smooth_terms <- predict(model_list$smooth, newdata = newdata_mock, type = "terms")
  smooth_log_mx <- rowSums(smooth_terms) 
  
  # 5. Combine them for the final predicted log(mx)
  # Ensure base_log_mx is treated as a flat numeric vector for safe addition
  final_pred <- as.numeric(base_log_mx) + as.numeric(smooth_log_mx)
  
  return(final_pred)
}

# 3. Smooth predictions with Kannisto model
# (Adapted to use the output of the new prediction function)
fit_kannisto <- function(model_list, newdata, limit_ages = c(80, 95)) {
  
  # Get predicted rates
  predicted_rates <- exp(predict_two_step_model(model_list, newdata = copy(newdata)))
  
  # Create model data
  model_data <- data.table(
    x = newdata$x,
    mx = predicted_rates
  )
  
  # Add logit of rates
  model_data[, logit_mx := log(mx / (1 - mx))]
  
  # Fit Kannisto model (Consider upgrading this to a Binomial GLM in the future)
  kannisto_fit <- lm(logit_mx ~ 1 + x, data = model_data[between(x, limit_ages[1], limit_ages[2])])
  
  return(kannisto_fit)
}

#' Fit a Flexible Relational Mortality Model (Counts)
#'
#' @description 
#' `singularize` converts abridged or sparse empirical mortality counts into a 
#' smooth, single-year mortality schedule. It achieves this by anchoring the 
#' empirical data to a standard demographic schedule and using a GAM to smooth 
#' the localized deviations.
#'
#' @details 
#' This function models the data-generating process directly using a Negative 
#' Binomial likelihood. Unlike models that rely on log-transformed rates, this 
#' approach naturally handles intervals with zero deaths and correctly weights 
#' age groups based on their exposure volume.
#' 
#' To prevent the flexible spline from absorbing the structural variance that 
#' should be explained by the standard curve (concurvity), the model fits 
#' sequentially. It first estimates the global trend against the standard using 
#' log(exposure) as an offset, and then passes that baseline fit as an additional 
#' offset to a penalized spline that models age-specific deviations.
#'
#' @param mx_data A `data.table` or `data.frame` containing the empirical data. 
#'   Must include columns `x` (age), `Dx` (death counts), and `Ex` (exposure).
#' @param mx_std_data A `data.table` containing the standard mortality schedule. 
#'   Must include columns `x` (age) and `mx_std` (standard mortality rate).
#' @param dof Integer. The dimension of the basis used to represent the smooth 
#'   term in the GAM (the `k` argument in `s()`). Controls maximum flexibility 
#'   of the residuals. Default is 10.
#' @param skip_open Logical. If TRUE, the highest age in `mx_data` is treated 
#'   as NA (its deaths are ignored) and the model extrapolates over it. Default is FALSE.
#' @param skip_ages A numeric vector of length 2, or a list of such vectors, 
#'   defining lower and upper bounds of ages to exclude (treat `Dx` as NA) during 
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
    mx_data, # Now expecting Dx and Ex columns
    mx_std_data,
    dof = 10,
    skip_open = FALSE,
    skip_ages = NULL,
    ages = c(0.14, 1.5:109.5, 111.5),
    smoooth_old_ages = FALSE,
    kannisto_limit_ages = c(80, 95),
    kannisto_extrapolate = 85
) {
  
  data <- merge(mx_data, mx_std_data, by = c('x'), all.x = TRUE)
  data[, log_mx_std := log(mx_std)]
  
  if (skip_open == TRUE) {
    data[x == max(x), Dx := NA]
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
        data[x %between% age_range, Dx := NA]
      } else {
        warning("Elements in skip_ages must be numeric vectors of length 2 (e.g., c(min, max)). Ignoring malformed input.")
      }
    }
  }
  
  model_list <- fit_two_step_model(data, dof = dof)
  
  new_data <- data.table(x = ages)
  new_data <- merge(new_data, mx_std_data, by = 'x')
  new_data <- merge(new_data, mx_data, by = 'x', all.x = TRUE)
  new_data[, log_mx_std := log(mx_std)]
  
  predicted_rates <- exp(predict_two_step_model(model_list, newdata = copy(new_data)))
  
  if (smoooth_old_ages == TRUE) {
    kannisto_fit <- fit_kannisto(model_list, new_data, kannisto_limit_ages)
    extrap_idx <- new_data$x > kannisto_extrapolate
    predicted_rates[extrap_idx] <- inv_logit(predict(kannisto_fit, newdata = new_data[extrap_idx]))
  }
  
  new_data[, mx_s := predicted_rates]
  return(new_data)
}

inv_logit <- function(x) {
  exp(x) / (1 + exp(x))
}