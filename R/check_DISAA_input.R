#' Validate Inputs for DISSA Analysis
#'
#' Ensures that the user has supplied either a `phyloseq_obj` or both
#' `counts` and `covariates`, that `formula` is a valid formula, and that all
#' numeric and logical parameters fall within their required ranges.
#'
#' @param phyloseq_obj A phyloseq object containing counts and sample data.
#'   If `NULL`, both `counts` and `covariates` must be provided.
#' @param counts Optional numeric matrix of raw counts with dimensions
#'   \emph{samples} \eqn{\times} \emph{taxa}. Required when \code{phyloseq_obj} is \code{NULL}.
#' @param covariates Optional data frame of sample-level covariates referenced in \code{formula}.
#'   Required when \code{phyloseq_obj} is \code{NULL}.
#' @param formula A model formula (e.g. `~ group + age`).
#' @param conf Numeric scalar in (0,1), confidence level for interval estimation.
#' @param lambda_beta Non-negative numeric scalar, penalty for regression betas.
#' @param zero_adj Logical flag; if `TRUE`, zero-adjustment is applied.
#' @param outlier_replace Logical flag; if `TRUE`, outliers are replaced.
#' @param outlier_conf Numeric scalar in (0,1), threshold for outlier detection.
#' @param beta_bound Numeric vector of length 2, increasing bounds for betas.
#' @param theta_bound Numeric vector of length 2, increasing positive bounds for dispersions.
#' @param eta_bound Numeric vector of length 2, increasing bounds for zero‐inflation.
#' @param max_iter Positive integer, maximum number of optimization iterations.
#' @param tol Positive numeric, convergence tolerance.
#'
#' @return Invisibly returns `TRUE` if all inputs pass validation; otherwise
#'   throws an error describing the violation.
#'
#' @noRd
check_DISSA_inputs = function(phyloseq_obj, counts, covariates, formula,
                              conf, lambda_beta, zero_adj, outlier_replace, outlier_conf,
                              beta_bound, theta_bound, eta_bound, max_iter, tol){

  # Check input source: phyloseq_obj OR (counts AND covariates) must be provided
  if(is.null(phyloseq_obj)){
    if(missing(counts) || missing(covariates)){
      stop("User should provide either a 'phyloseq_obj', or both 'counts' and 'covariates'.")
    }
    # Check that counts and covariates have matching sample size
    if(nrow(counts) != nrow(covariates)){
      stop("'counts' and 'covariates' must have the same number of rows (samples).")
    }
  }

  # Check formula is provided and correctly formatted
  if(missing(formula)){
    stop("A model 'formula' must be provided, e.g., ~ group + age.")
  }
  if(!inherits(formula, "formula")){
    stop("The 'formula' argument must be a valid formula object, e.g., ~ group + age.")
  }

  # Validate parameter values
  if(!is.numeric(conf) || length(conf) != 1 || conf <= 0 || conf >= 1){
    stop("'conf' must be a single numeric value between 0 and 1 (exclusive).")
  }
  if(!is.numeric(lambda_beta) || length(lambda_beta) != 1 || lambda_beta < 0){
    stop("'lambda_beta' must be a single non-negative numeric value.")
  }
  if(!is.logical(zero_adj) || length(zero_adj) != 1){
    stop("'zero_adj' must be a single logical value (TRUE or FALSE).")
  }
  if(!is.logical(outlier_replace) || length(outlier_replace) != 1){
    stop("'outlier_replace' must be a single logical value (TRUE or FALSE).")
  }
  if(!is.numeric(outlier_conf) || length(outlier_conf) != 1 || outlier_conf <= 0 || outlier_conf >= 1){
    stop("'outlier_conf' must be a single numeric value between 0 and 1 (exclusive).")
  }
  if(!is.numeric(beta_bound) || length(beta_bound) != 2 || beta_bound[1] >= beta_bound[2]){
    stop("'beta_bound' must be a numeric vector of length 2 with increasing values, e.g., c(-10, 10).")
  }
  if(!is.numeric(theta_bound) || length(theta_bound) != 2 || theta_bound[1] >= theta_bound[2] || any(theta_bound <= 0)){
    stop("'theta_bound' must be a positive numeric vector of length 2 with increasing values, e.g., c(0.01, 10).")
  }
  if(!is.numeric(eta_bound) || length(eta_bound) != 2 || eta_bound[1] >= eta_bound[2]){
    stop("'eta_bound' must be a numeric vector of length 2 with increasing values, e.g., c(-10, 10).")
  }
  if(!is.numeric(max_iter) || length(max_iter) != 1 || max_iter <= 0 || max_iter != floor(max_iter)){
    stop("'max_iter' must be a single positive integer.")
  }
  if(!is.numeric(tol) || length(tol) != 1 || tol <= 0){
    stop("'tol' must be a single positive numeric value.")
  }
}
