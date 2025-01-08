#' Update Eta Values
#'
#' Updates the `eta` parameter for the MZINB model. This function uses
#' optimization to maximize the likelihood of the `eta` parameter for
#' each observation in the data. The `eta` parameter represents the
#' logarithm of the mean of non-zero counts, adjusted for covariates.
#'
#' @param data A matrix of count data with rows representing samples and
#'   columns representing features (genes, taxa, etc.).
#' @param X A matrix of covariates with rows corresponding to samples and
#'   columns representing covariates.
#' @param beta A matrix of regression coefficients for each feature and
#'   covariate.
#' @param theta A vector of dispersion parameters for each feature.
#' @param tau A matrix of posterior probabilities for each sample and feature
#'   indicating whether a value is zero-inflated.
#' @param v A binary matrix indicating whether a value is zero-inflated (1 if
#'   zero, 0 otherwise).
#' @param eta A vector of current `eta` values for each sample.
#' @param eta_bound The upper bound for the `eta` parameter during optimization.
#'
#' @return A vector of updated `eta` values for each sample.
#'
#' @details
#' The function performs the following steps for each sample:
#' - Defines a likelihood function for `eta`, based on the observed counts,
#'   covariates, regression coefficients, and dispersion parameters.
#' - Optimizes the likelihood function using the `nloptr` package with
#'   specified bounds on `eta`.
#' - Updates the `eta` parameter for the sample based on the optimization
#'   result.
#'
#' This function is designed to be used iteratively as part of the
#' Expectation-Maximization (EM) algorithm for the MZINB model.
#'
#' @examples
#' # Example usage
#' updated_eta <- update_eta(
#'   data = counts,
#'   X = covariates,
#'   beta = beta_matrix,
#'   theta = theta_vector,
#'   tau = tau_matrix,
#'   v = v_matrix,
#'   eta = eta_vector,
#'   eta_bound = 10
#' )
#'
#' @export
update_eta <- function(data, X, beta, theta, tau, v, eta, eta_bound) {
  n <- nrow(data)

  for (i in seq_len(n)) {
    # Define likelihood function for eta
    li <- function(eta_i) {
      mu_i <- exp(eta_i + beta %*% as.matrix(X[i, ]))
      -sum(-((1 - tau[i, ] * v[i, ]) * theta + data[i, ]) * log(mu_i + theta) + data[i, ] * log(mu_i))
    }

    # Optimize eta using nloptr
    res <- nloptr::nloptr(
      x0 = as.numeric(eta[i]),
      eval_f = li,
      lb = c(-10),  # Lower bound for eta
      ub = c(eta_bound),  # Upper bound for eta
      opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1.0e-6)
    )

    # Update eta for the current sample
    eta[i] <- res$solution[1]
  }

  eta
}
