#' Update Eta Values for Each Sample
#'
#' The `update_eta` function optimizes the sample-specific `eta` parameter using the `optim()` function from the `stats` package.
#' This parameter represents the logarithm of the mean of non-zero counts, adjusted for covariates.
#'
#' @param params A list containing:
#'   - `data`: A matrix of count data (n x m) where rows represent samples and columns represent features.
#'   - `X`: A matrix of covariates (n x p).
#'   - `beta`: A matrix of regression coefficients (m x p).
#'   - `theta`: A vector of dispersion parameters (m).
#'   - `tau`: A matrix of posterior probabilities for zero inflation (n x m).
#'   - `v`: A binary matrix indicating zero counts (n x m).
#'   - `eta`: A vector of current `eta` values (n).
#'   - `eta_bound`: The upper bound for `eta` during optimization.
#'
#' @return A numeric vector of updated `eta` values for each sample.
#'
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#' @importFrom nloptr nloptr
#' @importFrom stats optim
#' @useDynLib DISAA
#'
#'
#' @export
update_eta <- function(params) {
  data <- params$data
  X <- params$X
  beta <- params$beta
  theta <- params$theta
  tau <- params$tau
  v <- params$v
  eta <- params$eta
  eta_bound <- params$eta_bound

  n <- nrow(data)  # Number of samples

  for (i in seq_len(n)) {
    li <- function(eta_i) {
      mu_i <- exp(eta_i + as.numeric(beta %*% as.matrix(X[i, ])))
      return(-sum(-((1 - tau[i, ] * v[i, ]) * theta + data[i, ]) * log(mu_i + theta) + data[i, ] * log(mu_i)))
    }

    res = nloptr(x0 = as.numeric(eta[i]), eval_f = li, lb = c(-10), ub = c(eta_bound),
                 opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1.0e-6))

    eta[i] <- res$solution[1]
  }

  return(eta)
}
