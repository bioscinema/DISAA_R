#' Perform the E-Step of the EM Algorithm
#'
#' The `e_step` function updates the `mu` and `tau` parameters during the Expectation-Maximization (EM) algorithm.
#' This step calculates the expected value of the latent variables based on current parameter estimates.
#'
#' @param mu A numeric matrix (n x m) of mean estimates for the count data.
#' @param tau A numeric matrix (n x m) of posterior probabilities for zero inflation.
#' @param data A numeric matrix (n x m) of observed count data.
#' @param X A numeric matrix (n x p) of covariates for the samples.
#' @param eta A numeric vector (n) of sample-specific offsets.
#' @param beta A numeric matrix (m x p) of covariate effects for each feature.
#' @param pi A numeric vector (m) of zero-inflation probabilities for each feature.
#' @param theta A numeric vector (m) of dispersion parameters for each feature.
#' @param v A binary matrix (n x m) indicating zero counts (1 if zero, 0 otherwise).
#'
#' @return A list containing updated values of `mu` and `tau`.
#' @examples
#' # Example usage:
#' # e_results <- e_step(mu, tau, data, X, eta, beta, pi, theta, v)
#' @export
e_step <- function(mu, tau, data, X, eta, beta, pi, theta, v) {
  for (i in seq_len(nrow(data))) {
    for (j in seq_len(ncol(data))) {
      mu[i, j] <- exp(eta[i] + sum(X[i, ] * beta[j, ]))
      tau[i, j] <- pi[j] / (pi[j] + (1 - pi[j]) * (theta[j] / (mu[i, j] + theta[j]))^theta[j])
    }
  }
  list(mu = mu, tau = tau)
}


#' Update Parameters During the M-Step
#'
#' The `update_parameters` function updates the `pi`, `theta`, and `beta` parameters during the Maximization step of the EM algorithm.
#' This step optimizes the likelihood of the model with respect to these parameters.
#'
#' @param data A numeric matrix (n x m) of observed count data.
#' @param X A numeric matrix (n x p) of covariates for the samples.
#' @param eta A numeric vector (n) of sample-specific offsets.
#' @param beta A numeric matrix (m x p) of covariate effects for each feature.
#' @param theta A numeric vector (m) of dispersion parameters for each feature.
#' @param pi A numeric vector (m) of zero-inflation probabilities for each feature.
#' @param tau A numeric matrix (n x m) of posterior probabilities for zero inflation.
#' @param v A binary matrix (n x m) indicating zero counts (1 if zero, 0 otherwise).
#' @param lambda A regularization parameter for controlling model complexity.
#' @param beta_bound A numeric vector specifying the bounds for `beta`.
#' @param theta_bound A numeric vector specifying the bounds for `theta`.
#'
#' @return A list containing updated values for `pi`, `theta`, and `beta`.
#' @examples
#' # Example usage:
#' # m_results <- update_parameters(data, X, eta, beta, theta, pi, tau, v, lambda, beta_bound, theta_bound)
#' @export

update_parameters <- function(data, X, eta, beta, theta, pi, tau, v, lambda, beta_bound, theta_bound) {
  for (j in seq_len(ncol(data))) {
    # Update pi
    pi[j] <- sum(tau[, j] * v[, j]) / nrow(data)

    # Update theta and beta
    par_j <- c(theta[j], beta[j, ])
    lj <- function(par_j) {
      theta_j <- par_j[1]
      beta_j <- par_j[-1]
      mu_j <- exp(eta + X %*% as.matrix(beta_j))

      lambda / theta_j^2 - sum(
        (1 - tau[, j] * v[, j]) * theta_j * log(theta_j) -
          ((1 - tau[, j] * v[, j]) * theta_j + data[, j]) * log(mu_j + theta_j) +
          data[, j] * log(mu_j) +
          lgamma(data[, j] + theta_j) -
          lgamma(theta_j)
      )
    }

    res <- nloptr::nloptr(
      x0 = as.numeric(par_j),
      eval_f = lj,
      lb = c(theta_bound[1], rep(beta_bound[1], length(beta[j, ]))),
      ub = c(theta_bound[2], rep(beta_bound[2], length(beta[j, ]))),
      opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1.0e-6)
    )

    theta[j] <- res$solution[1]
    beta[j, ] <- res$solution[-1]
  }
  list(pi = pi, theta = theta, beta = beta)
}

