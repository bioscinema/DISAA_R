#' Perform the E-Step of the EM Algorithm
#'
#' The `e_step` function performs the Expectation step (E-step) of the Expectation-Maximization (EM) algorithm for the Multivariate Zero-Inflated Negative Binomial (MZINB) model.
#' In this step, the function computes the expected values of the latent variables (`mu` and `tau`) based on the current parameter estimates.
#'
#' @param mu A numeric matrix (n x m) representing the mean values of the count data. This matrix will be updated during the E-step.
#' @param tau A numeric matrix (n x m) representing the posterior probabilities for zero inflation. This matrix will be updated during the E-step.
#' @param data A numeric matrix (n x m) of observed count data, where rows represent samples and columns represent features (e.g., genes or taxa).
#' @param X A numeric matrix (n x p) of covariates for the samples, where rows correspond to samples and columns represent covariates.
#' @param eta A numeric vector (length n) of sample-specific offsets.
#' @param beta A numeric matrix (m x p) of covariate effects for each feature.
#' @param pi A numeric vector (length m) of zero-inflation probabilities for each feature.
#' @param theta A numeric vector (length m) of dispersion parameters for each feature.
#' @param v A binary matrix (n x m) indicating zero counts (1 if zero, 0 otherwise).
#'
#' @return A list containing the updated values for `mu` and `tau`:
#' \describe{
#'   \item{mu}{A numeric matrix (n x m) of updated mean values for the count data.}
#'   \item{tau}{A numeric matrix (n x m) of updated posterior probabilities for zero inflation.}
#' }
#'
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



#' Update Parameters During the M-Step (Rcpp Version)
#'
#' This function updates the zero-inflation probabilities (\code{pi}), dispersion parameters (\code{theta}),
#' and regression coefficients (\code{beta}) during the Maximization (M) step of the EM algorithm for the
#' Multivariate Zero-Inflated Negative Binomial (MZINB) model. It leverages Rcpp for improved performance.
#'
#' @param data A numeric matrix of observed count data (n x m), where rows correspond to samples and columns to features.
#' @param X A numeric matrix (n x p) of covariates.
#' @param eta A numeric vector of length \code{n} representing sample-specific offsets.
#' @param beta A numeric matrix (m x p) of regression coefficients.
#' @param theta A numeric vector of length \code{m} representing dispersion parameters.
#' @param pi A numeric vector of length \code{m} containing the current zero-inflation probabilities.
#' @param tau A numeric matrix (n x m) of posterior probabilities for zero inflation.
#' @param v A binary matrix (n x m) indicating zero counts (1 if zero, 0 otherwise).
#' @param lambda A numeric regularization parameter controlling model complexity.
#' @param beta_bound A numeric vector of length 2 specifying the lower and upper bounds for the \code{beta} parameters.
#' @param theta_bound A numeric vector of length 2 specifying the lower and upper bounds for \code{theta}.
#' @param p An integer representing the number of covariates.
#'
#' @return A list containing updated values for \code{pi}, \code{theta}, and \code{beta}.
#'
#' @details
#' This function is part of the M-step of the EM algorithm for the MZINB model. It optimizes the negative log-likelihood
#' function for each feature (gene) using NLopt (via the \code{nloptr} package) to update the parameters.
#' This version is implemented in R but is intended to be used in conjunction with Rcpp-generated code for improved efficiency.
#'
#' @import Rcpp
#' @importFrom nloptr nloptr
#' @useDynLib DISAA, .registration = TRUE
#'
#' @export
update_parameters <- function(data, X, eta, beta, theta, pi, tau, v, lambda, beta_bound, theta_bound, p) {
  for (j in seq_len(ncol(data))) {
    # Update pi
    pi[j] <- sum(tau[, j] * v[, j]) / nrow(data)
    
    # Update theta and beta
    par_j <- c(theta[j], beta[j, ])
    
    # Define the log-likelihood function for optimization
    lj <- function(par_j) {
      theta_j <- par_j[1]
      beta_j <- par_j[-1]
      mu_j <- exp(eta + X %*% as.matrix(beta_j))
      
      # Calculate the negative log-likelihood
      return(lambda / theta_j^2 - sum(
        (1 - tau[, j] * v[, j]) * theta_j * log(theta_j) -
          ((1 - tau[, j] * v[, j]) * theta_j + data[, j]) * log(mu_j + theta_j) +
          data[, j] * log(mu_j) +
          lgamma(data[, j] + theta_j) - lgamma(theta_j)
      ))
    }
    
    # Use nloptr for optimization
    res = nloptr(x0 = as.numeric(par_j), eval_f = lj,
                 lb = c(theta_bound[1], rep(beta_bound[1], p)),
                 ub = c(theta_bound[2], rep(beta_bound[2], p)),
                 opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1.0e-6))
    
    # Update theta and beta
    theta[j] <- res$solution[1]
    beta[j, ] <- res$solution[-1]
  }
  
  # Return the updated values
  list(pi = pi, theta = theta, beta = beta)
}




