#' Compute P-Values for Beta Coefficients
#'
#' The `compute_pvalues` function computes the p-values for the `beta` coefficients using the Fisher Information Matrix.
#'
#' @import Rcpp
#' @import RcppArmadillo
#' @import MASS
#' @import Matrix
#' @useDynLib DISAA, .registration = TRUE
#'
#' @param data A numeric matrix (n x m) of observed count data.
#' @param beta A numeric matrix (m x p) of covariate effects for each feature.
#' @param theta A numeric vector (length m) of dispersion parameters for each feature.
#' @param mu A numeric matrix (n x m) of computed mean values for the count data.
#' @param tau A numeric matrix (n x m) of posterior probabilities for zero inflation.
#' @param v A binary matrix (n x m) indicating zero counts (1 if zero, 0 otherwise).
#' @param lambda A regularization parameter for controlling model complexity.
#'
#' @return A numeric vector of adjusted p-values for the `beta` coefficients.
#'
#' @export
compute_pvalues <- function(data, beta, theta, mu, tau, v, lambda) {
  m <- ncol(data)  # Number of features
  p <- ncol(beta)  # Number of covariates

  # Initialize structures for Fisher Information Matrix components
  I_beta <- Matrix::Diagonal(0)  # Sparse diagonal matrix
  I_theta <- numeric(m)
  I_theta_beta <- matrix(0, nrow = m, ncol = m * p)

  # Compute Fisher Information Matrix components
  for (j in seq_len(m)) {
    I_beta_j <- matrix(0, nrow = p, ncol = p)
    I_theta_j <- 0
    I_theta_beta_j <- matrix(0, nrow = p, ncol = 1)

    for (i in seq_len(nrow(data))) {
      mu_ij <- mu[i, j]
      theta_j <- theta[j]
      tau_ij <- tau[i, j]
      v_ij <- v[i, j]
      x_i <- as.matrix(beta[j, ])

      # Update Fisher Information components
      I_beta_j <- I_beta_j + ((data[i, j] + (1 - tau_ij * v_ij) * theta_j) * mu_ij * theta_j / (mu_ij + theta_j)^2) * x_i %*% t(x_i)
      I_theta_j <- I_theta_j + (data[i, j] * theta_j + (1 - tau_ij * v_ij) * mu_ij^2) / (theta_j * (mu_ij + theta_j)^2)
      I_theta_beta_j <- I_theta_beta_j + (((1 - tau_ij * v_ij) * mu_ij - data[i, j]) * mu_ij / (mu_ij + theta_j)^2) * x_i
    }

    I_theta[j] <- I_theta_j - 6 * lambda / theta[j]^4
    I_beta <- as.matrix(Matrix::bdiag(I_beta, I_beta_j))
    I_theta_beta[j, ((j - 1) * p + 1):(j * p)] <- I_theta_beta_j
  }

  # Assemble the Fisher Information Matrix
  I_theta <- diag(I_theta)
  I <- rbind(
    cbind(I_theta, I_theta_beta),
    cbind(t(I_theta_beta), I_beta)
  )

  # Invert the Fisher Information Matrix to get the covariance matrix
  cov <- MASS::ginv(I)

  # Extract covariance of beta coefficients
  cov_beta <- cov[(nrow(cov) - m * p + 1):nrow(cov), (nrow(cov) - m * p + 1):ncol(cov)]

  # Compute p-values for each beta coefficient
  pvalues <- numeric(m)
  for (j in seq_len(m)) {
    pvalues[j] <- 1 - stats::pchisq(beta[j, 2]^2 / cov_beta[2 * j, 2 * j], df = 1)
  }

  # Adjust p-values for multiple testing
  pvalues <- stats::p.adjust(pvalues, method = "BH")

  return(pvalues)
}


#' Compute Mu for MZINB Model
#'
#' This function computes the `mu` matrix based on the current values of `eta` and `beta`.
#'
#' @param data The count data matrix.
#' @param X The covariates matrix.
#' @param eta The `eta` vector.
#' @param beta The `beta` matrix.
#' @return The computed `mu` matrix.
compute_mu <- function(data, X, eta, beta) {
  n <- nrow(data)
  m <- ncol(data)
  mu <- matrix(0, n, m)

  for (i in seq_len(n)) {
    for (j in seq_len(m)) {
      mu[i, j] <- exp(eta[i] + sum(X[i, ] * beta[j, ]))
    }
  }

  return(mu)
}

