#' Compute P-Values for Beta Coefficients
#'
#' The `compute_pvalues` function computes the p-values for the `beta` coefficients using the Fisher Information Matrix.
#' It extracts the necessary model components from a list of parameters and then constructs, inverts, and uses the Fisher
#' Information Matrix to test the significance of the regression coefficients.
#'
#' @import Rcpp
#' @import RcppArmadillo
#' @import MASS
#' @import Matrix
#' @useDynLib DISAA, .registration = TRUE
#'
#' @param param A list containing the following elements:
#' \describe{
#'   \item{data}{A numeric matrix (n x m) of observed count data.}
#'   \item{X}{A numeric matrix (n x p) representing the design matrix.}
#'   \item{beta}{A numeric matrix (m x p) of covariate effects for each feature.}
#'   \item{theta}{A numeric vector (length m) of dispersion parameters for each feature.}
#'   \item{mu}{A numeric matrix (n x m) of computed mean values for the count data.}
#'   \item{tau}{A numeric matrix (n x m) of posterior probabilities for zero inflation.}
#'   \item{v}{A binary matrix (n x m) indicating zero counts (1 if zero, 0 otherwise).}
#'   \item{lambda}{A regularization parameter controlling model complexity.}
#'   \item{n}{The number of samples.}
#'   \item{m}{The number of features (e.g., genes).}
#'   \item{p}{The number of covariates.}
#' }
#'
#' @return A numeric vector of adjusted p-values for the `beta` coefficients.
#'
#' @export
compute_pvalues <- function(param) {
  # Extract elements from param list
  data   <- param$data      # n x m count matrix
  X      <- param$X         # n x p design matrix
  beta   <- param$beta      # m x p coefficient matrix
  theta  <- param$theta     # vector of length m (per-feature dispersion)
  mu     <- param$mu        # n x m fitted means
  tau    <- param$tau       # n x m weight matrix
  v      <- param$v         # n x m indicator matrix for zeros
  lambda <- param$lambda    # penalty parameter
  
  # Dimensions
  n <- param$n
  m <- param$m
  p <- param$p
  
  #### Initialize Fisher Information components ####
  # For feature-level parameter theta (length m) and sample-level parameter eta (length n)
  I_theta <- numeric(m)       # will become an m x m diagonal block
  I_eta   <- numeric(n)       # will become an n x n diagonal block
  
  # Cross-terms between theta and eta (m x n)
  I_theta_eta <- matrix(0, nrow = m, ncol = n)
  
  # For beta (regression coefficients): we will build a list of m blocks (each p x p)
  I_beta_list <- vector("list", m)
  
  # Cross-terms between theta and beta (m x (m*p))
  I_theta_beta <- matrix(0, nrow = m, ncol = m * p)
  
  # Cross-terms between eta and beta (n x (m*p))
  I_eta_beta <- matrix(0, nrow = n, ncol = m * p)
  
  #### Compute I_eta (information for eta) for each sample i ####
  for (i in 1:n) {
    new_I_eta <- 0
    for (j in 1:m) {
      new_I_eta <- new_I_eta +
        ((data[i, j] + (1 - tau[i, j] * v[i, j]) * theta[j]) *
           mu[i, j] * theta[j] / (mu[i, j] + theta[j])^2)
    }
    I_eta[i] <- new_I_eta
  }
  
  #### Compute feature-level Fisher Information for each gene j ####
  for (j in 1:m) {
    new_I_theta    <- 0
    new_I_beta     <- matrix(0, nrow = p, ncol = p)
    new_I_theta_beta <- matrix(0, nrow = p, ncol = 1)
    
    for (i in 1:n) {
      # Information for theta (gene j): note the addition of the trigamma terms
      new_I_theta <- new_I_theta +
        ((data[i, j] * theta[j] + (1 - tau[i, j] * v[i, j]) * mu[i, j]^2) /
           (theta[j] * (mu[i, j] + theta[j])^2)) +
        (trigamma(data[i, j] + theta[j]) - trigamma(theta[j]))
      
      # Information for beta (gene j): use the covariate row X[i, ]
      new_I_beta <- new_I_beta +
        ((data[i, j] + (1 - tau[i, j] * v[i, j]) * theta[j]) *
           mu[i, j] * theta[j] / (mu[i, j] + theta[j])^2) *
        (X[i, ] %*% t(X[i, ]))
      
      # Cross-term between theta and beta for gene j
      new_I_theta_beta <- new_I_theta_beta +
        (((1 - tau[i, j] * v[i, j]) * mu[i, j] - data[i, j]) *
           mu[i, j] / (mu[i, j] + theta[j])^2) *
        X[i, ]
      
      # For each sample i, store the cross-term between eta and beta in the block for gene j.
      I_eta_beta[i, ((j - 1) * p + 1):(j * p)] <-
        ((data[i, j] + (1 - tau[i, j] * v[i, j]) * theta[j]) *
           mu[i, j] * theta[j] / (mu[i, j] + theta[j])^2) * X[i, ]
      
      # Also fill in the (j,i) element for the cross-term between theta and eta
      I_theta_eta[j, i] <- ((1 - tau[i, j] * v[i, j]) * mu[i, j] - data[i, j]) *
        mu[i, j] / (mu[i, j] + theta[j])^2
    }
    # Adjust theta information with the penalty term
    I_theta[j] <- new_I_theta - 6 * lambda / theta[j]^4
    # Save the beta information for gene j in the list
    I_beta_list[[j]] <- new_I_beta
    # Place the cross-term between theta and beta into the proper columns for gene j
    I_theta_beta[j, ((j - 1) * p + 1):(j * p)] <- as.vector(new_I_theta_beta)
  }
  
  #### Assemble full Fisher Information matrix ####
  # Convert I_theta and I_eta into diagonal matrices.
  I_theta_mat <- diag(I_theta, nrow = m, ncol = m)
  I_eta_mat   <- diag(I_eta, nrow = n, ncol = n)
  # Combine the list of I_beta blocks into a block-diagonal matrix.
  I_beta_mat <- as.matrix(Matrix::bdiag(I_beta_list))
  
  # Order the parameters as: first θ (length m), then η (length n), then β (length m*p).
  # Assemble the full Fisher information matrix accordingly.
  I_full <- rbind(
    cbind(I_theta_mat, I_theta_eta, I_theta_beta),
    cbind(t(I_theta_eta), I_eta_mat, I_eta_beta),
    cbind(t(I_theta_beta), t(I_eta_beta), I_beta_mat)
  )
  
  # Invert the full Fisher information matrix to get the covariance matrix.
  cov_matrix <- MASS::ginv(I_full)
  
  # Extract the covariance block corresponding to β.
  # The β parameters occupy the last m*p rows and columns.
  cov_beta <- cov_matrix[(m + n + 1):(m + n + m * p), (m + n + 1):(m + n + m * p)]
  
  #### Compute p-values for each gene ####
  # For each gene j, test the null hypothesis H0: beta[j,2] == 0.
  # The index for the second coefficient in the j-th block is (j - 1)*p + 2.
  pvalues <- numeric(m)
  for (j in 1:m) {
    index <- (j - 1) * p + 2  # index in the block for gene j
    test_stat <- (beta[j, 2]^2) / cov_beta[index, index]
    pvalues[j] <- 1 - stats::pchisq(test_stat, df = 1)
  }
  # Adjust the p-values for multiple testing (using the Benjamini–Hochberg method)
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

