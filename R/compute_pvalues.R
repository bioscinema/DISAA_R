#' Compute Wald Test p-values for Each Covariate Across Taxa
#'
#' This function constructs and inverts the observed Fisher information matrix
#' for a log-normal model and computes Wald test p-values for each covariate
#' (excluding intercept and optionally a zero-adjust term) across all taxa.
#'
#' @param data A numeric matrix of count data (samples × taxa).
#' @param X A numeric design matrix for the covariates (samples × covariates).
#' @param mu A numeric matrix of expected means (same shape as `data`).
#' @param tau A numeric matrix of zero-inflation probabilities (samples × taxa).
#' @param v A numeric matrix of indicator variables for zero-inflation (samples × taxa).
#' @param theta A numeric vector of dispersion parameters for each taxon.
#' @param beta A numeric matrix of regression coefficients (taxa × covariates).
#' @param lambda A numeric value controlling the regularization strength on `theta`.
#' @param zero_adj A logical or numeric indicator specifying whether a zero-adjustment variable is present in `X`.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{I}{The full observed Fisher information matrix.}
#'   \item{cov}{The inverse of the Fisher information matrix (covariance matrix).}
#'   \item{pvalues_all}{A matrix of adjusted p-values for each covariate across all taxa.}
#' }
#'
#' @details
#' The function computes second derivatives (Hessian blocks) for `theta`, `eta`, and `beta`,
#' as well as cross-derivatives among them to build the full information matrix.
#' Wald statistics are calculated using the squared ratio of estimated beta coefficients
#' to their variances, and p-values are adjusted using the Benjamini-Hochberg method.
#'
#' @importFrom Matrix bdiag
#' @importFrom MASS ginv
#' @importFrom stats pchisq p.adjust
#'
#' @keywords internal
#' @noRd
compute_single_var_pvalues = function(data, X, mu, tau, v, theta, beta, lambda, zero_adj){
  
  n = nrow(data)
  m = ncol(data)
  p = ncol(X)
  
  I_theta = numeric(m)
  I_eta = numeric(n)
  I_theta_eta = matrix(0, m, n)
  I_beta = new("dgCMatrix", Dim = as.integer(c(0, 0)))  
  I_theta_beta = matrix(0, m, m * p)
  I_eta_beta = matrix(0, n, m * p)
  
  for(i in 1:n){
    new_I_eta = 0
    for(j in 1:m){
      temp = (data[i, j] + (1 - tau[i, j] * v[i, j]) * theta[j])
      new_I_eta = new_I_eta + (temp * mu[i, j] * theta[j]) / (mu[i, j] + theta[j])^2
    }
    I_eta[i] = new_I_eta
  }
  
  for(j in 1:m){
    new_I_theta = 0
    new_I_beta = matrix(0, p, p)
    new_I_theta_beta = numeric(p)
    
    for(i in 1:n){
      mu_ij = mu[i, j]
      data_ij = data[i, j]
      tau_ij = tau[i, j]
      v_ij = v[i, j]
      X_i = X[i, ]
      
      denom = (mu_ij + theta[j])^2
      common_term = data_ij + (1 - tau_ij * v_ij) * theta[j]
      
      ### Diagonal Hessian block for theta
      new_I_theta = new_I_theta -
        (data_ij * theta[j] + (1 - tau_ij * v_ij) * mu_ij^2) / (theta[j] * denom) +
        trigamma(data_ij + theta[j]) - trigamma(theta[j])
      
      ### Hessian block for beta
      weight = (common_term * mu_ij * theta[j]) / denom
      new_I_beta = new_I_beta + weight * (X_i %*% t(X_i))
      
      ### Cross derivative between theta and beta
      cross_term = ((1 - tau_ij * v_ij) * mu_ij - data_ij) * mu_ij / denom
      new_I_theta_beta = new_I_theta_beta + cross_term * X_i
      
      ### Cross derivative between eta and beta
      I_eta_beta[i, ((j - 1) * p + 1):(j * p)] = weight * X_i
      
      ### Cross derivative between theta and eta
      I_theta_eta[j, i] = cross_term
    }
    
    I_theta[j] = new_I_theta + 6 * lambda / theta[j]^4
    I_beta = as.matrix(bdiag(I_beta, new_I_beta))
    I_theta_beta[j, ((j - 1) * p + 1):(j * p)] = new_I_theta_beta
  }
  
  ### Convert I_theta and I_eta to diagonal matrices
  I_theta = diag(I_theta)
  I_eta = diag(I_eta)
  
  ### Combine all into full observed information matrix
  I = rbind(
    cbind(I_theta, I_theta_eta, I_theta_beta),
    cbind(t(I_theta_eta), I_eta, I_eta_beta),
    cbind(t(I_theta_beta), t(I_eta_beta), I_beta)
  )
  
  ### Extract beta-related submatrix
  cov = ginv(I)
  cov_beta = cov[(nrow(cov) - m * p + 1):nrow(cov), (nrow(cov) - m * p + 1):nrow(cov)]
  
  ### Initialize p-value matrix
  adjusted_start = 2 + as.numeric(zero_adj)
  pvalue_cols = p - 1 - as.numeric(zero_adj)
  pvalues_all = matrix(0, m, pvalue_cols)
  
  ### Compute p-values for each covariate (excluding intercept and zero-adjust column if applicable)
  for(var in adjusted_start:p){
    var_idx = var - 1 - as.numeric(zero_adj)
    
    ### Vectorized computation of p-values for all taxa
    beta_vec = beta[, var]
    cov_diag = sapply(1:m, function(j) cov_beta[p * (j - 1) + var, p * (j - 1) + var])
    
    chi_sq = beta_vec^2 / cov_diag
    raw_pvals = 1 - pchisq(chi_sq, df = 1)
    
    ### Adjust for multiple testing
    pvalues_all[, var_idx] = stats::p.adjust(raw_pvals, method = "BH")
  }
  
  return(list(I = I, cov = cov, pvalues_all = pvalues_all))
}