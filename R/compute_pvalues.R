#' Compute P-values for DISSA Model Parameters (Internal)
#'
#' Internal function to compute observed information matrix and pvalues.
#'
#' @param data Numeric matrix of observed counts (n samples × m taxa).
#' @param X Numeric design matrix of covariates (n samples × p predictors).
#' @param mu Numeric matrix of expected counts (n × m) from the EM fit.
#' @param tau Numeric matrix of posterior zero-probabilities (n × m) from the EM fit.
#' @param v Binary matrix (n × m) indicating zero observations (1 = zero count).
#' @param theta Numeric vector of dispersion parameters (length m).
#' @param beta Numeric matrix of regression coefficients (m taxa × p predictors).
#' @param lambda_beta Numeric penalty parameter for beta regularization.
#' @param zero_adj Logical flag indicating whether zero-adjustment was applied.
#' @param formula Model formula object specifying the effects tested (e.g., ~ group).
#' @param mu_theta Numeric prior mean for log-dispersion (from EB update).
#' @param sigma2_theta Numeric prior variance for log-dispersion (from EB update).
#' @param multi_class_factors Character vector of factor variables with >2 levels.
#' @param multi_class_dummy_list Named list mapping factor names to dummy column indices in X.
#' @param X_colnames Character vector of column names for X (length p).
#' @param conf Numeric scalar in (0,1), confidence level for p-value adjustment.
#' @param alpha Numeric scalar in (0,1), target false discovery rate for multiple testing.
#'
#' @return A list with components:
#'   \describe{
#'     \item{I}{Observed information matrix (dimensions (m+n+mp)²).}
#'     \item{cov}{Generalized inverse of I, the covariance matrix of all parameters.}
#'     \item{pvalues_all}{Matrix of adjusted p-values for each tested contrast (columns) across taxa (rows).}
#'   }
#'
#'
#' @noRd
compute_pvalues = function(data, X, mu, tau, v, theta, beta, lambda_beta, zero_adj, formula,
                           mu_theta, sigma2_theta, multi_class_factors, multi_class_dummy_list, X_colnames,
                           conf = 0.95, alpha = 0.1){

  n = nrow(data)
  m = ncol(data)
  p = ncol(X)

  I_theta = numeric(m)
  I_eta = numeric(n)
  I_theta_eta = matrix(0, m, n)
  I_beta = new("dgCMatrix", Dim = as.integer(c(0, 0)))
  I_theta_beta = matrix(0, m, m * p)
  I_eta_beta = matrix(0, n, m * p)

  # for(i in 1:n){
  #   new_I_eta = 0
  #   for(j in 1:m){
  #     temp = (data[i, j] + (1 - tau[i, j] * v[i, j]) * theta[j])
  #     new_I_eta = new_I_eta + (temp * mu[i, j] * theta[j]) / (mu[i, j] + theta[j])^2
  #   }
  #   I_eta[i] = new_I_eta
  # }
  # assume data, tau, v, mu are n x m; theta is length m
  Theta <- matrix(theta, nrow = nrow(data), ncol = ncol(data), byrow = TRUE)

  temp <- data + (1 - tau * v) * Theta
  num  <- temp * mu * Theta
  den  <- (mu + Theta)^2

  I_eta <- rowSums(num / den)


  for(j in 1:m){
    new_I_theta = 0
    new_I_beta = matrix(2 * lambda_beta, p, p)
    new_I_beta[1, 1] = 0
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

    I_theta[j] = new_I_theta + 1 / sigma2_theta * (1 - (log(theta[j]) - mu_theta)) / (theta[j]^2)
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

  ### G-inverse
  cov = ginv(I)

  ### Get variable index for contrast matrix, by default give all others vs reference
  C_index = list()
  pvalues_name = c()
  for(var in all.vars(formula)){
    if(var %in% multi_class_factors){

      var_index = which(multi_class_dummy_list[[var]])
      for(var_i in 2:length(var_index)){
        var_name = paste0(X_colnames[var_index[var_i]], " vs ", X_colnames[var_index[1]])
        pvalues_name = c(pvalues_name, var_name)
        C_index[[var_name]] = c(var_index[var_i], var_index[1])
      }

    } else {
      C_index[[var]] = which(var == X_colnames)
      pvalues_name = c(pvalues_name, var)
    }
  }

  ### Initialize p-value matrix
  pvalues_all = NULL

  for(var_name in pvalues_name){
    var_index = C_index[[var_name]]
    if(length(var_index) == 1){
      C_beta = 1
      C_cov = matrix(0, nrow = 1, ncol = p)
      C_cov[1, var_index] = 1
    }else{
      C_beta = c(1, -1)
      C_cov = matrix(0, nrow = 1, ncol = p)
      C_cov[1, var_index[1]] = -1
      C_cov[1, var_index[2]] = 1
    }
    pvals = vapply(1:m, function(j) {

      beta_tmp = sum(beta[j, var_index] * C_beta)
      C = cbind(matrix(0, nrow = 1, ncol = m + n + p * (j - 1)), C_cov, matrix(0, nrow = 1, ncol = p * (m - j)))
      stat = beta_tmp^2 / (C %*% cov %*% t(C))
      1 - pf(stat, df1 = 1, df2 = n - p)
    }, numeric(1))
    pvalues_all = cbind(pvalues_all, pvalues_adjust(pvals, conf, alpha))

  }

  colnames(pvalues_all) = pvalues_name

  return(list(I = I, cov = cov, pvalues_all = pvalues_all))
}
