#' Run EM Optimization for Zero-Inflated Negative Binomial Model
#'
#' This function performs an Expectation-Maximization (EM) optimization for a 
#' zero-inflated negative binomial (ZINB) model with sample-level and taxa-level covariates.
#'
#' @param data A matrix of observed count data (samples × taxa).
#' @param X A design matrix of sample-level covariates (samples × covariates).
#' @param eta A numeric vector of sample-level intercepts.
#' @param v A binary matrix indicating observed non-zero support (samples × taxa).
#' @param tau A numeric matrix of latent dropout probabilities (samples × taxa).
#' @param beta A matrix of taxa-level covariate effects (taxa × covariates).
#' @param theta A numeric vector of dispersion parameters (length = number of taxa).
#' @param pi A numeric vector of zero-inflation probabilities (length = number of taxa).
#' @param mu A matrix of mean parameters (samples × taxa).
#' @param lambda Regularization penalty for dispersion parameters.
#' @param beta_bound A numeric vector of lower and upper bounds for beta coefficients (length 2).
#' @param theta_bound A numeric vector of lower and upper bounds for theta (length 2).
#' @param eta_bound A numeric vector of lower and upper bounds for eta (length 2).
#' @param update_index Optional integer vector indicating which taxa to update (default: all taxa).
#' @param max_iter Maximum number of EM iterations.
#' @param tol Tolerance for convergence of the likelihood.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{beta}{Estimated beta coefficients (taxa × covariates).}
#'   \item{theta}{Estimated dispersion parameters for each taxon.}
#'   \item{pi}{Estimated zero-inflation probabilities.}
#'   \item{eta}{Updated intercept vector.}
#'   \item{tau}{Updated dropout probability matrix.}
#'   \item{mu}{Estimated mean matrix.}
#' }
#'
#' @importFrom nloptr nloptr
#' @export
run_em_optimization <- function(data, X, eta, v, tau, beta, theta, pi, mu, 
                               lambda, beta_bound, theta_bound, eta_bound, 
                               update_index = NULL, max_iter, tol){
  
  n = nrow(data)
  m = ncol(data)
  p = ncol(X)
  old_likelihood = 1
  
  ### If no specific taxa indices are provided, update all taxa
  if(is.null(update_index)){
    update_index = 1:m
  }
  
  for(iter in 1:max_iter){
    ### E-step
    for(i in 1:n){
      eta_i = eta[i]
      X_i = X[i, ]
      for(j in update_index){
        mu[i, j] = exp(eta_i + sum(X_i * beta[j, ]))
        denom = pi[j] + (1 - pi[j]) * (theta[j] / (mu[i, j] + theta[j]))^theta[j]
        tau[i, j] = pi[j] / denom
      }
    }
    
    ### M-step
    for(j in update_index){
      tau_v_j = tau[, j] * v[, j]
      pi[j] = sum(tau_v_j) / n
      
      par_j = c(theta[j], beta[j, ])
      
      lj = function(par_j){
        theta_j = par_j[1]
        beta_j = par_j[-1]
        mu_j = exp(eta + X %*% beta_j)
        penalty = lambda / theta_j^2
        neg_log_lik = -sum(
          (1 - tau_v_j) * theta_j * log(theta_j) -
            ((1 - tau_v_j) * theta_j + data[, j]) * log(mu_j + theta_j) +
            data[, j] * log(mu_j) + lgamma(data[, j] + theta_j) - lgamma(theta_j)
        )
        return(penalty + neg_log_lik)
      }
      
      res = nloptr(
        x0 = par_j, eval_f = lj,
        lb = c(theta_bound[1], rep(beta_bound[1], p)),
        ub = c(theta_bound[2], rep(beta_bound[2], p)),
        opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = tol)
      )
      
      theta[j] = res$solution[1]
      beta[j, ] = res$solution[-1]
    }
    
    ### Update eta and recompute tau
    new_likelihood = 0
    for(i in 1:n){
      li = function(eta_i){
        mu_i = exp(eta_i + beta %*% X[i, ])
        -sum(-((1 - tau[i, ] * v[i, ]) * theta + data[i, ]) * log(mu_i + theta) + data[i, ] * log(mu_i))
      }
      
      res = nloptr(
        x0 = eta[i], eval_f = li,
        lb = eta_bound[1], ub = eta_bound[2],
        opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = tol)
      )
      
      eta[i] = res$solution
      new_likelihood = new_likelihood + res$objective
      
      ### Update tau after eta update
      eta_i = eta[i]
      X_i = X[i, ]
      for(j in update_index){
        mu[i, j] = exp(eta_i + sum(X_i * beta[j, ]))
        denom = pi[j] + (1 - pi[j]) * (theta[j] / (mu[i, j] + theta[j]))^theta[j]
        tau[i, j] = pi[j] / denom
      }
    }
    
    if(abs((new_likelihood - old_likelihood) / old_likelihood) < tol && iter > max_iter){
      break
    }else{
      old_likelihood = new_likelihood
    }
  }
  
  return(list(beta = beta, theta = theta, pi = pi, eta = eta, tau = tau, mu = mu))
}