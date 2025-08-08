#' EM Optimization for DISSA Model (Internal)
#'
#' Performs the Expectation-Maximization (EM) algorithm to estimate parameters of a
#' zero-inflated negative binomial model for count data in the DISSA pipeline.
#'
#' @param data Numeric matrix of counts (n samples × m taxa).
#' @param X Numeric design matrix of covariates (n samples × p predictors).
#' @param eta Numeric vector of sample-specific intercepts (length n).
#' @param v Binary matrix indicating observed zeros (n × m), where 1 denotes a zero count.
#' @param tau Numeric matrix of posterior probabilities (n × m), initialized by the user.
#' @param beta Numeric regression coefficient matrix (m taxa × p predictors), initial values.
#' @param theta Numeric vector of dispersion parameters (length m), initial values.
#' @param pi Numeric vector of zero-inflation probabilities (length m), initial values.
#' @param mu Numeric matrix of expected counts (n × m), initial values.
#' @param lambda_beta Numeric penalty parameter for `beta` shrinkage.
#' @param beta_bound Numeric vector of length 2 with lower and upper bounds for `beta`.
#' @param theta_bound Numeric vector of length 2 with lower and upper bounds for `theta`.
#' @param eta_bound Numeric vector of length 2 with lower and upper bounds for `eta`.
#' @param update_index Integer vector of taxa indices to update; defaults to all taxa.
#' @param max_iter Positive integer, maximum number of EM iterations.
#' @param tol Positive numeric, convergence tolerance on relative log-likelihood change.
#'
#' @return A list with components:
#'   \describe{
#'     \item{beta}{Updated regression coefficients (m × p matrix).}
#'     \item{theta}{Estimated dispersion parameters (length m).}
#'     \item{pi}{Estimated zero-inflation probabilities (length m).}
#'     \item{eta}{Estimated sample-specific intercepts (length n).}
#'     \item{tau}{Final posterior zero-probabilities (n × m matrix).}
#'     \item{mu}{Final expected counts (n × m matrix).}
#'     \item{mu_theta}{Empirical Bayes prior mean for log(dispersion).}
#'     \item{sigma2_theta}{Empirical Bayes prior variance for log(dispersion).}
#'   }
#'
#' @noRd
run_em_optimization = function(data, X, eta, v, tau, beta, theta, pi, mu,
                               lambda_beta, beta_bound, theta_bound, eta_bound,
                               update_index = NULL, max_iter, tol){
  n = nrow(data)
  m = ncol(data)
  p = ncol(X)
  old_likelihood = 1

  # initialize EB hyperparameters for theta only
  mu_theta = 0          # prior mean for log(theta)
  sigma2_theta = 1      # prior variance for log(theta)

  if(is.null(update_index)) update_index = seq_len(m)

  for(iter in 1:max_iter){
    ### E-step
    # for(i in 1:n){
    #   eta_i = eta[i]
    #   X_i = X[i, ]
    #   for(j in update_index){
    #     mu[i, j] = exp(eta_i + sum(X_i * beta[j, ]))
    #     denom = pi[j] + (1 - pi[j]) * (theta[j] / (mu[i, j] + theta[j]))^theta[j]
    #     tau[i, j] = pi[j] / denom
    #   }
    # }
    # --- E-STEP (vectorized) ---
    mu  <- exp( outer(eta, rep(1, m)) + X %*% t(beta) )

    # Build column-constant matrices for theta, pi
    Theta <- matrix(theta, n, m, byrow = TRUE)
    Pi    <- matrix(pi,    n, m, byrow = TRUE)

    # denom_ij = pi_j + (1 - pi_j) * (theta_j/(mu_ij + theta_j))^{theta_j}
    base  <- Theta / (mu + Theta)
    denom <- Pi + (1 - Pi) * exp( sweep(log(base), 2, theta, `*`) )  # exp(theta_j * log(base_ij))
    tau   <- Pi / denom

    ### M-step for beta & theta
    for(j in update_index){
      tau_v_j = tau[, j] * v[, j]
      # pi[j] = sum(tau_v_j) / n
      pi = colMeans(tau*v)

      par_j = c(theta[j], beta[j, ])
      lj = function(par_j){
        theta_j = par_j[1]
        beta_j = par_j[-1]
        mu_j = exp(eta + X %*% beta_j)
        # fixed penalty for beta, EB penalty for theta
        penalty_beta = lambda_beta * sum(beta_j[-1]^2)
        penalty_theta = ((log(theta_j) - mu_theta)^2) / (2 * sigma2_theta)
        neg_log_lik = -sum(
          (1 - tau_v_j) * theta_j * log(theta_j) -
            ((1 - tau_v_j) * theta_j + data[, j]) * log(mu_j + theta_j) +
            data[, j] * log(mu_j) + lgamma(data[, j] + theta_j) - lgamma(theta_j)
        )
        penalty_beta + penalty_theta + neg_log_lik
      }

      res = nloptr(
        x0 = par_j, eval_f = lj,
        lb = c(theta_bound[1], rep(beta_bound[1], p)),
        ub = c(theta_bound[2], rep(beta_bound[2], p)),
        opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = tol)
      )

      theta[j] = res$solution[1]
      beta[j, ] = res$solution[-1]
    }

    ### M-step for eta and tau update
    new_likelihood = 0
    for(i in 1:n){
      li = function(eta_i){
        mu_i = exp(eta_i + beta %*% X[i, ])
        -sum(-((1 - tau[i, ] * v[i, ]) * theta + data[i, ]) * log(mu_i + theta) + data[i, ] * log(mu_i))
      }
      res = nloptr(
        x0 = eta[i], eval_f = li,
        lb = eta_bound[1], ub = eta_bound[2],
        opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = tol)
      )
      eta[i] = res$solution
      new_likelihood = new_likelihood + res$objective

      eta_i = eta[i]
      X_i = X[i, ]
      for(j in update_index){
        mu[i, j] = exp(eta_i + sum(X_i * beta[j, ]))
        denom = pi[j] + (1 - pi[j]) * (theta[j] / (mu[i, j] + theta[j]))^theta[j]
        tau[i, j] = pi[j] / denom
      }
    }

    ### EB hyperparameter update for theta only
    logs = log(theta[update_index])
    mu_theta = mean(logs)
    sigma2_theta = mean((logs - mu_theta)^2)

    if(abs((new_likelihood - old_likelihood) / old_likelihood) < tol) break
    old_likelihood <- new_likelihood
  }

  return(list(beta = beta, theta = theta, pi = pi, eta = eta, tau = tau, mu = mu, mu_theta = mu_theta, sigma2_theta = sigma2_theta))
}
