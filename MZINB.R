library(nloptr)
library(Matrix)
library(MASS)

MZINB = function(counts, covariates, lambda = 0.5,
                 beta_bound = c(-10, 10), theta_bound = c(1e-2, 10), eta_bound = 10, max_iter = 100){
  
  ### Initialization
  data = as.matrix(counts)
  row.names(data) = NULL
  colnames(data) = NULL
  X = covariates
  row.names(X) = NULL
  colnames(X) = NULL
  n = dim(data)[1]
  m = dim(data)[2]
  p = dim(X)[2]
  pi = rep(0.5, m)
  eta = apply(data, 1, function(x) median(log(x[x != 0])))
  eta = ifelse(is.na(eta), 0, eta)
  eta_bound = max(max(eta), eta_bound)
  beta = matrix(1, m, p)
  mu = matrix(0, n, m)
  theta = rep(1, m)
  v = 1 * (data == 0)
  tau = v
  
  ### Optimize
  for(iter in 1:max_iter){
    ### E step
    for(i in 1:n){
      for(j in 1:m){
        mu[i, j] = exp(eta[i] + sum(X[i, ] * beta[j, ]))
        tau[i, j] = pi[j] / (pi[j] + (1 - pi[j]) * (theta[j] / (mu[i, j] + theta[j])) ^ theta[j])
      }
    }
    
    for(j in 1:m){
      
      ### Update pi
      pi[j] = sum(tau[, j] * v[, j]) / n
      
      ### Update theta and beta
      par_j = c(theta[j], beta[j, ])
      
      lj = function(par_j){
        theta_j = par_j[1]
        beta_j = par_j[2:length(par_j)]
        mu_j = exp(eta + X %*% as.matrix(beta_j))
        
        return(lambda / theta_j^2 -sum((1 - tau[, j] * v[, j]) * theta_j * log(theta_j) - ((1 - tau[, j] * v[, j]) * theta_j + data[, j]) * log(mu_j + theta_j) +  data[, j] * log(mu_j) + lgamma(data[, j] + theta_j) - lgamma(theta_j)))
      }
      
      res = nloptr(x0 = as.numeric(par_j), eval_f = lj, lb = c(theta_bound[1], rep(beta_bound[1], p)), ub = c(theta_bound[2], rep(beta_bound[2], p)),
                   opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1.0e-6))
      
      theta[j] = res$solution[1]
      beta[j, ] = res$solution[2:length(par_j)]
      
    }
    
    ### Update eta
    for(i in 1:n){
      for(j in 1:m){
        mu[i, j] = exp(eta[i] + sum(X[i, ] * beta[j, ]))
        tau[i, j] = pi[j] / (pi[j] + (1 - pi[j]) * (theta[j] / (mu[i, j] + theta[j])) ^ theta[j])
      }
    }
    for(i in 1:n){
      
      li = function(eta_i){
        mu_i = exp(eta_i + beta %*% as.matrix(X[i, ]))
        return(-sum(-((1 - tau[i, ] * v[i, ]) * theta + data[i, ]) * log(mu_i + theta) + data[i, ] * log(mu_i)))
      }
      
      res = nloptr(x0 = as.numeric(eta[i]), eval_f = li, lb = c(-10), ub = c(eta_bound),
                   opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1.0e-6))
      
      eta[i] = res$solution[1]
      
    }
  }
  
  ### Compute mu
  for(i in 1:n){
    for(j in 1:m){
      mu[i, j] = exp(eta[i] + sum(X[i, ] * beta[j, ]))
    }
  }
  
  ### Compute pvalues
  I_theta = c()
  I_eta = c()
  I_theta_eta = matrix(0, m, n)
  I_beta = new("dgCMatrix", Dim = as.integer(c(0, 0)))
  I_theta_beta = matrix(0, m, m*p)
  I_eta_beta = matrix(0, n, m*p)
  
  for(i in 1:n){
    new_I_eta = 0
    for(j in 1:m){
      new_I_eta = new_I_eta + ((data[i, j] + (1 - tau[i, j] * v[i, j]) * theta[j]) * mu[i, j] * theta[j] / (mu[i, j] + theta[j]) ^ 2)
    }
    I_eta = c(I_eta, new_I_eta) 
  }
  for(j in 1:m){
    new_I_theta = 0
    new_I_beta = matrix(0, p, p)
    new_I_theta_beta = matrix(0, p, 1)
    for(i in 1:n){
      new_I_theta = new_I_theta + (data[i, j] * theta[j] + (1 - tau[i, j] * v[i, j]) * mu[i, j] ^ 2) / (theta[j] * (mu[i, j] + theta[j]) ^ 2) + trigamma(data[i, j] + theta[j]) - trigamma(theta[j])
      new_I_beta = new_I_beta + ((data[i, j] + (1 - tau[i, j] * v[i, j]) * theta[j]) * mu[i, j] * theta[j] / (mu[i, j] + theta[j]) ^ 2) * X[i, ] %*% t(X[i, ])
      new_I_theta_beta = new_I_theta_beta + (((1 - tau[i, j] * v[i, j]) * mu[i, j] - data[i, j]) * mu[i, j] / (mu[i, j] + theta[j]) ^ 2) * X[i, ]
      I_eta_beta[i, ((j-1)*p+1):(j*p)] =  ((data[i, j] + (1 - tau[i, j] * v[i, j]) * theta[j]) * mu[i, j] * theta[j] / (mu[i, j] + theta[j]) ^ 2) * X[i, ]
      I_theta_eta[j, i] = ((1 - tau[i, j] * v[i, j]) * mu[i, j] - data[i, j]) * mu[i, j] / (mu[i, j] + theta[j]) ^ 2
    }
    I_theta = c(I_theta, new_I_theta - 6 * lambda / theta[j]^4)
    I_beta = as.matrix(bdiag(I_beta, new_I_beta))
    I_theta_beta[j, ((j-1)*p+1):(j*p)] = new_I_theta_beta
  }
  I_theta = diag(I_theta)
  I_eta = diag(I_eta)
  
  I = rbind(
    cbind(I_theta, I_theta_eta, I_theta_beta),
    cbind(t(I_theta_eta), I_eta, I_eta_beta),
    cbind(t(I_theta_beta), t(I_eta_beta), I_beta)
  )
  cov = ginv(I)
  cov_beta = cov[(dim(cov)[1] - m*p + 1):dim(cov)[1], (dim(cov)[1] - m*p + 1):dim(cov)[1]]
  pvalues = c()
  for(j in 1:m){
    pvalues = c(pvalues, 1 - pchisq(beta[j, 2] ^ 2 / cov_beta[2*j, 2*j], 1))
  }
  pvalues = stats::p.adjust(pvalues, "BH")
  
  return(list(pi = pi, theta = theta, eta = eta, beta = beta, mu = mu, pvalues = pvalues))
  
}


N = 100 # number of samples
M = 200 # number of genes

pi = 0.5
theta = 1
pi_real = rep(pi, M)
theta_real = rep(theta, M)

### Generate eta
eta_real = runif(N, min = -3, max = 3)

### Generate beta, 30% significant
beta_real = cbind(rep(1, M), rbinom(N, 2, 0.3))

### Generate X
X = cbind(rep(1, N), rbinom(N, 1, 0.5)) 

### Generate zero-inflated negative binomial
counts = matrix(0, N, M)
mu = matrix(0, N, M)
for(i in 1:N){
  for(j in 1:M){
    mu_ij = exp(eta_real[i] + 1 + sum(X[i, ] * beta_real[j, ]))
    mu[i, j] = mu_ij
    if(runif(1) <= 1 - pi_real[j]) counts[i, j] = rnbinom(1, mu = mu_ij, size = theta_real[j])
  }
}

### Fit
out = MZINB(counts, X, lambda = 0.1, max_iter = 30)

### Results
mean(out$pi)
mean(out$theta)
mean(out$pvalues < 0.05)
