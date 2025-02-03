// compute_pvalues.cpp
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// A simple approximation of the trigamma function (the derivative of the digamma)
// using recurrence for small x and an asymptotic expansion for larger x.
double my_trigamma(double x) {
  double result = 0.0;
  // Use recurrence until x is sufficiently large.
  while (x < 5.0) {
    result += 1.0 / (x * x);
    x += 1.0;
  }
  double r = 1.0 / x;
  double r2 = r * r;
  // Asymptotic expansion: psi1(x) ~ 1/x + 1/(2x^2) + 1/(6x^3) - 1/(30x^5) + 1/(42x^7)
  result += 0.5 * r2 + (1.0 + r2 * (1.0/6.0 - r2 * (1.0/30.0 - r2/42.0)))/x;
  return result;
}

// [[Rcpp::export]]
NumericVector compute_pvalues(List param) {
  // Extract elements from the input list
  NumericMatrix data_R = param["data"];   // n x m count matrix
  NumericMatrix X_R = param["X"];           // n x p design matrix
  NumericMatrix beta_R = param["beta"];     // m x p coefficient matrix
  NumericVector theta_R = param["theta"];   // vector of length m (dispersion)
  NumericMatrix mu_R = param["mu"];         // n x m fitted means
  NumericMatrix tau_R = param["tau"];       // n x m posterior probabilities for zero inflation
  NumericMatrix v_R = param["v"];           // n x m indicator matrix for zeros
  double lambda = as<double>(param["lambda"]);  // penalty parameter
  
  // Dimensions (provided in the param list)
  int n = as<int>(param["n"]);  // number of samples
  int m = as<int>(param["m"]);  // number of features (e.g., genes)
  int p = as<int>(param["p"]);  // number of covariates
  
  // Convert R objects to Armadillo objects (without copying data)
  arma::mat data(data_R.begin(), n, m, false);
  arma::mat X(X_R.begin(), n, p, false);
  arma::mat beta(beta_R.begin(), m, p, false);
  arma::vec theta(theta_R.begin(), theta_R.size(), false);
  arma::mat mu(mu_R.begin(), n, m, false);
  arma::mat tau(tau_R.begin(), n, m, false);
  arma::mat v(v_R.begin(), n, m, false);
  
  // Allocate Fisher Information components
  arma::vec I_eta(n, fill::zeros);         // Will become an n x n diagonal block
  arma::vec I_theta(m, fill::zeros);         // Will become an m x m diagonal block
  
  // Cross-terms between theta and eta (m x n)
  arma::mat I_theta_eta(m, n, fill::zeros);
  
  // For beta: build a vector of m matrices (each p x p)
  std::vector<arma::mat> I_beta_list(m, arma::mat(p, p, fill::zeros));
  
  // Cross-terms between theta and beta (m x (m*p))
  arma::mat I_theta_beta(m, m * p, fill::zeros);
  
  // Cross-terms between eta and beta (n x (m*p))
  arma::mat I_eta_beta(n, m * p, fill::zeros);
  
  // --- Compute I_eta for each sample i ---
  for (int i = 0; i < n; i++) {
    double new_I_eta = 0.0;
    for (int j = 0; j < m; j++) {
      double d = data(i, j);
      double mu_ij = mu(i, j);
      double tau_ij = tau(i, j);
      double v_ij = v(i, j);
      double theta_j = theta[j];
      double denom = mu_ij + theta_j;
      double term = (d + (1 - tau_ij * v_ij) * theta_j) * mu_ij * theta_j / (denom * denom);
      new_I_eta += term;
    }
    I_eta[i] = new_I_eta;
  }
  
  // --- Compute Fisher Information for each feature (gene) j ---
  for (int j = 0; j < m; j++) {
    double new_I_theta = 0.0;
    arma::mat new_I_beta(p, p, fill::zeros);
    arma::vec new_I_theta_beta(p, fill::zeros);
    
    for (int i = 0; i < n; i++) {
      double d = data(i, j);
      double mu_ij = mu(i, j);
      double tau_ij = tau(i, j);
      double v_ij = v(i, j);
      double theta_j = theta[j];
      double denom = mu_ij + theta_j;
      
      // Use our own trigamma approximation instead of an external dependency
      double trigamma_term = my_trigamma(d + theta_j) - my_trigamma(theta_j);
      
      double term_theta = (d * theta_j + (1 - tau_ij * v_ij) * (mu_ij * mu_ij)) / (theta_j * (denom * denom)) + trigamma_term;
      new_I_theta += term_theta;
      
      double factor = (d + (1 - tau_ij * v_ij) * theta_j) * mu_ij * theta_j / (denom * denom);
      for (int a = 0; a < p; a++) {
        for (int b = 0; b < p; b++) {
          new_I_beta(a, b) += factor * X(i, a) * X(i, b);
        }
      }
      
      double factor_theta_beta = ((1 - tau_ij * v_ij) * mu_ij - d) * mu_ij / (denom * denom);
      for (int a = 0; a < p; a++) {
        new_I_theta_beta[a] += factor_theta_beta * X(i, a);
      }
      
      // Store the cross-term between eta and beta for gene j into the appropriate block
      for (int a = 0; a < p; a++) {
        I_eta_beta(i, j * p + a) = factor * X(i, a);
      }
      
      // Fill in the (j,i) element for the cross-term between theta and eta
      I_theta_eta(j, i) = ((1 - tau_ij * v_ij) * mu_ij - d) * mu_ij / (denom * denom);
    }
    // Adjust theta information with the penalty term
    I_theta[j] = new_I_theta - 6 * lambda / std::pow(theta[j], 4);
    // Save I_beta for gene j
    I_beta_list[j] = new_I_beta;
    // Place the cross-term between theta and beta into the appropriate block
    for (int a = 0; a < p; a++) {
      I_theta_beta(j, j * p + a) = new_I_theta_beta[a];
    }
  }
  
  // --- Assemble the full Fisher Information matrix ---
  arma::mat I_theta_mat = diagmat(I_theta);  // m x m
  arma::mat I_eta_mat = diagmat(I_eta);        // n x n
  
  // Build block-diagonal matrix for I_beta (dimensions: (m*p) x (m*p))
  arma::mat I_beta_mat(m * p, m * p, fill::zeros);
  for (int j = 0; j < m; j++) {
    I_beta_mat.submat(j * p, j * p, j * p + p - 1, j * p + p - 1) = I_beta_list[j];
  }
  
  // Total dimension of the full Fisher Information matrix:
  // Parameters are ordered as [theta (length m); eta (length n); beta (length m*p)]
  int total = m + n + m * p;
  arma::mat I_full(total, total, fill::zeros);
  
  // Fill in blocks:
  I_full.submat(0, 0, m - 1, m - 1) = I_theta_mat;                   // theta block
  I_full.submat(0, m, m - 1, m + n - 1) = I_theta_eta;                 // theta-eta cross
  I_full.submat(0, m + n, m - 1, total - 1) = I_theta_beta;            // theta-beta cross
  I_full.submat(m, 0, m + n - 1, m - 1) = I_theta_eta.t();              // eta-theta cross
  I_full.submat(m, m, m + n - 1, m + n - 1) = I_eta_mat;                // eta block
  I_full.submat(m, m + n, m + n - 1, total - 1) = I_eta_beta;           // eta-beta cross
  I_full.submat(m + n, 0, total - 1, m - 1) = I_theta_beta.t();         // beta-theta cross
  I_full.submat(m + n, m, total - 1, m + n - 1) = I_eta_beta.t();       // beta-eta cross
  I_full.submat(m + n, m + n, total - 1, total - 1) = I_beta_mat;         // beta block
  
  // Invert the full Fisher Information matrix to obtain the covariance matrix
  arma::mat cov_matrix = pinv(I_full);
  
  // Extract the covariance block corresponding to beta parameters.
  // Beta parameters occupy the last m*p rows and columns.
  arma::mat cov_beta = cov_matrix.submat(m + n, m + n, total - 1, total - 1);
  
  // --- Compute p-values for each gene ---
  // For each gene j, test the null hypothesis H0: beta[j,2] == 0.
  // The index for the second coefficient in the j-th block is (j - 1)*p + 2 in R's 1-indexing,
  // which corresponds to index j*p + 1 in 0-indexing.
  NumericVector pvalues(m);
  for (int j = 0; j < m; j++) {
    int index = j * p + 1;  // second coefficient for gene j (0-indexed)
    double var_beta = cov_beta(index, index);
    double beta_val = beta(j, 1);  // beta[j,2] in R (1-indexed) is column 1 here
    double test_stat = (beta_val * beta_val) / var_beta;
    double pval = 1.0 - R::pchisq(test_stat, 1.0, 1, 0);
    pvalues[j] = pval;
  }
  
  return pvalues;
}
