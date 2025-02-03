#include <RcppArmadillo.h>
#include <cmath> // for std::exp, std::log, std::lgamma

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Function to update the eta parameter for each sample
// [[Rcpp::export]]
arma::vec update_eta(const arma::mat& data, const arma::mat& X, const arma::mat& beta,
                     const arma::vec& theta, const arma::mat& tau, const arma::mat& v,
                     arma::vec eta, double eta_bound) {

  int n = data.n_rows;  // Number of samples
  int m = data.n_cols;  // Number of features

  for (int i = 0; i < n; ++i) {
    // Define the likelihood function for eta_i
    auto likelihood_function = [&](double eta_i) {
      arma::vec mu_i = exp(eta_i + X.row(i) * beta.t());
      double neg_log_likelihood = 0.0;

      for (int j = 0; j < m; ++j) {
        double term1 = -((1 - tau(i, j) * v(i, j)) * theta[j] + data(i, j)) * log(mu_i[j] + theta[j]);
        double term2 = data(i, j) * log(mu_i[j]);
        neg_log_likelihood += term1 + term2;
      }
      Rcpp::Rcout << "Negative log-likelihood: " << neg_log_likelihood << std::endl;
      return -neg_log_likelihood;  // Return negative log-likelihood
    };

    // Optimize eta using a simple gradient descent approach
    double eta_i = eta[i];
    double step_size = 0.01;
    int max_iter = 100;

    for (int iter = 0; iter < max_iter; ++iter) {
      double grad = (likelihood_function(eta_i + step_size) - likelihood_function(eta_i)) / step_size;
      eta_i -= step_size * grad;

      // Apply bounds
      if (eta_i > eta_bound) eta_i = eta_bound;
      if (eta_i < -10) eta_i = -10;

      // Convergence check
      if (std::abs(grad) < 1e-6) break;
    }

    // Update eta
    eta[i] = eta_i;
  }

  return eta;
}
