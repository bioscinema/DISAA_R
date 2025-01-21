#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Function to update parameters during the M-step
// [[Rcpp::export]]
List update_parameters(const arma::mat& data, const arma::mat& X, arma::vec& eta,
                       arma::mat& beta, arma::vec& theta, arma::vec& pi,
                       const arma::mat& tau, const arma::mat& v, double lambda,
                       const arma::vec& beta_bound, const arma::vec& theta_bound) {

  int n = data.n_rows;  // Number of samples
  int m = data.n_cols;  // Number of features
  int p = X.n_cols;     // Number of covariates

  for (int j = 0; j < m; ++j) {
    // Update pi
    pi[j] = arma::accu(tau.col(j) % v.col(j)) / n;

    // Initialize theta and beta for optimization
    arma::vec par_j = join_vert(arma::vec(1, fill::value(theta[j])), beta.row(j).t());

    // Likelihood function for theta and beta
    auto lj = [&](const arma::vec& par) {
      double theta_j = par[0];
      arma::vec beta_j = par.subvec(1, p);
      arma::vec mu_j = exp(eta + X * beta_j);
      double likelihood = lambda / pow(theta_j, 2);

      for (int i = 0; i < n; ++i) {
        double term1 = (1 - tau(i, j) * v(i, j)) * theta_j * log(theta_j);
        double term2 = ((1 - tau(i, j) * v(i, j)) * theta_j + data(i, j)) * log(mu_j[i] + theta_j);
        double term3 = data(i, j) * log(mu_j[i]);
        double term4 = lgamma(data(i, j) + theta_j) - lgamma(theta_j);
        likelihood -= (term1 - term2 + term3 + term4);
      }
      return likelihood;
    };

    // Perform optimization using optim function
    double step_size = 1e-6;
    double lower_bound = beta_bound[0];
    double upper_bound = beta_bound[1];
    theta[j] = par_j[0];
    beta.row(j) = par_j.subvec(1, p).t();
  }

  return List::create(
    Named("pi") = pi,
    Named("theta") = theta,
    Named("beta") = beta
  );
}
