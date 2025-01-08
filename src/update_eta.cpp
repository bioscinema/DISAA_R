#include <RcppArmadillo.h>
#include <Rmath.h> // For log() functions

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Define the objective function for optimization
double eta_objective(double eta_i, const arma::rowvec& X_i, const arma::mat& beta, const arma::vec& theta, const arma::rowvec& tau_i, const arma::rowvec& v_i, const arma::rowvec& data_i) {
  arma::vec mu_i = exp(eta_i + beta * X_i.t()); // Calculating mu for given eta_i
  double sum_val = 0.0;

  for (arma::uword j = 0; j < data_i.n_elem; ++j) {
    sum_val += -(((1 - tau_i[j] * v_i[j]) * theta[j] + data_i[j]) * log(mu_i[j] + theta[j]) - data_i[j] * log(mu_i[j]));
  }

  return sum_val;
}

// [[Rcpp::export]]
List update_eta_cpp(List params) {
  int n = params["n"];
  int m = params["m"];
  int p = params["p"];
  arma::mat X = as<arma::mat>(params["X"]);
  arma::mat beta = as<arma::mat>(params["beta"]);
  arma::vec theta = as<arma::vec>(params["theta"]);
  arma::mat tau = as<arma::mat>(params["tau"]);
  arma::mat v = as<arma::mat>(params["v"]);
  arma::mat data = as<arma::mat>(params["data"]);
  arma::vec eta = as<arma::vec>(params["eta"]);

  // Loop over each row in the dataset to update eta[i]
  for (int i = 0; i < n; ++i) {
    arma::rowvec X_i = X.row(i);
    arma::rowvec tau_i = tau.row(i);
    arma::rowvec v_i = v.row(i);
    arma::rowvec data_i = data.row(i);

    // Initial value for eta_i
    double eta_i = eta[i];

    // Optimization - using a basic gradient-free approach (simple step-based optimization)
    double lb = -10;
    double ub = 10;
    double best_eta = eta_i;
    double min_val = eta_objective(eta_i, X_i, beta, theta, tau_i, v_i, data_i);

    // Simple optimization by searching within bounds with a small step
    for (double test_eta = lb; test_eta <= ub; test_eta += 0.1) {
      double obj_val = eta_objective(test_eta, X_i, beta, theta, tau_i, v_i, data_i);
      if (obj_val < min_val) {
        min_val = obj_val;
        best_eta = test_eta;
      }
    }

    // Update eta[i]
    eta[i] = best_eta;
  }

  // Update the params list with the new eta values
  params["eta"] = eta;

  return params;
}
