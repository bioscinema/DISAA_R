#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List compute_pvalues_cpp(List params) {
  // Extract values from params list
  int n = params["n"];
  int m = params["m"];
  int p = params["p"];
  arma::mat X = as<arma::mat>(params["X"]);
  arma::mat beta = as<arma::mat>(params["beta"]);
  arma::vec theta = as<arma::vec>(params["theta"]);
  arma::mat tau = as<arma::mat>(params["tau"]);
  arma::mat v = as<arma::mat>(params["v"]);
  arma::mat data = as<arma::mat>(params["data"]);
  arma::mat mu = as<arma::mat>(params["mu"]);

  // Initializing Fisher information elements
  arma::vec I_theta(m, fill::zeros);
  arma::vec I_eta(n, fill::zeros);
  arma::mat I_theta_eta(m, n, fill::zeros);
  arma::mat I_theta_beta(m, m * p, fill::zeros);
  arma::mat I_eta_beta(n, m * p, fill::zeros);
  arma::mat I_beta = arma::zeros<arma::mat>(m * p, m * p);

  // Compute the Fisher information elements
  for (int i = 0; i < n; ++i) {
    double new_I_eta = 0.0;

    for (int j = 0; j < m; ++j) {
      double mu_ij = mu(i, j);
      double theta_j = theta[j];
      double data_ij = data(i, j);
      double tau_v = (1.0 - tau(i, j) * v(i, j));

      new_I_eta += (data_ij + tau_v * theta_j) * mu_ij * theta_j / std::pow(mu_ij + theta_j, 2);
    }

    I_eta[i] = new_I_eta;
  }

  for (int j = 0; j < m; ++j) {
    double new_I_theta = 0.0;
    arma::mat new_I_beta = arma::zeros<arma::mat>(p, p);
    arma::vec new_I_theta_beta = arma::zeros<arma::vec>(p);

    for (int i = 0; i < n; ++i) {
      double mu_ij = mu(i, j);
      double theta_j = theta[j];
      double data_ij = data(i, j);
      double tau_v = (1.0 - tau(i, j) * v(i, j));
      double denom = std::pow(mu_ij + theta_j, 2);

      new_I_theta += (data_ij * theta_j + tau_v * mu_ij * mu_ij) / (theta_j * denom) +
        R::trigamma(data_ij + theta_j) - R::trigamma(theta_j);

      new_I_beta += (data_ij + tau_v * theta_j) * mu_ij * theta_j / denom * X.row(i).t() * X.row(i);

      new_I_theta_beta += ((tau_v * mu_ij - data_ij) * mu_ij / denom) * X.row(i).t();

      I_eta_beta(i, span(j * p, (j + 1) * p - 1)) = ((data_ij + tau_v * theta_j) * mu_ij * theta_j / denom) * X.row(i).t();
      I_theta_eta(j, i) = ((tau_v * mu_ij - data_ij) * mu_ij / denom);
    }

    I_theta[j] = new_I_theta;
    I_beta.submat(j * p, j * p, (j + 1) * p - 1, (j + 1) * p - 1) = new_I_beta;
    I_theta_beta(j, span(j * p, (j + 1) * p - 1)) = new_I_theta_beta.t();
  }

  // Construct Fisher Information matrix
  arma::mat I_theta_diag = diagmat(I_theta);
  arma::mat I_eta_diag = diagmat(I_eta);

  arma::mat I_top = join_horiz(join_horiz(I_theta_diag, I_theta_eta), I_theta_beta);
  arma::mat I_mid = join_horiz(join_horiz(I_theta_eta.t(), I_eta_diag), I_eta_beta);
  arma::mat I_bot = join_horiz(join_horiz(I_theta_beta.t(), I_eta_beta.t()), I_beta);

  arma::mat I = join_vert(join_vert(I_top, I_mid), I_bot);

  // Calculate the covariance matrix
  arma::mat cov = pinv(I); // Using pseudo-inverse for numerical stability
  arma::mat cov_beta = cov.submat((cov.n_rows - m * p), (cov.n_cols - m * p), cov.n_rows - 1, cov.n_cols - 1);

  // Compute p-values
  arma::vec pvalues(m, fill::zeros);
  for (int j = 0; j < m; ++j) {
    pvalues[j] = 1.0 - R::pchisq(std::pow(beta(j, 1), 2) / cov_beta(2 * j, 2 * j), 1, 1, 0);
  }

  // Update params with p-values
  params["pvalues"] = pvalues;

  return params;
}

