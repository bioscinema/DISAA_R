#include <RcppArmadillo.h>
#include <cmath> // For math functions like pow and log

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec compute_pvalues(const arma::mat& data,
                          const arma::mat& beta,
                          const arma::vec& theta,
                          const arma::mat& mu,
                          const arma::mat& tau,
                          const arma::mat& v,
                          double lambda) {
  int m = data.n_cols;  // Number of features
  int p = beta.n_cols;  // Number of covariates

  // Initialize structures for Fisher Information Matrix components
  arma::sp_mat I_beta(p * m, p * m);  // Sparse matrix for I_beta
  arma::vec I_theta(m, fill::zeros);
  arma::mat I_theta_beta(m, p * m, fill::zeros);

  // Compute Fisher Information Matrix components
  for (int j = 0; j < m; ++j) {
    arma::mat I_beta_j(p, p, fill::zeros);
    double I_theta_j = 0.0;
    arma::vec I_theta_beta_j(p, fill::zeros);

    for (int i = 0; i < data.n_rows; ++i) {
      double mu_ij = mu(i, j);
      double theta_j = theta[j];
      double tau_ij = tau(i, j);
      double v_ij = v(i, j);
      arma::rowvec x_i = beta.row(j);

      // Update Fisher Information components
      I_beta_j += ((data(i, j) + (1 - tau_ij * v_ij) * theta_j) * mu_ij * theta_j / std::pow(mu_ij + theta_j, 2)) * (x_i.t() * x_i);
      I_theta_j += (data(i, j) * theta_j + (1 - tau_ij * v_ij) * std::pow(mu_ij, 2)) / (theta_j * std::pow(mu_ij + theta_j, 2));
      I_theta_beta_j += (((1 - tau_ij * v_ij) * mu_ij - data(i, j)) * mu_ij / std::pow(mu_ij + theta_j, 2)) * x_i.t();
    }

    I_theta[j] = I_theta_j - 6 * lambda / std::pow(theta[j], 4);
    I_beta.submat(j * p, j * p, (j + 1) * p - 1, (j + 1) * p - 1) = I_beta_j;
    I_theta_beta.submat(j, j * p, j, (j + 1) * p - 1) = I_theta_beta_j.t();
  }

  // Assemble the Fisher Information Matrix
  arma::mat I_theta_diag = arma::diagmat(I_theta);
  arma::mat I_upper = join_horiz(I_theta_diag, I_theta_beta);
  arma::mat I_lower = join_horiz(I_theta_beta.t(), arma::mat(I_beta));  // Convert sparse matrix to dense
  arma::mat I = join_vert(I_upper, I_lower);

  // Invert the Fisher Information Matrix to get the covariance matrix
  arma::mat cov = pinv(I);

  // Extract covariance of beta coefficients
  arma::mat cov_beta = cov.submat(m, m, cov.n_rows - 1, cov.n_cols - 1);

  // Compute p-values for each beta coefficient
  arma::vec pvalues(m, fill::zeros);
  for (int j = 0; j < m; ++j) {
    double stat = std::pow(beta(j, 1), 2) / cov_beta(2 * j, 2 * j);  // Beta[j, 2]^2
    pvalues[j] = 1 - R::pchisq(stat, 1, 1, 0);
  }

  // Adjust p-values for multiple testing using BH method
  Function p_adjust("p.adjust");
  pvalues = as<arma::vec>(p_adjust(Named("p", pvalues), Named("method", "BH")));

  return pvalues;
}
