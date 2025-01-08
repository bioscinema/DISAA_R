#' DISAA: MZINB Model for Analyzing Count Data
#'
#' The `disaa` function implements the Multivariate Zero-Inflated Negative Binomial (MZINB) model for analyzing count data,
#' particularly useful for datasets with a high prevalence of zero counts.
#' It leverages the Expectation-Maximization (EM) algorithm to iteratively optimize parameters, providing a flexible
#' and robust framework for modeling count data with covariates.
#'
#' @param counts A numeric matrix (n x m) of count data, where rows represent samples and columns represent features.
#'               It is expected to contain non-negative integers.
#' @param covariates A numeric matrix (n x p) of covariates for the samples, where rows correspond to samples and columns correspond to covariates.
#' @param lambda A regularization parameter for controlling model complexity. Default is 0.5.
#' @param beta_bound A numeric vector specifying the lower and upper bounds for the `beta` parameter. Default is `c(-10, 10)`.
#' @param theta_bound A numeric vector specifying the lower and upper bounds for the `theta` parameter. Default is `c(1e-2, 10)`.
#' @param eta_bound A numeric value specifying the upper bound for the `eta` parameter. Default is `10`.
#' @param max_iter An integer specifying the maximum number of iterations for the EM algorithm. Default is `100`.
#'
#' @return A list containing the following elements:
#' \item{pi}{A numeric vector of estimated zero-inflation probabilities for each feature.}
#' \item{theta}{A numeric vector of estimated dispersion parameters for each feature.}
#' \item{eta}{A numeric vector of sample-specific offsets.}
#' \item{beta}{A numeric matrix (m x p) of covariate effects for each feature.}
#' \item{mu}{A numeric matrix (n x m) of estimated mean values for the count data.}
#' \item{pvalues}{A numeric vector of adjusted p-values for the significance of covariate effects.}
#'
#' @examples
#' # Example usage:
#' # Assuming `counts` is an n x m count matrix and `covariates` is an n x p covariate matrix
#' # result <- disaa(counts, covariates)
#'
#' @details
#' The MZINB model assumes that the observed counts are generated from a mixture of zero-inflated and Negative Binomial (NB) components.
#' The model is suitable for analyzing high-dimensional count data with covariates, such as those from microbiome or RNA-seq studies.
#'
#' @export

#' Main MZINB Function
#'
#' This function performs the MZINB analysis using an Expectation-Maximization algorithm.
#'
#' @param counts A matrix of count data (n x m).
#' @param covariates A matrix of covariates for the samples (n x p).
#' @param lambda Regularization parameter for the model.
#' @param beta_bound Bounds for the `beta` parameter.
#' @param theta_bound Bounds for the `theta` parameter.
#' @param eta_bound Bound for the `eta` parameter.
#' @param max_iter The maximum number of iterations for the EM algorithm.
#' @return A list containing optimized values for `pi`, `theta`, `eta`, `beta`, `mu`, and computed p-values.
#' @export
disaa <- function(counts, covariates, lambda = 0.5, beta_bound = c(-10, 10), theta_bound = c(1e-2, 10), eta_bound = 10, max_iter = 100) {
  params <- initialize_mzinb(counts, covariates, lambda, beta_bound, theta_bound, eta_bound)

  for (iter in seq_len(max_iter)) {
    e_results <- e_step(params$mu, params$tau, params$data, params$X, params$eta, params$beta, params$pi, params$theta, params$v)
    params$mu <- e_results$mu
    params$tau <- e_results$tau

    m_results <- update_parameters(params$data, params$X, params$eta, params$beta, params$theta, params$pi, params$tau, params$v, params$lambda, params$beta_bound, params$theta_bound)
    params$pi <- m_results$pi
    params$theta <- m_results$theta
    params$beta <- m_results$beta

    params$eta <- update_eta(params)
  }

  params$mu <- compute_mu(params$data, params$X, params$eta, params$beta)

  pvalues <- compute_pvalues(params$data, params$beta, params$theta, params$mu, params$tau, params$v, params$lambda)

  return(list(pi = params$pi, theta = params$theta, eta = params$eta, beta = params$beta, mu = params$mu, pvalues = pvalues))
}


