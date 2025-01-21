#' DISAA: MZINB Model for Analyzing Count Data
#'
#' The `disaa` function implements the Multivariate Zero-Inflated Negative Binomial (MZINB) model for analyzing high-dimensional count data.
#' It uses the Expectation-Maximization (EM) algorithm to iteratively optimize model parameters, allowing the analysis of datasets with a high prevalence of zero counts.
#'
#' @param counts A numeric matrix (n x m) of count data, where rows represent samples and columns represent features. The matrix is expected to contain non-negative integers.
#' @param covariates A numeric matrix (n x p) of covariates for the samples, where rows correspond to samples and columns correspond to covariates.
#' @param lambda A regularization parameter for controlling model complexity. The default value is 0.5.
#' @param beta_bound A numeric vector specifying the lower and upper bounds for the `beta` parameter. The default is `c(-10, 10)`.
#' @param theta_bound A numeric vector specifying the lower and upper bounds for the `theta` parameter. The default is `c(1e-2, 10)`.
#' @param eta_bound A numeric value specifying the upper bound for the `eta` parameter. The default is 10.
#' @param max_iter An integer specifying the maximum number of iterations for the EM algorithm. The default is 100.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{pi}{A numeric vector of estimated zero-inflation probabilities for each feature.}
#'   \item{theta}{A numeric vector of estimated dispersion parameters for each feature.}
#'   \item{eta}{A numeric vector of sample-specific offsets.}
#'   \item{beta}{A numeric matrix (m x p) of covariate effects for each feature.}
#'   \item{mu}{A numeric matrix (n x m) of estimated mean values for the count data.}
#'   \item{pvalues}{A numeric vector of adjusted p-values for the significance of covariate effects.}
#' }
#'
#' @details
#' The `disaa` function iteratively updates the model parameters using the Expectation-Maximization (EM) algorithm:
#' \itemize{
#'   \item **E-Step**: Computes the expected values of the latent variables based on current parameter estimates.
#'   \item **M-Step**: Updates the model parameters to maximize the likelihood of the observed data.
#' }
#' The model is suitable for analyzing count data with covariates, such as those from microbiome or RNA-seq studies, which often have a high proportion of zero counts.
#'
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


