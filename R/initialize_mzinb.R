#' Initialize Parameters for the MZINB Model
#'
#' The `initialize_mzinb` function initializes all parameters and data structures required for the MZINB model.
#' It sets up matrices and vectors for iterative optimization using the EM algorithm.
#'
#' @param counts A numeric matrix (n x m) of count data. Rows represent samples, and columns represent features.
#' @param covariates A numeric matrix (n x p) of covariates for the samples.
#' @param lambda A regularization parameter for controlling model complexity. Default is 0.5.
#' @param beta_bound A numeric vector specifying the lower and upper bounds for the `beta` parameter. Default is `c(-10, 10)`.
#' @param theta_bound A numeric vector specifying the lower and upper bounds for the `theta` parameter. Default is `c(1e-2, 10)`.
#' @param eta_bound A numeric value specifying the upper bound for the `eta` parameter. Default is `10`.
#'
#' @return A list containing initialized values for all parameters and variables required in the EM algorithm.
#' @examples
#' # Example usage:
#' # params <- initialize_mzinb(counts, covariates)
#' @export

initialize_mzinb <- function(counts, covariates, lambda = 0.5, beta_bound = c(-10, 10), theta_bound = c(1e-2, 10), eta_bound = 10) {
  data <- as.matrix(counts)
  X <- as.matrix(covariates)

  n <- dim(data)[1]
  m <- dim(data)[2]
  p <- dim(X)[2]

  pi <- rep(0.5, m)
  eta <- apply(data, 1, function(x) median(log(x[x != 0])))
  eta <- ifelse(is.na(eta), 0, eta)
  eta_bound <- max(max(eta), eta_bound)

  beta <- matrix(1, m, p)
  mu <- matrix(0, n, m)
  theta <- rep(1, m)
  v <- 1 * (data == 0)
  tau <- v

  list(data = data, X = X, n = n, m = m, p = p, pi = pi, eta = eta, eta_bound = eta_bound, beta = beta, mu = mu, theta = theta, v = v, tau = tau, lambda = lambda, beta_bound = beta_bound, theta_bound = theta_bound)
}

