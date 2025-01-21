#' Initialize Parameters for the MZINB Model
#'
#' The `initialize_mzinb` function sets up the initial parameters for the Multivariate Zero-Inflated Negative Binomial (MZINB) model.
#' It converts the count and covariate data into matrices and initializes values for `pi`, `eta`, `beta`, `mu`, `theta`, `v`, and `tau`.
#'
#' @param counts A numeric matrix (n x m) of count data, where rows represent samples and columns represent features.
#' @param covariates A numeric matrix (n x p) of covariates for the samples, where rows correspond to samples and columns correspond to covariates.
#' @param lambda A regularization parameter for controlling model complexity. Default is `0.5`.
#' @param beta_bound A numeric vector specifying the lower and upper bounds for the `beta` parameter. Default is `c(-10, 10)`.
#' @param theta_bound A numeric vector specifying the lower and upper bounds for the `theta` parameter. Default is `c(1e-2, 10)`.
#' @param eta_bound A numeric value specifying the upper bound for the `eta` parameter. Default is `10`.
#'
#' @return A list containing initialized parameters:
#'   - `data`: Count data matrix.
#'   - `X`: Covariates matrix.
#'   - `n`: Number of samples.
#'   - `m`: Number of features.
#'   - `p`: Number of covariates.
#'   - `pi`: Initial zero-inflation probabilities.
#'   - `eta`: Initial `eta` values computed from the non-zero counts.
#'   - `eta_bound`: Updated upper bound for `eta`.
#'   - `beta`: Initial regression coefficients.
#'   - `mu`: Initial mean values for the count data.
#'   - `theta`: Initial dispersion parameters.
#'   - `v`: Binary matrix indicating zero counts.
#'   - `tau`: Matrix of initial posterior probabilities for zero inflation.
#'   - `lambda`: Regularization parameter.
#'   - `beta_bound`: Bounds for the `beta` parameter.
#'   - `theta_bound`: Bounds for the `theta` parameter.
#'
#' @importFrom stats median
#'
#' @examples
#' \donttest{
#' # Example usage
#' counts <- matrix(rpois(100, lambda = 10), nrow = 10)
#' covariates <- matrix(runif(20), nrow = 10)
#' params <- initialize_mzinb(counts, covariates)
#' }
#'
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

  list(
    data = data, X = X, n = n, m = m, p = p, pi = pi, eta = eta, eta_bound = eta_bound,
    beta = beta, mu = mu, theta = theta, v = v, tau = tau,
    lambda = lambda, beta_bound = beta_bound, theta_bound = theta_bound
  )
}
