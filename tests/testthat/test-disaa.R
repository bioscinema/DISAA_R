# tests/testthat/test-disaa.R

test_that("disaa function works correctly", {
  set.seed(123)  # Ensure reproducibility

  # Simulate data
  N <- 100  # Number of samples
  M <- 200  # Number of genes

  pi <- 0.5
  theta <- 1
  pi_real <- rep(pi, M)
  theta_real <- rep(theta, M)

  # Generate eta
  eta_real <- runif(N, min = -3, max = 3)

  # Generate beta, 30% significant
  beta_real <- matrix(rbinom(M * 2, 2, 0.3), nrow = M, ncol = 2)

  # Generate X (covariates)
  X <- cbind(rep(1, N), rbinom(N, 1, 0.5))

  # Generate zero-inflated negative binomial counts
  mu <- exp(outer(eta_real, rep(1, M)) + X %*% t(beta_real))
  counts <- matrix(rnbinom(N * M, mu = mu, size = theta_real), nrow = N, ncol = M)
  zero_mask <- matrix(runif(N * M) <= pi, nrow = N, ncol = M)
  counts[zero_mask] <- 0

  # Run the disaa function
  out <- disaa(counts, X, max_iter = 30)

  # Validate the output
  expect_type(out, "list")
  expect_true(all(c("pi", "theta", "eta", "beta", "mu", "pvalues") %in% names(out)))
  expect_equal(length(out$pi), M)
  expect_equal(length(out$theta), M)
  expect_equal(length(out$eta), N)
  expect_equal(dim(out$beta), c(M, ncol(X)))
  expect_equal(dim(out$mu), dim(counts))
})
