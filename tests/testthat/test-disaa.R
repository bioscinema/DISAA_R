# tests/testthat/test-disaa.R

test_that("disaa() runs with counts+covariates (no phyloseq)", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("nloptr")
  skip_if_not_installed("fastDummies")

  set.seed(123)

  # ---- simulate counts (samples x taxa) ----
  N <- 30
  M <- 60
  grp <- rbinom(N, 1, 0.5)
  covariates <- data.frame(grp = grp)

  p <- 2
  eta_true  <- runif(N, -1.0, 1.0)
  beta_true <- cbind(
    rnorm(M, 0, 0.2),                       # intercept
    rbinom(M, 1, 0.3) * rnorm(M, 1.0, 0.3)  # some signal on grp
  )
  MU <- exp(outer(eta_true, rep(1, M)) + as.matrix(cbind(1, grp)) %*% t(beta_true))
  theta_true <- rep(1, M)

  counts <- matrix(rnbinom(N * M, mu = MU, size = theta_true), nrow = N, ncol = M)
  # zero inflation
  pi_true <- 0.4
  counts[matrix(runif(N * M) < pi_true, nrow = N)] <- 0

  # ---- run disaa() WITHOUT phyloseq ----
  out <- disaa(
    phyloseq_obj   = NULL,
    counts         = counts,
    covariates     = covariates,
    formula        = ~ grp,
    conf           = 0.95,
    lambda_beta    = 0.05,
    zero_adj       = FALSE,
    outlier_replace= FALSE,
    beta_bound     = c(-10, 10),
    theta_bound    = c(1e-2, 10),
    eta_bound      = c(-10, 10),
    max_iter       = 3,      # small for CI speed
    tol            = 1e-4,
    seed           = 42
  )

  expect_type(out, "list")
  expect_true(all(c("aux","res") %in% names(out)))

  res <- out$res
  expect_true(all(c("pi","theta","eta","beta","mu","pvalues","significant") %in% names(res)))

  expect_equal(length(res$pi),    M)
  expect_equal(length(res$theta), M)
  expect_equal(length(res$eta),   N)
  expect_equal(dim(res$beta),     c(M, 2))   # intercept + grp
  expect_equal(dim(res$mu),       c(N, M))

  expect_true(is.matrix(res$pvalues))
  expect_equal(nrow(res$pvalues), M)
  expect_equal(nrow(res$significant), M)
})

test_that("disaa() also runs with a phyloseq object (optional)", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("nloptr")
  skip_if_not_installed("fastDummies")
  skip_if_not_installed("phyloseq")

  set.seed(123)

  N <- 20
  M <- 30
  grp <- rbinom(N, 1, 0.5)
  meta <- data.frame(grp = grp)
  rownames(meta) <- paste0("S", seq_len(N))

  eta_true  <- runif(N, -1.0, 1.0)
  beta_true <- cbind(rnorm(M, 0, 0.2), rbinom(M, 1, 0.3) * rnorm(M, 1.0, 0.3))
  MU <- exp(outer(eta_true, rep(1, M)) + as.matrix(cbind(1, grp)) %*% t(beta_true))
  theta_true <- rep(1, M)

  counts <- matrix(rnbinom(N * M, mu = MU, size = theta_true), nrow = N, ncol = M)
  counts[matrix(runif(N * M) < 0.4, nrow = N)] <- 0

  rownames(counts) <- rownames(meta)
  colnames(counts) <- paste0("T", seq_len(M))

  phy <- phyloseq::phyloseq(
    phyloseq::otu_table(counts, taxa_are_rows = FALSE),
    phyloseq::sample_data(meta)
  )

  out <- disaa(
    phyloseq_obj   = phy,
    counts         = NULL,
    covariates     = NULL,
    formula        = ~ grp,
    conf           = 0.95,
    lambda_beta    = 0.05,
    zero_adj       = FALSE,
    outlier_replace= FALSE,
    beta_bound     = c(-10, 10),
    theta_bound    = c(1e-2, 10),
    eta_bound      = c(-10, 10),
    max_iter       = 3,
    tol            = 1e-4,
    seed           = 42
  )

  expect_type(out, "list")
  expect_true(all(c("aux","res") %in% names(out)))
})

