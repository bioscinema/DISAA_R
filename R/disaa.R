#' Differential Abundance Detection via Log-Normal Modeling under Varying Support (DISSA)
#'
#' Performs differential abundance analysis using a zero-inflated negative binomial (ZINB)
#' model with a log-normal latent structure. The method supports covariate adjustment,
#' optional outlier replacement, and multiple testing correction for both binary and
#' multi-class covariates. Input can be provided directly as a count matrix and a covariate
#' data frame, or extracted from a \code{phyloseq} object.
#'
#' @details
#' Let \eqn{Y \in \mathbb{N}^{n \times m}} be the counts (samples \eqn{\times} taxa) and
#' \eqn{X \in \mathbb{R}^{n \times p}} the design matrix (including an intercept).
#' The model assumes, for each taxon \eqn{j} and sample \eqn{i},
#' \deqn{Y_{ij} \sim \text{ZINB}\big(\mu_{ij}, \theta_j, \pi_j\big), \qquad
#' \log \mu_{ij} = \eta_i + x_i^\top \beta_j,}
#' where \eqn{\eta_i} are sample-specific log scaling factors, \eqn{\beta_j} are taxon-specific
#' regression coefficients, \eqn{\theta_j} is the NB overdispersion, and \eqn{\pi_j} is the
#' zero-inflation probability.
#'
#' For inference on a contrast between two levels of a covariate, the function constructs
#' a contrast vector \eqn{C_\beta = (1, -1)} (or the appropriate multi-class contrast) and uses
#' the observed information to form a Wald statistic. For taxon \eqn{j}, the test statistic is
#' \deqn{ \mathrm{stat}_j \;=\; \frac{(\beta_{j,\text{lvl1}} - \beta_{j,\text{lvl2}})^2}
#'              {\,C \; \mathrm{cov} \; C^\top\,}, }
#' which is compared to an \eqn{F(1, n - p)} reference under the null. Raw p-values are then
#' adjusted for multiple comparisons according to the requested level.
#'
#' Optionally, a sample-wise zero proportion can be included as an additional covariate
#' (\code{zero_adj = TRUE}). Outlier replacement (if enabled) uses a NB GLM fitted on
#' a subset of samples and flags outliers via Cook's distance at the specified confidence.
#'
#' @param phyloseq_obj A \code{phyloseq} object containing an OTU table and sample metadata.
#'   If provided, \code{counts} and \code{covariates} are extracted from it (taxa are assumed
#'   to be columns/samples are rows; rows are transposed if needed).
#' @param counts Optional numeric matrix of raw counts with dimensions
#'   \emph{samples} \eqn{\times} \emph{taxa}. Required when \code{phyloseq_obj} is \code{NULL}.
#' @param covariates Optional data frame of sample-level covariates referenced in \code{formula}.
#'   Required when \code{phyloseq_obj} is \code{NULL}.
#' @param formula Right-hand-sided formula specifying covariates of interest (e.g., \code{~ Group + Age}).
#' @param conf Confidence level used to determine significance (default \code{0.95}).
#' @param lambda_beta Nonnegative regularization strength applied to \eqn{\beta} blocks
#'   in the information matrix (default \code{0.5}).
#' @param zero_adj Logical; if \code{TRUE}, augments the design with a sample-wise zero proportion
#'   covariate (default \code{FALSE}).
#' @param outlier_replace Logical; if \code{TRUE}, enables outlier replacement using NB-GLM
#'   diagnostics (default \code{FALSE}).
#' @param outlier_conf Confidence level for outlier detection via Cook's distance (default \code{0.99}).
#' @param beta_bound Numeric length-2 vector giving lower/upper bounds for regression coefficients
#'   (default \code{c(-10, 10)}).
#' @param theta_bound Numeric length-2 vector giving bounds for overdispersion \eqn{\theta_j}
#'   (default \code{c(1e-2, 10)}).
#' @param eta_bound Numeric length-2 vector giving bounds for \eqn{\eta_i}
#'   (default \code{c(-10, 10)}).
#' @param max_iter Integer; maximum EM iterations.
#' @param tol Convergence tolerance for the EM inner optimizers and overall loop (default \code{1e-6}).
#' @param seed Integer random seed for initialization and any stochastic steps (default \code{42}).
#'
#' @return An object of class \code{"DISSA"}: a list with components
#' \describe{
#'   \item{res}{Main results: \code{pi}, \code{theta}, \code{eta}, \code{beta}, \code{mu},
#'     \code{pvalues}, and \code{significant}.}
#'   \item{aux}{Auxiliary elements used internally: data matrices, design matrix,
#'     observed information blocks, and factor expansion bookkeeping.}
#' }
#'
#' @section Output shapes:
#' \itemize{
#'   \item \code{pi}, \code{theta}: length-\eqn{m} (taxa).
#'   \item \code{eta}: length-\eqn{n} (samples).
#'   \item \code{beta}: \eqn{m \times p}.
#'   \item \code{mu}: \eqn{n \times m}.
#'   \item \code{pvalues}: \eqn{m \times q} matrix (taxa \eqn{\times} contrasts).
#' }
#'
#' @note The function expects \code{counts} to be samples-by-taxa. If a \code{phyloseq} object
#'   has taxa as rows, they are transposed internally.
#'
#'
#' @importFrom methods new
#' @importFrom Matrix bdiag
#' @importFrom MASS ginv glm.nb
#' @importFrom nloptr nloptr
#' @importFrom stats pf qf median cooks.distance
#' @importFrom phyloseq taxa_are_rows
#' @importFrom fastDummies dummy_cols
#'
#' @export
disaa = function(phyloseq_obj, counts, covariates,formula, conf = 0.95,
                 lambda_beta = 0.05, zero_adj = FALSE, outlier_replace = FALSE, outlier_conf = 0.99,
                 beta_bound = c(-10, 10), theta_bound = c(1e-2, 10), eta_bound = c(-10, 10), max_iter=50, tol = 1e-6,seed=42){
  set.seed(seed)
  ### Check all the inputs are correct
  check_DISSA_inputs(phyloseq_obj, counts, covariates, formula,
                     conf, lambda_beta, zero_adj, outlier_replace, outlier_conf,
                     beta_bound, theta_bound, eta_bound, max_iter, tol)

  #################################################################################################################################
  # Initialization
  #################################################################################################################################

  ### Extract counts and covariates from phyloseq_obj
  if(!is.null(phyloseq_obj)){
    counts = phyloseq_obj@otu_table@.Data
    if(taxa_are_rows(phyloseq_obj)){
      counts = t(counts)
    }
    covariates = phyloseq_obj@sam_data
  }

  ### Get variables in the formula
  X = covariates[, all.vars(formula), drop = FALSE]

  ### Identify columns that are neither numeric nor factor
  invalid_cols = names(X)[!(sapply(X, is.numeric) | sapply(X, is.factor))]
  if(length(invalid_cols) > 0){
    stop(paste("Error: The following columns are neither numeric nor factor:", paste(invalid_cols, collapse = ", ")))
  }

  ### Deal with factor columns
  factors = names(X)[sapply(X, function(x) is.factor(x))]
  multi_class_factors = names(X)[sapply(X, function(x) is.factor(x))]
  multi_class_dummy_list = list()
  if(length(factors) > 0){
    X = dummy_cols(X, select_columns = factors, remove_first_dummy = FALSE, remove_selected_columns = TRUE)
    for(factor_col in multi_class_factors){
      multi_class_dummy_list[[factor_col]] = c(FALSE, grepl(paste0("^", factor_col, "_"), colnames(X)))
      if(zero_adj){
        multi_class_dummy_list[[factor_col]] = c(FALSE, multi_class_dummy_list[[factor_col]])
      }
    }
  }

  ### If no factors, add intercept
  #if(length(factors) == 0){
  Intercept = rep(1, nrow(X))
  X = cbind(Intercept, X)

  ### Save rownames and colnames, change to matrix
  data = as.matrix(counts)
  data_colnames = colnames(data)
  row.names(data) = NULL; colnames(data) = NULL;

  ### Adjust for zero proportion
  if(zero_adj){
    pct0 = function(x){1 - sum(x == 0) / length(x)}
    zero_adj_ = apply(data, 1, pct0)
    X = cbind(zero_adj_, X)
  }


  X_rownames = row.names(X); X_colnames = colnames(X);
  row.names(X) = NULL; colnames(X) = NULL;
  X = as.matrix(X)

  ### Initialize all parameters
  n = dim(data)[1]
  m = dim(data)[2]
  p = dim(X)[2]

  pi = rep(0.5, m)
  eta = apply(data, 1, function(x) median(log(x[x != 0])))
  eta = ifelse(is.na(eta), 0, eta)
  eta = ifelse(eta > eta_bound[2], eta_bound[2] - tol, eta)
  eta = ifelse(eta < eta_bound[1], eta_bound[1] + tol, eta)

  beta = matrix(1, m, p)
  mu = matrix(0, n, m)
  theta = rep(1, m)
  v = 1 * (data == 0)
  tau = v

  #################################################################################################################################
  # Estimation
  #################################################################################################################################

  ### Estimation
  est = run_em_optimization(data = data, X = X, eta = eta, v = v, tau = tau, beta = beta, theta = theta, pi = pi, mu = mu,
                            lambda_beta = lambda_beta, beta_bound = beta_bound, theta_bound = theta_bound, eta_bound = eta_bound,
                            max_iter = max_iter, tol = tol)
  beta = est$beta; theta = est$theta; pi = est$pi; eta = est$eta; tau = est$tau; mu = est$mu;
  mu_theta = est$mu_theta; sigma2_theta = est$sigma2_theta;


  ### Compute covariance matrix and pvalues
  res = compute_pvalues(data, X, mu, tau, v, theta, beta, lambda_beta, zero_adj, formula,
                        mu_theta, sigma2_theta, multi_class_factors, multi_class_dummy_list, X_colnames)
  I = res$I; cov = res$cov; pvalues_all = res$pvalues_all;

  ### Replace outliers and re-estimate
  if(outlier_replace){

    ### Get the non-significant taxas
    nonsig = which(apply(pvalues_all < 1 - conf, 1, max) == 0)

    counts_copy = counts

    ### Replace outliers for non-significant taxa
    for (j in nonsig){

      ### Adjusted counts by removing eta effect
      adj_y_j = round(data[, j] - exp(eta))

      ### Identify positive counts and likely non-inflated samples
      pos_j = which(adj_y_j > 0)
      order_j = order(tau[, j], decreasing = TRUE)
      non_pos_j = setdiff(order_j, pos_j)

      ### Estimate number of non-inflated samples
      noninflated_j = round(n * (1 - pi[j]) - length(pos_j))

      ### Select sample indices
      if(length(non_pos_j) > 0 && noninflated_j <= length(non_pos_j)){
        index_j = union(pos_j, non_pos_j[seq_len(noninflated_j)])
      }else{
        index_j = order_j
      }

      ### Subsample adjusted counts
      adj_y_sub = adj_y_j[index_j]
      n_j = length(adj_y_sub)
      if(mean(adj_y_sub > 0) < 0.2 || n_j < 8) next

      ### Construct model input
      newdata_j = data.frame(
        y = ifelse(adj_y_sub > 0, adj_y_sub, 0),
        X = X[index_j, 2:p, drop = FALSE]
      )

      ### Fit NB model and detect outliers
      tryCatch({
        suppressWarnings({
          nb_model_j = glm.nb(y ~ ., data = newdata_j)
          cutoff = qf(outlier_conf, df1 = p - 1, df2 = n_j - p + 1)
          outliers_j = which(cooks.distance(nb_model_j) > cutoff)

          if(length(outliers_j) > 0){
            outliers_idx = index_j[outliers_j]
            keep_idx = setdiff(index_j, outliers_idx)

            if(length(keep_idx) > 0){
              ### Replace outliers with mean + eta
              counts_copy[outliers_idx, j] = round(mean(counts_copy[keep_idx, j]) + exp(eta[keep_idx]))
            }
          }
        })
      }, error = function(e) {
        ### Model fitting failed — skip
      })
    }

    ### Re-estimation
    data = as.matrix(counts_copy)
    est = run_em_optimization(data = data, X = X, eta = eta, v = v, tau = tau, beta = beta, theta = theta, pi = pi, mu = mu,
                              lambda_beta = lambda_beta, beta_bound = beta_bound, theta_bound = theta_bound, eta_bound = eta_bound,
                              update_index = nonsig, max_iter = max_iter, tol = tol)
    beta = est$beta; theta = est$theta; pi = est$pi; eta = est$eta; tau = est$tau; mu = est$mu;
    mu_theta = est$mu_theta; sigma2_theta = est$sigma2_theta;

    ### Re-compute covariance matrix and pvalues
    res = compute_pvalues(data, X, mu, tau, v, theta, beta, lambda_beta, zero_adj, formula,
                          mu_theta, sigma2_theta, multi_class_factors, multi_class_dummy_list, X_colnames)
    I = res$I; cov = res$cov; pvalues_all = res$pvalues_all;

  }


  #################################################################################################################################
  # Finalize results
  #################################################################################################################################

  ### Assign names
  names(pi) = data_colnames; names(theta) = data_colnames; names(eta) = X_rownames;
  colnames(beta) = X_colnames; row.names(beta) = data_colnames;
  colnames(mu) = data_colnames; row.names(mu) = X_rownames;
  row.names(pvalues_all) = data_colnames;

  ### Get significant values
  significant = pvalues_all < (1-conf)

  ### Construct return list
  aux = list(data = data, X = X, beta = beta, I = I, cov = cov, multi_class_factors = multi_class_factors,
             multi_class_dummy_list = multi_class_dummy_list, X_colnames = X_colnames)
  res = list(pi = pi, theta = theta, eta = eta, beta = beta, mu = mu, pvalues = pvalues_all, significant = significant)
  out = list(aux = aux, res = res)
  class(out) = "DISSA"
  return(out)


}


