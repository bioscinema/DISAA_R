#' Differential Abundance Detection via Log-Normal Modeling under Varying Support (DISSA)
#'
#' Performs differential abundance analysis using a zero-inflated negative binomial model with log-normal structure.
#' This function supports covariate adjustment, outlier replacement, and multiple testing correction for both binary and multi-class covariates.
#'
#' @param phyloseq_obj A `phyloseq` object containing OTU table and sample metadata. If provided, `counts` and `covariates` will be extracted from this object.
#' @param counts A numeric matrix of raw count data (samples x taxa). Only required if `phyloseq_obj` is not provided.
#' @param covariates A data frame of sample-level covariates. Only required if `phyloseq_obj` is not provided.
#' @param formula A right-hand-sided formula specifying the covariates of interest (e.g., `~ Group + Age`).
#' @param conf Confidence level for identifying significant taxa (default = 0.95).
#' @param lambda Regularization strength for penalizing the inverse squared overdispersion (default = 0.5).
#' @param zero_adj Logical; if `TRUE`, adjusts for the sample-wise zero proportion in the model (default = `FALSE`).
#' @param outlier_replace Logical; if `TRUE`, performs outlier replacement based on deviance residuals and Cook’s distance (default = `FALSE`).
#' @param outlier_conf Confidence level used in outlier detection (default = 0.99).
#' @param beta_bound Bounds for coefficients (default = `c(-10, 10)`).
#' @param theta_bound Bounds for overdispersion parameters (default = `c(1e-2, 10)`).
#' @param eta_bound Bounds for log-scaling factors (default = `c(-10, 10)`).
#' @param max_iter Maximum number of iterations for the EM algorithm (default = 50).
#' @param tol Convergence tolerance for EM optimization (default = 1e-6).
#'
#' @return A list containing:
#' \describe{
#'   \item{pi}{Estimated dropout probabilities for each taxon.}
#'   \item{theta}{Estimated overdispersion parameters.}
#'   \item{eta}{Estimated sample-specific log-scaling factors.}
#'   \item{beta}{Coefficient matrix for covariate effects.}
#'   \item{mu}{Estimated expected counts under the model.}
#'   \item{pvalues}{Adjusted p-values for all covariates.}
#'   \item{significant}{Logical matrix indicating significance under FDR or Bonferroni correction.}
#'   \item{multiple_pvalues}{(Optional) Multi-group global and pairwise p-values for multi-class covariates.}
#'   \item{multiple_significant}{(Optional) Significance matrix for multi-group tests.}
#' }
#'
#'
#' @export
DISSA = function(phyloseq_obj, counts, covariates, formula, conf = 0.95, 
                 lambda = 0.5, zero_adj = FALSE, outlier_replace = FALSE, outlier_conf = 0.99, 
                 beta_bound = c(-10, 10), theta_bound = c(1e-2, 10), eta_bound = c(-10, 10), max_iter = 50, tol = 1e-6){
  
  ### Check all the inputs are correct
  check_DISSA_inputs(phyloseq_obj, counts, covariates, formula,
                     conf, lambda, zero_adj, outlier_replace, outlier_conf,
                     beta_bound, theta_bound, eta_bound, max_iter, tol)
  
  #################################################################################################################################
  # Initialization
  #################################################################################################################################
  
  ### Extract counts and covariates from phyloseq_obj
  if(!is.null(phyloseq_obj)){
    counts = phyloseq_obj@otu_table@.Data
    if(taxa_are_rows(ps)){
      counts = t(counts)
    }
    covariates = phyloseq_obj@sam_data
  }
  
  ### Get variables in the formula, add intercept
  Intercept = rep(1, nrow(covariates))
  X = cbind(Intercept, covariates[, all.vars(formula), drop = FALSE])
  
  ### Identify columns that are neither numeric nor factor
  invalid_cols = names(X)[!(sapply(X, is.numeric) | sapply(X, is.factor))]
  if(length(invalid_cols) > 0){
    stop(paste("Error: The following columns are neither numeric nor factor:", paste(invalid_cols, collapse = ", ")))
  }
  
  ### Deal with factor columns
  factors = names(X)[sapply(X, function(x) is.factor(x))]
  multi_class_factors = names(X)[sapply(X, function(x) is.factor(x) && length(unique(x)) > 2)]
  multi_class_dummy_list = list()
  reference_groups = list()
  if(length(factors) > 0){
    for(factor_col in multi_class_factors){
      reference_groups[[factor_col]] = levels(X[[factor_col]])[1]
    }
    X = dummy_cols(X, select_columns = factors, remove_first_dummy = TRUE, remove_selected_columns = TRUE)
    for(factor_col in multi_class_factors){
      multi_class_dummy_list[[factor_col]] = grepl(paste0("^", factor_col, "_"), colnames(X))
    } 
  }
  
  ### Save rownames and colnames, change to matrix
  data = as.matrix(counts)
  data_colnames = colnames(data)
  row.names(data) = NULL; colnames(data) = NULL;
  X_rownames = row.names(X); X_colnames = colnames(X);
  
  ### Adjust for zero proportion
  if(zero_adj){
    pct0 = function(x){1 - sum(x == 0) / length(x)}
    X = cbind(X[, 1], apply(data, 1, pct0), X[, 2:dim(X)[2]])
  }
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
                            lambda = lambda, beta_bound = beta_bound, theta_bound = theta_bound, eta_bound = eta_bound, 
                            max_iter = max_iter, tol = tol)
  beta = est$beta; theta = est$theta; pi = est$pi; eta = est$eta; tau = est$tau; mu = est$mu;
  
  ### Compute covariance matrix and pvalues
  res = compute_single_var_pvalues(data, X, mu, tau, v, theta, beta, lambda, zero_adj)
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
                              lambda = lambda, beta_bound = beta_bound, theta_bound = theta_bound, eta_bound = eta_bound, 
                              update_index = nonsig, max_iter = max_iter, tol = tol)
    beta = est$beta; theta = est$theta; pi = est$pi; eta = est$eta; tau = est$tau; mu = est$mu;
    
    ### Re-compute covariance matrix and pvalues
    res = compute_single_var_pvalues(data, X, mu, tau, v, theta, beta, lambda, zero_adj)
    I = res$I; cov = res$cov; pvalues_all = res$pvalues_all;
    
  }
  
  #################################################################################################################################
  # Multigroup test
  #################################################################################################################################
  
  if(length(multi_class_dummy_list) > 0){ 
    
    multiple_pvalues = NULL
    multiple_significant = NULL
    multiple_correction = c()
    multiple_colnames = c()
    
    ### Multiple test
    for(i in 1:length(multi_class_dummy_list)){
      
      col_num = multi_class_dummy_list[[i]]
      col_num = if (zero_adj) which(c(FALSE, col_num)) else which(col_num)
      
      factor_name = names(multi_class_dummy_list)[i]
      reference = reference_groups[[i]]
      others = sub("^.*_", "", X_colnames[multi_class_dummy_list[[i]]])
      len = length(others)
      pairs = len * (len + 1) / 2
      
      ### Global test
      pvals_global = vapply(1:m, function(j) {
        beta_tmp = beta[j, col_num]
        inv_cov_tmp = I[p * (j - 1) + col_num, p * (j - 1) + col_num]
        stat = t(beta_tmp) %*% inv_cov_tmp %*% beta_tmp
        1 - pf(stat, df1 = length(col_num), df2 = n - p)
      }, numeric(1))
      
      multiple_pvalues = cbind(multiple_pvalues, stats::p.adjust(pvals_global, "BH"))
      multiple_correction = c(multiple_correction, 1)
      multiple_colnames = c(multiple_colnames, factor_name)
      
      ### Paired experimental groups vs reference
      multiple_pvalues = cbind(multiple_pvalues, pvalues_all[, col_num - 1 - as.numeric(zero_adj)])
      multiple_colnames = c(multiple_colnames, paste0(factor_name, ": ", reference, " vs ", others))
      
      ### Paired among experimental groups
      for(v1 in 1:(len-1)){
        for(v2 in (v1+1):len){
          R = matrix(0, nrow = 1, ncol = len)
          R[1, v1] = 1; R[1, v2] = -1
          pvals_pair = vapply(1:m, function(j) {
            beta_tmp = beta[j, col_num]
            cov_tmp = cov[p * (j - 1) + col_num, p * (j - 1) + col_num]
            stat = t(R %*% beta_tmp) %*% ginv(R %*% cov_tmp %*% t(R)) %*% (R %*% beta_tmp)
            1 - pf(stat, df1 = 1, df2 = n - p)
          }, numeric(1))
          multiple_pvalues = cbind(multiple_pvalues, stats::p.adjust(pvals_pair, "BH"))
          multiple_colnames = c(multiple_colnames, paste0(factor_name, ": ", others[v1], " vs ", others[v2]))
        }
      }
      
      multiple_correction = c(multiple_correction, rep(1 / pairs, pairs))
    }
    
    ### Assign names
    colnames(multiple_pvalues) = multiple_colnames; row.names(multiple_pvalues) = data_colnames;
    ### Bonferroni correction for multi-class 
    alpha_multiple = multiple_pvalues
    alpha_multiple[, ] = 1 - conf
    for(i in 1:length(multiple_correction)){
      alpha_multiple[, i] = alpha_multiple[, i] * multiple_correction[i]
    }
    ### Get significant values
    multiple_significant = multiple_pvalues < alpha_multiple
  }
  
  #################################################################################################################################
  # Finalize results
  #################################################################################################################################
  
  ### Remove beta for added coveriate
  if(zero_adj){
    beta = cbind(beta[, 1], beta[, 3:p])
  }
  
  ### Assign names
  names(pi) = data_colnames; names(theta) = data_colnames; names(eta) = X_rownames; 
  colnames(beta) = X_colnames; row.names(beta) = data_colnames;
  colnames(mu) = data_colnames; row.names(mu) = X_rownames;
  colnames(pvalues_all) = X_colnames[2:length(X_colnames)]; row.names(pvalues_all) = data_colnames;
  
  ### Significant 
  alpha_all = pvalues_all
  alpha_all[, ] = 1 - conf
  ### Bonferroni correction for multi-class 
  if(length(multi_class_dummy_list) > 0){
    for(i in 1:length(multi_class_dummy_list)){
      col = X_colnames[multi_class_dummy_list[[i]]]
      alpha_all[, col] = alpha_all[, col] / length(col)
    } 
  }
  ### Get significant values
  significant = pvalues_all < alpha_all
  
  ### Return
  if(length(multi_class_dummy_list) > 0){
    return(list(pi = pi, theta = theta, eta = eta, beta = beta, mu = mu, pvalues = pvalues_all, significant = significant,
                multiple_pvalues = multiple_pvalues, multiple_significant = multiple_significant))
  }else{
    return(list(pi = pi, theta = theta, eta = eta, beta = beta, mu = mu, pvalues = pvalues_all, significant = significant))    
  }
  
}


