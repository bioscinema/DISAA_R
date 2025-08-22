#' Compute Contrast-Specific P-Values for DISSA Model
#'
#' Evaluates hypothesis tests for a specified pairwise contrast of a categorical
#' predictor within the DISSA model, returning adjusted p-values for each
#' taxon based on the Wald statistic and F-distribution.
#'
#' @param DISSA_obj A list returned by \code{run_em_optimization()} containing
#'   fitted model components and the observed information matrix in \code{aux}.
#' @param contrast Character vector of length three: \code{c("varname", "value1", "value2")},
#'   specifying the factor variable and the two levels to contrast (level "value1" vs "value2").
#' @param conf Numeric scalar in (0,1), confidence level for p-value adjustment (default = 0.95).
#' @param alpha Numeric scalar in (0,1), target false discovery rate for multiple testing (default = 0.1).
#'
#' @return A numeric vector of length equal to the number of taxa, containing
#'   FDR-adjusted p-values for the specified contrast across all taxa.
#'
#' @details
#' Internally, the function constructs a contrast vector \code{C_beta} = (1, -1) for the
#' two specified levels and embeds it within the full parameter covariance matrix
#' obtained from \code{DISSA_obj$aux$cov}.
#'
#' For each taxon, the test statistic is computed as
#' \deqn{stat_j = \frac{(\beta_{j,value1} - \beta_{j,value2})^2}{C \; cov \; C^T}}{
#'   stat_j = (beta_j,value1 - beta_j,value2)^2 / (C cov C^T),
#' }
#' which follows an F(1, n - p) distribution under the null hypothesis.
#' Raw p-values are then adjusted for multiple comparisons using the chosen FDR level.
#'
#' @noRd
compute_pvalues_contrast = function(DISSA_obj, contrast = c("varname", "value1", "value2"), conf = 0.95, alpha = 0.1){

  ### Extract values
  aux = DISSA_obj$aux
  data = aux$data; X = aux$X; X_colnames = aux$X_colnames; beta = aux$beta; cov = aux$cov;
  n = nrow(data)
  m = ncol(data)
  p = ncol(X)

  ### Get contrast value location
  varname = contrast[1]; value1 = contrast[2]; value2 = contrast[3];
  name1 = paste0(varname, '_', value1); name2 = paste0(varname, '_', value2);
  index1 = which(X_colnames == name1); index2 = which(X_colnames == name2);

  ### Compute pvalues
  C_beta = c(1, -1)
  C_cov = matrix(0, nrow = 1, ncol = p)
  C_cov[1, index1] = -1
  C_cov[1, index2] = 1
  pvals = vapply(1:m, function(j) {

    beta_tmp = sum(beta[j, c(index1, index2)] * C_beta)
    C = cbind(matrix(0, nrow = 1, ncol = m + n + p * (j - 1)), C_cov, matrix(0, nrow = 1, ncol = p * (m - j)))
    stat = beta_tmp^2 / (C %*% cov %*% t(C))
    1 - pf(stat, df1 = 1, df2 = n - p)

  }, numeric(1))
  pvalues_adjust(pvals, conf, alpha)

}

#' Adjust P-Values
#'
#' Performs Benjamini–Hochberg FDR adjustment on a vector of p-values, then
#' iteratively removes the largest raw p-values until the proportion of
#' adjusted p-values below the significance threshold meets the target FDR.
#'
#' @param pvals Numeric vector of raw p-values.
#' @param conf Numeric scalar in (0,1), confidence level for interval cutoff
#'   (default 0.95).
#' @param alpha Numeric scalar in (0,1), desired false discovery rate (default 0.1).
#'
#' @return A numeric vector of the same length as \code{pvals}, containing the
#'   adjusted p-values after iterative selection.
#'
#' @details
#' 1. Apply \code{p.adjust(..., method = "BH")} to the current set of p-values.
#' 2. Compute the proportion of adjusted p-values less than \code{1 - conf}.
#' 3. If that proportion is at least \code{alpha}, stop; otherwise, remove the
#'    largest raw p-value from consideration and repeat.
#'
#' @noRd
pvalues_adjust = function(pvals, conf = 0.95, alpha = 0.1) {

  padj = pvals
  keep = rep(TRUE, length(pvals))
  index = 1:length(pvals)

  repeat {

    padj[keep] = stats::p.adjust(pvals[keep], method = "BH")
    prop = mean(padj < (1 - conf))
    if (prop >= alpha | sum(keep) == 0) break

    ### Remove the largest original p-value
    idx_to_remove = index[keep][which.max(pvals[keep])]
    keep[idx_to_remove] = FALSE
  }

  return(padj)

}
