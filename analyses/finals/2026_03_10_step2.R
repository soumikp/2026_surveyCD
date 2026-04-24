################################################################################
# Weighted PC-Stable Algorithm for Ordinal Data
#
# Modifications:
# 1. Swapped "Normal Scores" for true "Weighted Polychoric Correlations" 
#    using the wCorr package to properly handle discrete ordinal thresholds.
# 2. Added Matrix::nearPD to ensure the correlation matrix is positive-definite.
################################################################################

# Required packages
# install.packages(c("wCorr", "pcalg", "Matrix"))
pacman::p_load(wCorr, pcalg, Matrix)

# ------------------------------------------------------------------------------
# 1. Survey-Weighted Polychoric Correlation Matrix
# ------------------------------------------------------------------------------
weighted_copula_cormat <- function(data, weights) {
  n <- nrow(data)
  p <- ncol(data)
  varnames <- colnames(data)
  if (is.null(varnames)) varnames <- paste0("V", seq_len(p))
  
  cormat <- matrix(1, p, p)
  colnames(cormat) <- varnames
  rownames(cormat) <- varnames
  
  cat("Computing survey-weighted polychoric correlation matrix...\n")
  
  # Compute pairwise weighted polychoric correlations
  for (i in seq_len(p - 1)) {
    for (j in (i + 1):p) {
      # The wCorr package handles the survey weights and the ordinal nature natively
      rho <- wCorr::weightedCorr(x = as.numeric(data[, i]), 
                                y = as.numeric(data[, j]), 
                                method = "polychoric", 
                                weight = weights)
      cormat[i, j] <- rho
      cormat[j, i] <- rho
    }
  }
  
  # Pairwise polychoric matrices are not mathematically guaranteed to be 
  # strictly positive-definite. pcalg requires a PD matrix.
  # We project it to the nearest Positive-Definite matrix just in case.
  if (!matrixcalc::is.positive.definite(cormat)) {
    cormat <- as.matrix(Matrix::nearPD(cormat, corr = TRUE)$mat)
  }
  
  return(cormat)
}

# ------------------------------------------------------------------------------
# 2. Weighted PC-Stable Algorithm
# ------------------------------------------------------------------------------
pc_stable_weighted <- function(data, weights, alpha = 0.01, verbose = FALSE) {
  
  # 1. Get the proper latent correlation matrix
  cormat <- weighted_copula_cormat(data, weights)
  
  # 2. Define the sufficient statistics for the Gaussian conditional independence test
  # Note: We use the effective sample size (n_eff) or raw N depending on your 
  # specific survey design effect. Using raw N is standard for pseudo-likelihoods.
  suffStat <- list(C = cormat, n = nrow(data))
  
  # 3. Run the standard PC algorithm using the Gauss CI test on the latent copula
  pc_fit <- pcalg::pc(suffStat = suffStat,
                      indepTest = gaussCItest,
                      alpha = alpha,
                      labels = colnames(data),
                      verbose = verbose)
  
  # 4. Extract the adjacency matrix (Skeleton / CPDAG)
  # pcalg uses: amat[i,j] = 1 means i -> j. (Opposite of OrdCD's internal logic)
  amat <- as(pc_fit@graph, "matrix")
  
  return(list(
    pc_fit = pc_fit,
    amat = amat,
    cormat = cormat
  ))
}

# Example extraction for swa-oBN step 2:
result <- pc_stable_weighted(data.frame(shep_clustered[,str_detect(colnames(shep_clustered), "o_|n_")]), 
                             shep_clustered |> pull("WEIGHT"))

saveRDS(result, file.path(here(), "analyses", "finals", "2026_03_10_step2_output.Rds"))
