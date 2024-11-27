mixBetaSum <- function(ab, w, seed = NULL, n_sim = 100000){
  # derive summaries for a mixture of beta densities via simulation
  # ab        list of shape parameters of beta components
  # w         vector of mixture weights
  # n_sim     nb of random samples from mixture from which summaries are derived

  if (!is.null(seed)) set.seed(seed)
  
  n_comp <- length(ab)
  ab_mat <- matrix(unlist(ab), ncol = 2, byrow = TRUE)
  
  if (length(w) < n_comp){
    w <- c(w, 1 - sum(w))
  }
  
  ## generate random samples from mixture
  comp <- sample.int(n_comp, size = n_sim, replace = TRUE, prob = w)
  sim <- rbeta(n_sim, ab_mat[comp, 1], ab_mat[comp, 2])

  ## derive summaries
  out <- c(mean = mean(sim),
           CI = quantile(sim, probs = c(0.025, 0.975))
           )

  return(out)
}
