Prior2Posterior <- function(prior, x, n){
  # prior    mixture of beta densities, function with 'params' attribute, which is a list containing 'ab' and 'w'
  # x, n     nb of events x out of n

  
  ## prior parameters
  ab_prior <- attr(prior, 'params')$ab
  w_prior <-  attr(prior, 'params')$w
  
  n_comp <- length(ab_prior)

  if (length(w_prior) < n_comp){
    w_prior <- c(w_prior, 1 - sum(w_prior))
  }
  
  ## posterior parameters
  ab_post <- lapply(ab_prior, function(u) c(u[1] + x, u[2] + n - x) )
  
  marg_lik <- sapply(ab_prior, function(u) dbetabinom.ab(x, n, u[1], u[2]))
  w_post <- marg_lik * w_prior / sum(marg_lik * w_prior)
  
  
  ## construct the posterior given as the output
  post <- function(x, ...){
    y <- mixBeta(x, ab = ab_post, w = w_post)
    return(y)
  }
  
  attributes(post) <- list(params = list(ab = ab_post, w = w_post))
  
  return(post)
}
