MAP.Binary2PriorML3C <- function(fit, ...){
  
  mix.3c <- fitdistr(x = fit$sims.list$p.pred,
                     densfun = function(x, a1, b1, a2, b2, a3, b3, w1, w2){
                       mixBeta(x, ab = list(c(a1, b1), c(a2, b2), c(a3, b3)), w = c(w1, w2))
                     },
                     start = list(
                       a1 = 2.2,
                       b1 = 18,
                       a2 = 14,
                       b2 = 115,
                       a3 = 1,
                       b3 = 3,
                       w1 = 0.5,
                       w2 = 0.42
                     ),
                     lower = c(0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001),
                     upper = c(Inf, Inf, Inf, Inf, Inf, Inf, 0.999, 0.999))
  
  ab <- list(c(mix.3c$estimate[['a1']], mix.3c$estimate[['b1']]),
             c(mix.3c$estimate[['a2']], mix.3c$estimate[['b2']]),
             c(mix.3c$estimate[['a3']], mix.3c$estimate[['b3']]))
  w <- c(mix.3c$estimate[['w1']], mix.3c$estimate[['w2']])
  
  prior <- function(x,...){
    y <- mixBeta(x, ab = ab, w = w)
    return(y)
  }
  
  attributes(prior) <- list(params = list(ab = ab, w = w))
  
  return(prior)
}
