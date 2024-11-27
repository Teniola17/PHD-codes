MAP.Binary2PriorML2C <- function(fit, ...){

  mix.2c <- fitdistr(x = fit$sims.list$p.pred,
                     densfun = function(x, a1, b1, a2, b2, w1){
                       mixBeta(x, ab = list(c(a1, b1), c(a2, b2)), w = w1)
                     },
                     start = list(
                       a1 = 5,
                       b1 = 45,
                       a2 = 1,
                       b2 = 5,
                       w1 = 0.8
                     ),
                     lower = c(0.001, 0.001, 0.001, 0.001, 0.001),
                     upper = c(Inf, Inf, Inf, Inf, 0.999))

  ab <- list(c(mix.2c$estimate[['a1']], mix.2c$estimate[['b1']]),
            c(mix.2c$estimate[['a2']], mix.2c$estimate[['b2']]))
  w <- mix.2c$estimate[['w1']]
  
  prior <- function(x,...){
    y <- mixBeta(x, ab = ab, w = w)
    return(y)
  }
  
  attributes(prior) <- list(params = list(ab = ab, w = w))
  
  return(prior)
}
