MAP.Binary2PriorRB <- function(fit, ...){

    prior <- function(x, ...){
        ## vectorized version with apply(...) leads to memory error (would need usage of memory.size)
        M <- length(fit$sims.list$mu)
        logit.x <- logit(x)
        y <- x*0
        for (i in 1:M) {
            y <- y + (1/M)*dnorm(logit.x, mean = fit$sims.list$mu[i], sd = fit$sims.list$tau[i])
        }
        y <- y*(1/x + 1/(1-x))
        return(y)
    }

    return(prior)
}
