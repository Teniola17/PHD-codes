MAP.Binary2PriorML <- function(fit, ...){
    mle.fit <- fitdistr(x = fit$sims.list$p.pred,
                        densfun = 'beta',
                        ...)

    prior <- function(x,...){
        y <- dbeta(x, shape1 = mle.fit$estimate[['shape1']], shape2 = mle.fit$estimate[['shape2']])
        return(y)
    }

    attributes(prior) <- list(params = mle.fit$estimate)

    return(prior)
}
