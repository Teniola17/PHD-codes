mixBeta <- function(x, ab, w=NULL){
  # ab list of vectors with elements a, b (shape parameters of each beta component)
    n.comp <- length(ab)
    if (n.comp == 1) {out <- dbeta(x, shape1 = ab[[1]][[1]], shape2 = ab[[1]][[2]])}
    else {
        if (length(w) == (n.comp - 1)) { w <- c(w, 1-sum(w)) }

        out <- 0
        for (i in 1:n.comp){ out <- out + w[i]*dbeta(x, shape1 = ab[[i]][[1]], shape2 = ab[[i]][[2]])}
    }
    return(out)
}
