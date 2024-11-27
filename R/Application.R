#' ## Section 4: Application (Ulcerative Colitis)
#' PLATFORM:  R 3.0.0 and WinBUGS 1.4.3 (run through wine) run on Mac OS X 10.9  
#' AUTHOR: Sandro Gsteiger, Institute of Social and Preventive Medicine, University of Bern, Bern, Switzerland   

#+, session_init
rm(list=ls())
set.seed(684324)

library(R2WinBUGS)
library(MASS)
library(ggplot2)
library(VGAM)

for (i in dir('src')){
  source(file.path('src', i))
}

sessionInfo()


#' ### Derive standard MAP prior from historical data

#' Read in the historical data
dhist <-  read.csv2('data_colitis.csv', stringsAsFactors=FALSE)
dhist

#' Do the MA and derive the predictive distribution.
wb_data <- list(N_std = nrow(dhist),
                y = dhist$y,
                n = dhist$n,
                prior_prec_tau = 1)

wb_init <- function() {
  list(logit_p = rnorm(wb_data$N_std, mean = -1, sd = 0.5),
       logit_p.pred = rnorm(1, mean = -1, sd = 0.5),
       mu = rnorm(1, mean = -1, sd = 0.5),
       tau = rexp(1, rate = 4)
  )  
}

fit <- bugs(
  data = wb_data,
  inits = wb_init,
  parameters.to.save = c("mu", "tau", "p.pred"),
  model.file = 'BUGS_model.txt',
  n.burnin = 1000, 
  n.iter = 401000,
  n.chains = 5,
  n.thin = 20,
  bugs.seed = 83457,
  # debug  = TRUE,
  # useWINE = TRUE, 
  # bugs.directory = "c:/Program Files/WinBUGS14/"
)

fit

#' Derive various approximations to the predictive distribution:
#' * the Rao-Blackwellized density estimate;
#' * the classical MAP prior: 1 component beta (ML estimation of shape parameters);
#' * 2 component beta mixture;
#' * 3 component beta mixture.
prior_rb <- MAP.Binary2PriorRB(fit) # Rao-Blackwellized density estimate (for logit(p.pred), then backtransformed to p.pred)

prior_1c <- MAP.Binary2PriorML(fit, start= list(shape1 = 1.7, shape2 = 12.3), lower = c(0.001, 0.001))
attributes(prior_1c)$params

prior_2c <- MAP.Binary2PriorML2C(fit)
attributes(prior_2c)$params

prior_3c <- MAP.Binary2PriorML3C(fit)
attributes(prior_3c)$params

#' Plot hist data along with the various approximations to the predictive distribution:  
#' * calculate the density values of the various estimates;  
#' * determine the KL div between RB and the other approximations; 
#' * do the plot.  

eps <- 0.001
dppred <- within(data.frame(p = seq(0, 1, eps)),
                 { y_rb = prior_rb(p)
                 y_1c = prior_1c(p)
                 y_2c = prior_2c(p)
                 y_3c = prior_3c(p)
                 })

dppred.l <- log(dppred[, c("y_rb", "y_1c", "y_2c", "y_3c")])
names(dppred.l) <- paste('l', names(dppred.l), sep='')

KL <- with(data.frame(dppred, dppred.l)[-c(1, nrow(dppred)),], {
  c(map = eps * sum( y_rb * (ly_rb - ly_1c) ),
    mix_2c = eps * sum( y_rb * (ly_rb - ly_2c) ),
    mix_3c = eps * sum( y_rb * (ly_rb - ly_3c) ))
})
round(KL, 3)

offset <- 3 + max(dppred$y_rb, na.rm=TRUE)
x.max <- 0.6
nt <- nrow(dhist)
plot(dppred$p, dppred$y_rb, type = 'l', lty = 1, lwd=4, axes = FALSE,
     xlab = 'Probability of remission', ylab = '',
     xlim = c(0, x.max),
     ylim = c(0, offset + nt + 1)
)
axis(1, pos=c(0, 0))
lines(dppred$p, dppred$y_1c, lty=3, lwd=4)
lines(dppred$p, dppred$y_2c, lty=2, lwd=4, col=grey(0.3))
lines(dppred$p, dppred$y_3c, lty=4, lwd=4, col=grey(0.6))
legend(0.535*x.max, max(dppred$y_rb, na.rm=TRUE), xjust=0, yjust=1,
       legend = c('1 component', '2 components', '3 components', 'Rao-Blackwellized'),
       col = c(1, grey(0.3), grey(0.6), 1),
       lty = c(3, 2, 4, 1),
       lwd = 2,
       bty = 'n')
text(x.max, max(dppred$y_rb, na.rm=TRUE)*1.03,
     paste(c('KL div', round(KL, 3)), collapse='\n'), xpd=TRUE,
     adj = c(1, 1))

dhist2 <- within(dhist, {
  est = y / n
  CI = BinaryExactCI(y, n)
})

for(i in 1:nt) {
  lines(x=dhist2[i, 'CI'], y=rep(offset + nt - i + 1, 2))
  points(x=dhist2[i, 'est'], y=c(offset + nt - i + 1))
  text(x=0.65*x.max, y=c(offset + nt - i + 1), labels=dhist2[i, 'study'], adj=c(0, 0.5))
  text(x=x.max, y=c(offset + nt - i + 1),
       labels=paste(dhist2[i, 'y'], '/', dhist2[i, 'n'], ' (', round(100*dhist2[i, 'est']), '%)', sep=''),
       adj=c(1, 0.5))
}

text(-0.1, offset + 1, 'Historical data', adj=c(0, 1), srt=90, xpd = TRUE)
text(-0.1, 0.3*offset, 'New trial', adj=c(0.5, 1), srt=90, xpd = TRUE)


#' Calculate summaries and the effective sample size (ESS) for the 3-component mixture prior.
round(
  mixBetaSum(attr(prior_3c, 'params')$ab, attr(prior_3c, 'params')$w), 
  2)

ess_3c <- ess(K = 3,
              w = c(attr(prior_3c, 'params')$w, 1 - sum(attr(prior_3c, 'params')$w)),
              alpha = matrix(unlist(attr(prior_3c, 'params')$ab), ncol = 2, byrow = TRUE)[,1],
              beta = matrix(unlist(attr(prior_3c, 'params')$ab), ncol = 2, byrow = TRUE)[,2],
              c = 100)
ess_3c


#' Determine the posterior from the 3C mixture prior and hypothetical outcomes.
dsums <- data.frame(y = c(0, 2, 5, 10, 15), n = 20)
dsums <- cbind(
  dsums, 
  round(t(
    apply(dsums, MAR = 1, FUN = function(u){
      post <- Prior2Posterior(prior_3c, x = u[['y']], n = u[['n']])
      sum <- mixBetaSum(attr(post, 'params')$ab, attr(post, 'params')$w)
      ess <- ess(K = length(attr(post, 'params')$ab),
                 w = attr(post, 'params')$w,
                 alpha = matrix(unlist(attr(post, 'params')$ab), ncol = 2, byrow = TRUE)[,1],
                 beta = matrix(unlist(attr(post, 'params')$ab), ncol = 2, byrow = TRUE)[,2],
                 c = 100)
      return(c(w = attr(post, 'params')$w, sum, ess = ess))
    })
  ), 2)
)
dsums

### End of script
