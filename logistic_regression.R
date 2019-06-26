# Bayesian logistic regression example

options(scipen = 999)

library(tidyverse)
library(rjags)

# https://www4.stat.ncsu.edu/~reich/ABA/code/GLM

con <- url("http://www4.stat.ncsu.edu/~reich/ABA/Code/gambia.RData")
load(con)

# Exploration -------------------------------------------------------------

dim(X) # 2035 observations of 5 variables
names(X) # age, netuse, treated, green and phc
X <- as.matrix(X)

Y <- pos # Response variable
n <- length(Y)

table(Y) # 1308 without and 727 with

mle <- glm(Y ~ X, family = "binomial") # Logistic regression using MLE
summary(mle)

b <- mle$coefficients

# Using rjags -------------------------------------------------------------

# Specify the JAGS model

logistic_model <- "model{

  # Likelihood

  for(i in 1:n){
    Y[i] ~ dbern(q[i])
    logit(q[i]) <- beta[1] + beta[2]*X[i, 1] + beta[3]*X[i, 2] +
                   beta[4]*X[i, 3] + beta[5]*X[i, 4] + beta[6]*X[i, 5]
  }

  # Prior

  for(j in 1:6){
    beta[j] ~ dnorm(0, 0.1)
  }

}"

dat <- list(Y = Y, n = n, X = X)

model <- jags.model(textConnection(logistic_model),
                    data = dat, n.chains = 3, quiet = TRUE)

update(model, 10000)

samp <- coda.samples(model, variable.names = c("beta"), 10000)

# From scratch ------------------------------------------------------------

X <- cbind(intercept = 1, X) # Add intercept column to design matrix

# The plan is to use Metropolis-within-Gibbs to sample from the posterior

# Set up the prior on beta (extended to include beta zero) as above
p <- 6
mu <- rep(0, p)
sigma <- rep(0.1, p)

# Classify to 1 with probability
q <- function(x, b) {
  exp(b %*% x) / (1 + exp(b %*% x)) 
}

# (proportional to) log posterior in the indep normals prior case
logpost <- function(b, X, mu, sigma) {
  logprior <- sum((b - mu)^2 / 2*sigma)
  nu <- apply(X, 1, function(x) b %*% x) # Vector of linear predictors
  loglike <- sum(-log(1 + exp(nu))) + sum(nu[Y == 1])
  logprior + loglike
}

# Function to update jth component of beta via RWMH step
mh <- function(b, X, mu, sigma, j, scale) {
  y <- b
  y[j] <- y[j] + rnorm(1, mean = 0, sd = scale[j])
  a <- exp(logpost(y, X, mu, sigma) - logpost(b, X, mu, sigma)) # Acceptance probability
  if(a > runif(1)) return(y) else return(b)
}

# Metropolis-within-Gibbs sampler (random scan)
mwg <- function(b0, X, mu, sigma, scale, nsim) {
  p <- length(b0)
  r <- array(NA, c(nsim, p))
  r[1, ] <- b0 # init
  for(i in 2:nsim) {
    j <- sample(1:p, 1)
    r[i, ] <- mh(r[i-1, ], X, mu, sigma, j, scale)
  }
  r <- as.data.frame(r)
  names(r) <- sprintf("b%d", 1:p)
  return(r)
}

trace_plots_simple <- function(x, nsim) { # Helper function
  a <- qplot(x = 1:nsim, y = x$b1, geom = "line")
  b <- qplot(x = 1:nsim, y = x$b2, geom = "line")
  c <- qplot(x = 1:nsim, y = x$b3, geom = "line")
  d <- qplot(x = 1:nsim, y = x$b4, geom = "line")
  e <- qplot(x = 1:nsim, y = x$b5, geom = "line")
  f <- qplot(x = 1:nsim, y = x$b6, geom = "line")
  cowplot::plot_grid(a, b, c, d, e, f, ncol = 3)
}

# trace_plots <- function(x, p, nsim) { # Helper function
#   for (i in 1:p) {
#     assign(paste0("plot", i), qplot(x = 1:nsim, y = x[[paste0("b", i)]], geom = "line"))
#   }
#   # This isn't working yet...
# }

test <- mwg(b, X, mu, sigma, scale = rep(0.1, 6), nsim = 10^4) # But this is starting at roughly the right values...
trace_plots_simple(test, 10^4) # Seems like b2 and b5 are more sticky, maybe reduce scale here? 

scale2 <- c(0.1, 0.0005, 0.1, 0.1, 0.005, 0.1) # Tried a few configurations here to see what looks best
test2 <- mwg(b, X, mu, sigma, scale = scale2, nsim = 10^4)
trace_plots_simple(test2, 10^4)

# Okay lets give a larger sample size and starting at zero a go

chains <- mwg(rep(0, 6), X, mu, sigma, scale = scale2, nsim = 10^5)
trace_plots_simple(chains, 10^5) # 1, 2, 5 could still be improved - not sure how - ask Murray, Pier?
saveRDS(chains, "results/chains.Rds")

colMeans(chains); b # Pretty close to the ML estimates, so the code is probably correct mostly

# Subset modelling --------------------------------------------------------

head(X) # Reminder about the variables we have
# Now try only using a subset of the predictors in each model

# |   | b1 | b2 | b3 | b4 | b5 | b6 |
# |---|----|----|----|----|----|----|
# | 1 | Y  | Y  | N  | N  | Y  | Y  |
# | 2 | Y  | N  | Y  | Y  | N  | Y  |

# i.e. link parameter is (b1, b6), coefficients corresponding to the intercept term and phc

X1 <- X[, c(1, 2, 5, 6)] # Design matrix 1: intercept, age, green and phc
p1 <- 4

X2 <- X[, c(1, 3, 4, 6)] # Design matrix 2: intercept netuse, treated and phc
p2 <- 4

# Priors same as in the joint modelling to start with (so consistent in the link parameter)
mu1 <- rep(0, p1)
sigma1 <- rep(0.1, p1)

mu2 <- rep(0, p2)
sigma2 <- rep(0.1, p2)

# Run the samplers
chains1 <- mwg(b0 = rep(0, 4), X = X1, mu = mu1, sigma = sigma1,
               scale = scale2[c(1, 2, 5, 6)], nsim = 10^5)
names(chains1) <- c("b1", "b2", "b5", "b6") # Correct the naming
saveRDS(chains1, "results/chains1.Rds")

chains2 <- mwg(b0 = rep(0, 4), X = X2, mu = mu2, sigma = sigma2,
               scale = scale2[c(1, 3, 4, 6)], nsim = 10^5)
names(chains2) <- c("b1", "b3", "b4", "b6") # Correct the naming
saveRDS(chains2, "results/chains2.Rds")