# Bayesian logistic regression example

options(scipen = 999)

library(tidyverse)
library(rjags)

# https://www4.stat.ncsu.edu/~reich/ABA/code/GLM

con <- url("http://www4.stat.ncsu.edu/~reich/ABA/Code/gambia.RData")
load(con)

dim(X) # 2035 observations of 5 variables
names(X) # age, netuse, treated, green and phc
X <- as.matrix(X)

Y <- pos # Response variable
n <- length(Y)

table(Y)

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

nu <- apply(X, 1, function(x) b %*% x) # Vector of linear predictors
qx <- apply(X, 1, function(x) q(x, b)) # Vector of q(x)

# (proportional to) log posterior
logpost <- function(b) {
  logprior <- sum((b - mu)^2 / 2*sigma)
  nu <- apply(X, 1, function(x) b %*% x)
  loglike <- sum(-log(1 + exp(nu))) + sum(nu[Y == 1])
  logprior + loglike
}

# Function to update jth component of beta via RWMH step
mh <- function(b, j, scale) {
  y <- b
  y[j] <- y[j] + rnorm(1, mean = 0, sd = scale[j])
  a <- exp(logpost(y) - logpost(b)) # Acceptance probability
  if(a > runif(1)) return(y) else return(b)
}

# Metropolis-within-Gibbs sampler (random scan)
wmg <- function(b0, nsim, scale) {
  p <- length(b0)
  r <- array(NA, c(nsim, p))
  r[1, ] <- b0 # init
  for(i in 2:nsim) {
    j <- sample(1:p, 1)
    r[i, ] <- mh(r[i-1, ], j, scale)
  }
  r <- as.data.frame(r)
  names(r) <- sprintf("b%d", 1:p)
  return(r)
}

start_time <- Sys.time()
test <- wmg(b, nsim = 10000, scale = rep(0.1, 6)) # But this is starting at roughly the right values...
end_time <- Sys.time()

end_time - start_time # 1.1 mins

trace_plots <- function(x, nsim) { # Helper function
  a <- qplot(x = 1:nsim, y = x$b1, geom = "line")
  b <- qplot(x = 1:nsim, y = x$b2, geom = "line")
  c <- qplot(x = 1:nsim, y = x$b3, geom = "line")
  d <- qplot(x = 1:nsim, y = x$b4, geom = "line")
  e <- qplot(x = 1:nsim, y = x$b5, geom = "line")
  f <- qplot(x = 1:nsim, y = x$b6, geom = "line")
  cowplot::plot_grid(a, b, c, d, e, f, ncol = 3)
}

trace_plots(test, 10000) # Seems like b2 and b5 are more sticky, maybe reduce scale here? 

scale2 <- c(0.1, 0.0005, 0.1, 0.1, 0.005, 0.1) # Tried a few configurations to see what looks good

test2 <- wmg(b, 10000, scale2)
trace_plots(test2, 10000)

# Okay lets give a larger sample size and starting at zero a go

chains <- wmg(rep(0, 6), 100000, scale2)

trace_plots(chains, 100000) # 1, 2, 5 could still be improved - not sure how - ask Murray, Pier?

colMeans(chains); b # Pretty close to the ML estimates, so it's probably getting the right answer

# Subset modelling --------------------------------------------------------

# Now let's try only using a subset of the predictors in each model

head(X) # Reminder about the variables we have

X1 <- X[, c(1, 2, 5, 6)] # intercept, age, green and phc
p1 <- 4

X2 <- x[, c(1, 3, 4, 6)] # intercept netuse, treated and phc
p2 <- 4