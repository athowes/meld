# Bayesian logistic regression example

# Set-up ------------------------------------------------------------------

options(scipen = 999)
pacman::p_load(tidyverse, rjags, geoR, bayesplot)

data(gambia)
X <- as.matrix(gambia[, c(4:8)])
X <- cbind(intercept = 1, X) # Add intercept column to design matrix
Y <- gambia[, 3] # Response variable
n <- length(Y)

mle <- glm(Y ~ X - 1, family = "binomial") # Logistic regression using MLE
summary(mle)

b <- mle$coefficients

# Full Model --------------------------------------------------------------

# The plan is to use Metropolis-within-Gibbs to sample from the posterior

# Set up the prior on beta (extended to include beta zero)
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
  loglike <- sum(nu[Y == 1]) + sum(-log(1 + exp(nu)))
  logprior + loglike
}

# Function to update jth component of beta via RWMH step
mh <- function(b, X, mu, sigma, j, scale) {
  y <- b
  y[j] <- y[j] + rnorm(1, mean = 0, sd = scale)
  a <- exp(logpost(y, X, mu, sigma) - logpost(b, X, mu, sigma)) # Acceptance probability
  if(a > runif(1)) return(y) else return(b)
}

# Metropolis-within-Gibbs sampler (random scan by default)
mwg <- function(b0, X, mu, sigma, scale, nsim, random_scan = TRUE) {
  p <- length(b0)
  r <- array(NA, c(nsim, p)) # For the chain
  r[1, ] <- b0 # Init chain
  s <- array(0, c(3, p)) # For acceptance rates
  for(i in 2:nsim) {
    if(random_scan) {
      j <- sample(1:p, 1) # Random scan
    } else {
      j <- (i %% p) + 1 # Systematic scan
    }
    s[1, j] <- s[1, j] + 1 # Update pick count
    r[i, ] <- mh(r[i-1, ], X, mu, sigma, j, scale[j])
    if(!identical(r[i, ], r[i-1, ])) s[2, j] <- s[2, j] + 1 # Update accept count
  }
  r <- as.data.frame(r)
  names(r) <- sprintf("b%d", 0:(p-1))
  s[3, ] <- s[2, ] / s[1, ] # Acceptance rates
  return(list("chain" = r, "accept" = s))
}

# Function to calculate the total acceptance rate retroactively (this may not be so useful)
accept_rate <- function(b, chains) {
  lag <- rbind(as.matrix(t(b)), as.matrix(chains[-nrow(chains),]))
  sum(rowSums(chains - lag) == 0) / nrow(chains)
}

# Scaling and Running -----------------------------------------------------

test <- mwg(b, X, mu, sigma, scale = rep(0.1, 6), nsim = 10^4) # Initialising at roughly the right values
mcmc_trace(test$chain)
test$accept # Plan is to adjust scale based on optimal acceptance rate 0.234

# Optimal scalings for each parameter should roughly be some constant multiplied by their standard error
se <- sqrt(diag((vcov(mle))))

test_se <- mwg(b, X, mu, sigma, scale = se, nsim = 10^4)
test_se$chain
test_se$accept

# Objective function to be optimised
g <- function(s, j) {
  v <- rep(0.1, 6)
  v[j] <- s
  accept <- mwg(b, X, mu, sigma, scale = v, nsim = 10^3)$accept
  (accept[3, j] - 0.234)^2
}

find_scale <- function(j) {
  # Optimise the (stochastic) objective function 5 times then average
  runs <- replicate(5, optimise(f = g, c(0, 0.1), j = j)$minimum)
  k <- mean(runs) / se[j]
  k * se # Guess at optimal scaling
}

opt_scale1 <- find_scale(1) # And so on..
opt_scale2 <- find_scale(2)
opt_scale3 <- find_scale(3)

# This seems not to be working well

my_guess <- c(0.175, 0.00025, 0.2, 0.3, 0.005, 0.2)

test2 <- mwg(b, X, mu, sigma, scale = my_guess, nsim = 10^4)
mcmc_trace(test2$chain)
test2$accept

# Larger sample size and starting at zero
full <- mwg(rep(0, 6), X, mu, sigma, scale = my_guess, nsim = 10^5)
mcmc_trace(full$chain) # 1, 2, 5 could still be improved but looking better
full$accept

saveRDS(full, "results/full_model.Rds")

colMeans(full$chain); b # Pretty close to the ML estimates, so the code is probably correct

# Subset Models -----------------------------------------------------------

head(X) # Reminder about the variables we have
# Only using a subset of the predictors in each model

# |   | Intercept | b1 | b2 | b3 | b4 | b5 |
# |---|-----------|----|----|----|----|----|
# | 1 |     Y     | Y  | N  | N  | Y  | Y  |
# | 2 |     Y     | N  | Y  | Y  | N  | Y  |

# i.e. link parameter is b5, coefficient corresponding to phc

X1 <- X[, c(1, 2, 5, 6)] # Design matrix 1: intercept, age, green and phc
p1 <- 4

X2 <- X[, c(1, 3, 4, 6)] # Design matrix 2: intercept, netuse, treated and phc
p2 <- 4

# GLM for model 1
mle1 <- glm(Y ~ X1 - 1, family = "binomial")
summary(mle1)

b1 <- mle1$coefficients

# GLM for model 2
mle2 <- glm(Y ~ X2 - 1, family = "binomial")
summary(mle2)

b2 <- mle2$coefficients

# Priors same as in the joint modelling to start with (so consistent in the link parameter)
mu1 <- rep(0, p1)
sigma1 <- rep(0.1, p1)

mu2 <- rep(0, p2)
sigma2 <- rep(0.1, p2)

# Run the samplers
model1 <- mwg(b0 = rep(0, 4), X = X1, mu = mu1, sigma = sigma1,
               scale = my_guess[c(1, 2, 5, 6)], nsim = 10^5)
names(model1$chain) <- c("b0", "b1", "b4", "b5") # Correct the naming
saveRDS(model1, "output/model1.Rds")

model2 <- mwg(b0 = rep(0, 4), X = X2, mu = mu2, sigma = sigma2,
               scale = my_guess[c(1, 3, 4, 6)], nsim = 10^5)
names(model2$chain) <- c("b0", "b2", "b3", "b5") # Correct the naming
saveRDS(model2, "output/model2.Rds")

mcmc_trace(model1$chain) # Traceplots
mcmc_hist(model1$chain) # Histograms
model1$accept
colMeans(model1$chain); b1

mcmc_trace(model2$chain) # Traceplots
mcmc_hist(model2$chain) # Histograms
model2$accept
colMeans(model2$chain); b2

# Markov Combination ------------------------------------------------------

logtarget <- function(b, s1, s2, X1, X2, mu1, mu2, sigma1, sigma2) {
  b1 <- b[s1]; b2 <- b[s2]
  as.numeric(logpost(b1, X1, mu1, sigma1) +
             logpost(b2, X2, mu2, sigma2) -
             sum((b1[c(1, 4)] - mu1[c(1, 4)])^2 / 2*sigma1[c(1, 4)]))
}

comb_mwg <- function(b0, s1, s2, X, mu, sigma, scale, nsim) { # Systematic-scan
  X1 <- X[, s1]; X2 <- X[, s2]

  mu1 <- mu[s1]; mu2 <- mu[s2]

  sigma1 <- sigma[s1]; sigma2 <- sigma[s2]

  p <- length(b0)
  r <- array(NA, c(nsim, p))
  r[1, ] <- b0 # Init chain

  order <- c(setdiff(s1, s2), setdiff(s2, s1), intersect(s1, s2))

  for (i in 2:nsim) {
    j <- order[((i - 2) %% p) + 1] # Index of the parameter to be updated

    if(j %in% setdiff(s1, s2)) {
      y <- r[i-1, ]
      y[j] <- y[j] + rnorm(1, mean = 0, sd = scale[j])
      a <- exp(logpost(y[s1], X1, mu1, sigma1) - logpost(b[s1], X1, mu1, sigma1))
      if(a > runif(1)) {
        r[i, ] <- y
      } else {
        r[i, ] <- r[i-1, ]
      }
    }

    if(j %in% setdiff(s2, s1)) {
      y <- r[i-1, ]
      y[j] <- y[j] + rnorm(1, mean = 0, sd = scale[j])
      a <- exp(logpost(y[s2], X2, mu2, sigma2) - logpost(b[s2], X2, mu2, sigma2))
      if(a > runif(1)) {
        r[i, ] <- y
      } else {
        r[i, ] <- r[i-1, ]
      }
    }

    if(j %in% intersect(s2, s1)) {
      y <- r[i-1, ]
      y[j] <- y[j] + rnorm(1, mean = 0, sd = scale[j])
      a <- exp(logtarget(y, s1, s2, X1, X2, mu1, mu2, sigma1, sigma2) -
               logtarget(b, s1, s2, X1, X2, mu1, mu2, sigma1, sigma2))
      if(a > runif(1)) {
        r[i, ] <- y
      } else {
        r[i, ] <- r[i-1, ]
      }
    }
  }
  r <- as.data.frame(r)
  names(r) <- sprintf("b%d", 0:(p-1))
  return(r)
}

b0 <- rep(0, 6)
s1 <- c(1, 2, 5, 6)
s2 <- c(1, 3, 4, 6)

test <- comb_mwg(b0, s1, s2, X, mu, sigma, scale = my_guess, nsim = 10^5)
saveRDS(test, "output/comb.Rds")
colMeans(test)
mcmc_trace(test)
mcmc_hist(test)
