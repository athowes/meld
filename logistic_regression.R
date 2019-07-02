# Bayesian logistic regression example

options(scipen = 999)
pacman::p_load(tidyverse, rjags, geoR)

data(gambia)
X <- as.matrix(gambia[, c(4:8)])
Y <- gambia[, 3] # Response variable
n <- length(Y)

mle <- glm(Y ~ X, family = "binomial") # Logistic regression using MLE
summary(mle)

b <- mle$coefficients

X <- cbind(intercept = 1, X) # Add intercept column to design matrix

# Full model --------------------------------------------------------------

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
  loglike <- sum(nu[Y == 1]) + sum(-log(1 + exp(nu)))
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
    r[i, ] <- mh(r[i-1, ], X, mu, sigma, j, scale)
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

# Scaling and running -----------------------------------------------------

trace_plots <- function(x, nsim) { # Helper function
  plot1 <- qplot(x = 1:nsim, y = x$chain$b1, geom = "line")
  plot2 <- qplot(x = 1:nsim, y = x$chain$b2, geom = "line")
  plot3 <- qplot(x = 1:nsim, y = x$chain$b3, geom = "line")
  plot4 <- qplot(x = 1:nsim, y = x$chain$b4, geom = "line")
  plot5 <- qplot(x = 1:nsim, y = x$chain$b5, geom = "line")
  plot6 <- qplot(x = 1:nsim, y = x$chain$b6, geom = "line")
  cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 3)
}

test <- mwg(b, X, mu, sigma, scale = rep(0.1, 6), nsim = 10^4) # Initialising at roughly the right values
trace_plots(test, 10^4)
test$accept # Plan is to adjust scale based on optimal acceptance rate 0.234

# Optimal scalings for each parameter should roughly be some constant multiplied by their standard error
se <- sqrt(diag((vcov(mle))))

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
trace_plots(test2, 10^4)
test2$accept

# Larger sample size and starting at zero
full <- mwg(rep(0, 6), X, mu, sigma, scale = my_guess, nsim = 10^5)
trace_plots(full, 10^5) # 1, 2, 5 could still be improved but looking better
full$accept

saveRDS(full, "results/full_model.Rds")

colMeans(full$chain); b # Pretty close to the ML estimates, so the code is probably correct

# Subset models -----------------------------------------------------------

head(X) # Reminder about the variables we have
# Only using a subset of the predictors in each model

# |   | Intercept | b1 | b2 | b3 | b4 | b5 |
# |---|-----------|----|----|----|----|----|
# | 1 |     Y     | Y  | N  | N  | Y  | Y  |
# | 2 |     Y     | N  | Y  | Y  | N  | Y  |

# i.e. link parameter is b5, coefficient corresponding to phc

X1 <- X[, c(1, 2, 5, 6)] # Design matrix 1: intercept, age, green and phc
p1 <- 4

X2 <- X[, c(1, 3, 4, 6)] # Design matrix 2: intercept netuse, treated and phc
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
names(model1$chain) <- c("b01", "b1", "b4", "b5") # Correct the naming
saveRDS(model1, "results/model1.Rds")

model2 <- mwg(b0 = rep(0, 4), X = X2, mu = mu2, sigma = sigma2,
               scale = my_guess[c(1, 3, 4, 6)], nsim = 10^5)
names(model2$chain) <- c("b02", "b2", "b3", "b5") # Correct the naming
saveRDS(model2, "results/model2.Rds")

colMeans(model1$chains); b1

colMeans(model2$chains); b2 # Still getting close: good!

# Markov combination ------------------------------------------------------

# 4.1. equation (*) can be calculated by logpost(model1) + logpost(model2) - logprior(beta5)
# This is only used for the link parameter update
logtarget <- function(b, X1, X2, mu1, sigma1, mu2, sigma2) {
  b1 <- b[c(1, 2, 5, 6)]
  b2 <- b[c(1, 3, 4, 6)]
  as.numeric(logpost(b1, X1, mu1, sigma1) + 
             logpost(b2, X2, mu2, sigma2) - 
             ((b[6] - mu1[4])^2 / 2*sigma1[4])) # Could use either mu1 or mu2 here
}

logtarget(b, X1, X2, mu1, sigma1, mu2, sigma2)