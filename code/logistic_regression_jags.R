# Bayesian logistic regression example with JAGS

options(scipen = 999)
pacman::p_load(tidyverse, rjags, geoR)

data(gambia)
X <- as.matrix(gambia[, c(4:8)])
dim(X) # 2035 observations of 5 variables
Y <- gambia[, 3] # Response variable
n <- length(Y)
table(Y) # 1308 without and 727 with

mle <- glm(Y ~ X, family = "binomial") # Logistic regression using MLE
summary(mle)

b <- mle$coefficients

# Specify the JAGS model
# https://www4.stat.ncsu.edu/~reich/ABA/code/GLM

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