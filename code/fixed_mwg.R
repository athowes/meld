library(mcmcse)

set.seed(1) # For reproducibility

# Prior on kappa
kappa0 <- -0.1
upsilon0 <- 0.5

# Prior on lambda
lambda0 <- 0
gamma0 <- 0.25

# Sample sizes
N2 <- 15
n2 <- 10

# True values for simulation
sim_kappa <- -0.25
sim_lambda <- -1
L <- sim_lambda - sim_kappa/2 # Lower
U <- sim_lambda + sim_kappa/2 # Upper

qC <- exp(L) / (1 + exp(L)) # Probabilities
qT <- exp(U) / (1 + exp(U))

# Simulate data
sC <- rbinom(n2, N2, qC)
sT <- rbinom(n2, N2, qT)

# Logposterior
logpost <- function(kappa, lambda, sC, sT) {
  L <- lambda - kappa/2
  U <- lambda + kappa/2
  qC <- exp(L) / (1 + exp(L))
  qT <- exp(U) / (1 + exp(U))
  sum1 <- sum(sC) * log(qC) + (2 + n2*N2 - sum(sC)) * log(1 - qC) + sum(sT) * log(qT) + (2 + n2*N2 - sum(sT)) * log(1 - qT)
  sum2 <- -0.5 * upsilon0 * (kappa - kappa0)^2 - 0.5 * gamma0 * (lambda - lambda0)^2
  return(sum1 + sum2)
}

# Function to update jth component of x via RWM step
mh <- function(x, j, scale) {
  # x = (kappa, lambda)
  # j in (1, 2)
  xs <- x # x is current, xs is x^star
  xs[j] <- xs[j] + rnorm(1, mean = 0, sd = scale)
  # Acceptance probability
  ls <- logpost(xs[1], xs[2], sC, sT)
  l <- logpost(x[1], x[2], sC, sT)
  a <- exp(ls - l)
  if(a > runif(1)) return(xs) else return(x)
}

# Metropolis-within-Gibbs sampler (random scan by default)
mwg <- function(x0, vscale, nsim, random_scan = TRUE) {
  p <- length(x0)
  r <- array(NA, c(nsim, p)) # For the chain
  r[1, ] <- x0 # Init chain
  s <- array(0, c(3, p)) # For acceptance rates
  for(i in 2:nsim) {
    if(random_scan) {
      j <- sample(1:p, 1) # Random scan
    } else {
      j <- (i %% p) + 1 # Systematic scan
    }
    s[1, j] <- s[1, j] + 1 # Update pick count
    r[i, ] <- mh(r[i-1, ], j, vscale[j])
    if(!identical(r[i, ], r[i-1, ])) s[2, j] <- s[2, j] + 1 # Update accept count
  }
  r <- as.data.frame(r)
  s[3, ] <- s[2, ] / s[1, ] # Acceptance rates
  return(list("chain" = r, "accept" = s))
}

x0 <-c(-0.1, 0) # Init at a-priori means
vscale <- c(0.55, 0.3) # Tuning to 0.44 (Neal and Roberts 2006)

run <- mwg(x0, vscale, nsim = 1000000) # First run

minESS <- minESS(p = 2, eps = .025, alpha = .025) # (Vats, Felgal and Jones 2015)
ESS <- multiESS(run$chain)
ESS/minESS # About 2.5x the required effective samples

saveRDS(c(sC, sT), "../output/data_fixed.Rds") # Save data
saveRDS(c(sim_kappa, sim_lambda), "../output/truth_fixed.Rds") # Ground truth
saveRDS(run, "../output/mcmc_fixed.Rds") # MCMC output

kappa_mean <- function(run) {
  colMeans(run$chain)[1]
}

fixed_means <- replicate(100, kappa_mean(mwg(x0, vscale, nsim = 1000000)))
saveRDS(fixed_means, "../output/fixed_means.Rds")
