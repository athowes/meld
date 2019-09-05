library(mcmcse)

set.seed(1) # For reproducibility

# Prior on kappa
kappa0 <- 0
upsilon0 <- 0.1

# Prior on upsilon
a <- 3
b <- 1

# Prior on mu
mu0 <- 0
tau0 <- 0.25

# Sample sizes
N1 <- 20
n1 <- 5

# True values for simulation
sim_kappa <- -0.25
sim_upsilon <- 4
sim_mu <- -1
sim_tau <- 2

# Simulate data
delta <- rnorm(n1, sim_kappa, sqrt(1/sim_upsilon))
mu <- rnorm(n1, sim_mu, sqrt(1/sim_tau))

L <- mu - delta/2 # Lower
U <- mu + delta/2 # Upper

pC <- exp(L) / (1 + exp(L)) # Probabilities
pT <- exp(U) / (1 + exp(U))

rC <- rbinom(n1, N1, pC)
rT <- rbinom(n1, N1, pT)

# Logposterior
logpost <- function(kappa, upsilon, pC, pT, rC, rT) {
  if(all(pT > 0 & pT < 1) & all(pC > 0 & pC < 1) & upsilon > 0) {
  omega <- (tau0^-1 + 0.25 * upsilon^-1)^-1
  sum1 <- sum(((rC - 1) * log(pC)) + ((N1 - rC - 1) * log(1 - pC))) + sum(((rT - 1) * log(pT)) + ((N1 - rT - 1) * log(1 - pT)))
  sum2 <- - 0.5 * omega * sum( (log(pC / (1 - pC)) - mu0 + kappa/2)^2 + (log(pT / (1 - pT)) - mu0 - kappa/2)^2 )
  sum3 <- 2 * n1 * omega - 0.5 * upsilon0 * (kappa - kappa0)^2 + (a - 1) * log(upsilon) - b * upsilon
  return(sum1 + sum2 + sum3)
  }
  else return(-Inf)
}

# Function to update jth component of x via RWM step
mh <- function(x, j, scale) {
  # x = (kappa, upsilon, pC, pT)
  # j in (1, ..., 2 + 2*n_1)
  xs <- x # x is current, xs is x^star
  xs[j] <- xs[j] + rnorm(1, mean = 0, sd = scale)
  # Acceptance probability
  ls <- logpost(xs[1], xs[2], xs[3:(2+n1)], xs[(3+n1):(2+2*n1)], rC, rT)
  l <- logpost(x[1], x[2], x[3:(2+n1)], x[(3+n1):(2+2*n1)], rC, rT)
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

# Geyer: no need for burn-in so long as you don't mind starting value in your chain http://users.stat.umn.edu/~geyer/mcmc/burn.html
x0 <-c(0, 1/3, rep(0.5, n1), rep(0.5, n1)) # Init at sensible values

# Tuning to 0.44 (Neal and Roberts 2006)
pC_scale <- c(0.3, 0.25, 0.3, 0.25, 0.2)
pT_scale <- c(0.15, 0.15, 0.25, 0.3, 0.2)
vscale <- c(3, 3.5, pC_scale, pT_scale)

run <- mwg(x0, vscale, nsim = 5000000)

minESS <- minESS(p = 12, eps = .025, alpha = .025) # (Vats, Felgal and Jones 2015)
ESS <- multiESS(run$chain)
ESS/minESS # About 1.5x the required effective samples

saveRDS(c(rC, rT), "../output/data_smith.Rds") # Save data
saveRDS(c(sim_kappa, sim_upsilon, pC, pT), "../output/truth_smith.Rds") # Ground truth
saveRDS(run, "../output/mcmc_smith.Rds") # MCMC output

# Replicates for posterior mean
kappa_mean <- function(run) {
  colMeans(run$chain)[1]
}

means <- replicate(100, kappa_mean(mwg(x0, vscale, nsim = 5000000)))
saveRDS(means, "../output/means_smith.Rds")
