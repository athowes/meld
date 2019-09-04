set.seed(1)

# Import data sources
y1 <- readRDS("../output/data_smith.Rds")
y2 <- readRDS("../output/data_fixed.Rds")

# Sample sizes
n1 <- 5
N1 <- 20

n2 <- 10
N2 <- 15

# Extract data
rC <- head(y1, n = n1)
rT <- tail(y1, n = n1)

sC <- head(y2, n = n2)
sT <- tail(y2, n = n2)

# Priors from Submodel 1
kappa01 <- 0 # Add the 1 index to distinguish from the kappa0 specified in Submodel 2
upsilon01 <- 0.1
a <- 3
b <- 1
mu0 <- 0
tau0 <- 0.25

# Priors form Submodel 2
kappa02 <- -0.25 # Likewise with the 2 index
upsilon02 <- 1
lambda0 <- 0
gamma0 <- 0.25

# Logposterior from Submodel 1
l1 <- function(kappa, upsilon, pC, pT, rC, rT) {
  if(all(pT > 0 & pT < 1) & all(pC > 0 & pC < 1) & upsilon > 0) {
    omega <- (tau0^-1 + 0.25 * upsilon^-1)^-1
    sum1 <- sum(((rC - 1) * log(pC)) + ((N1 - rC - 1) * log(1 - pC))) + sum(((rT - 1) * log(pT)) + ((N1 - rT - 1) * log(1 - pT)))
    sum2 <- - 0.5 * omega * sum( (log(pC / (1 - pC)) - mu0 + kappa/2)^2 + (log(pT / (1 - pT)) - mu0 - kappa/2)^2 )
    sum3 <- 2 * n1 * omega - 0.5 * upsilon01 * (kappa - kappa01)^2 + (a - 1) * log(upsilon) - b * upsilon
    return(sum1 + sum2 + sum3)
  }
  else return(-Inf)
}

# Logposterior from Submodel 2
l2 <- function(kappa, lambda, sC, sT) {
  L <- lambda - kappa/2
  U <- lambda + kappa/2
  qC <- exp(L) / (1 + exp(L))
  qT <- exp(U) / (1 + exp(U))
  sum1 <- sum(sC) * log(qC) + (2 + n2*N2 - sum(sC)) * log(1 - qC) + sum(sT) * log(qT) + (2 + n2*N2 - sum(sT)) * log(1 - qT)
  sum2 <- -0.5 * upsilon02 * (kappa - kappa02)^2 - 0.5 * gamma0 * (lambda - lambda0)^2
  return(sum1 + sum2)
}

# Parameter format: (kappa, upsilon, pC, pT, lambda)
x0 <-c(0, 1/3, rep(0.5, n1), rep(0.5, n1), 0)
vscale <- c(0.6, 3.5, 0.3, 0.25, 0.3, 0.25, 0.2, 0.15, 0.15, 0.25, 0.3, 0.2, 0.4)

meld <- function(x0, vscale, nsim) { # Systematic-scan
  p <- length(x0)
  r <- array(NA, c(nsim, p))
  r[1, ] <- x0 # Init chain
  s <- array(0, c(3, p)) # For acceptance rates

  for (i in 2:nsim) {

    j <- (i %% p) + 1 # Systematic scan
    s[1, j] <- s[1, j] + 1 # Update pick count

    if(j == 1) { # phi update
      x  <- r[i-1, ]
      xs <- r[i-1, ]
      xs[j] <- xs[j] + rnorm(1, mean = 0, sd = vscale[j])
      ls <- l1(xs[1], xs[2], xs[3:7], xs[8:12], rC, rT) + l2(xs[1], xs[13], sC, sT) # Product of Experts pooling
      l <- l1(x[1], x[2], x[3:7], x[8:12], rC, rT) + l2(x[1], x[13], sC, sT)
      a <- exp(ls - l)
      if(a > runif(1)) {
        r[i, ] <- xs
        s[2, j] <- s[2, j] + 1 # Update accept count
      } else {
        r[i, ] <- x
      }
    }
    if(j %in% 2:12) { # psi1 update
      x  <- r[i-1, ]
      xs <- r[i-1, ]
      xs[j] <- xs[j] + rnorm(1, mean = 0, sd = vscale[j])
      ls <- l1(xs[1], xs[2], xs[3:7], xs[8:12], rC, rT)
      l <- l1(x[1], x[2], x[3:7], x[8:12], rC, rT)
      a <- exp(ls - l)
      if(a > runif(1)) {
        r[i, ] <- xs
        s[2, j] <- s[2, j] + 1 # Update accept count
      } else {
        r[i, ] <- x
      }
    }
    if(j == 13) { # psi2 update
      x  <- r[i-1, ]
      xs <- r[i-1, ]
      xs[j] <- xs[j] + rnorm(1, mean = 0, sd = vscale[j])
      ls <- l2(xs[1], xs[13], sC, sT)
      l <- l2(x[1], x[13], sC, sT)
      a <- exp(ls - l)
      if(a > runif(1)) {
        r[i, ] <- xs
        s[2, j] <- s[2, j] + 1 # Update accept count
      } else {
        r[i, ] <- x
      }
    }
  }
  r <- as.data.frame(r)
  s[3, ] <- s[2, ] / s[1, ] # Acceptance rates
  return(list("chain" = r, "accept" = s))
}

run <- meld(x0, vscale, nsim = 5000000)

minESS <- minESS(p = 13, eps = .025, alpha = .025) # (Vats, Felgal and Jones 2015)
ESS <- multiESS(run$chain)

truth <- readRDS("../output/truth_smith.Rds") # Import ground truth from Submodel 1
saveRDS(c(truth, -1), "../output/truth_meld.Rds") # Append GT from S2 and save

saveRDS(run, "../output/mcmc_meld.Rds") # MCMC output
