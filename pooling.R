# Reproduce Figure 4 from Melding paper

library(ggplot2)
library(RColorBrewer)

display.brewer.pal(3, "Dark2")
brewer.pal(3, "Dark2")

dlin_pool <- function(x, mean, sd) { # Linear pool two densities, equal weighting
  0.5 * (dnorm(x, mean[1], sd[1]) + dnorm(x, mean[2], sd[2])) 
}

dlog_pool <- function(x, mean, sd) { # Logarithmic pool two densities, equal weighting
  k <- 2
  k * dnorm(x, mean[1], sd[1])^0.5 * dnorm(x, mean[2], sd[2])^0.5
} # Not sure how to calculate normalising constant k

ggplot(data = data.frame(x = c(-5, 5)), aes(x)) +
  stat_function(fun = dnorm, n = 201, args = list(mean = -1, sd = 1)) +
  stat_function(fun = dnorm, n = 201, args = list(mean = 1, sd = 1)) +
  stat_function(fun = dlin_pool, n = 201, args = list(mean = c(-1, 1), sd = c(1, 1)),
                geom = "area", fill = "#1B9E77", alpha = 0.5) +
  stat_function(fun = dlog_pool, n = 201, args = list(mean = c(-1, 1), sd = c(1, 1)),
                geom = "area", fill = "#D95F02", alpha = 0.5)