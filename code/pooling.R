# Reproduce Figure 4 from Melding paper (probably very inefficiently)

library(ggplot2)
library(wesanderson)

p <- c("#E69F00", "#56B4E9", "#009E73")

dlin_pool <- function(x, mean, sd, w, k) { # Linear pool two densities, equal weighting
  k * (w[1] * dnorm(x, mean[1], sd[1]) + w[2] * dnorm(x, mean[2], sd[2]))
}

dlog_pool <- function(x, mean, sd, w, k) { # Logarithmic pool two densities, equal weighting
  k * dnorm(x, mean[1], sd[1])^w[1] * dnorm(x, mean[2], sd[2])^w[2]
}

dPoE_pool <- function(x, mean, sd, k) { # PoE pool two densities, equal weighting
  k * dnorm(x, mean[1], sd[1]) * dnorm(x, mean[2], sd[2])
}

pool_plot <- function(mean, sd, w) {
  k1 <- 1/integrate(dPoE_pool, mean = c(-1, 1), sd = c(1, 1),
                    k = 1, lower = -25, upper = 25)$value
  k2 <- 1/integrate(dlog_pool, mean = c(-1, 1), sd = c(1, 1),
                    w = w, k = 1, lower = -25, upper = 25)$value
  k3 <- 1/integrate(dlin_pool, mean = c(-1, 1), sd = c(1, 1),
                    w = w, k = 1, lower = -25, upper = 25)$value
  ggplot(data = data.frame(x = c(-5, 5)), aes(x)) +
    stat_function(fun = dnorm, n = 201,
                  args = list(mean = mean[1], sd = sd[1]),
                  lty = "dashed", lwd = .3) +
    stat_function(fun = dnorm, n = 201,
                  args = list(mean = mean[2], sd = sd[2]),
                  lwd = .3) +
    stat_function(fun = dPoE_pool, n = 201, args = list(mean, sd, k1),
                  geom = "area", fill = p[3], alpha = 0.5) +
    stat_function(fun = dlog_pool, n = 201, args = list(mean, sd, w, k2),
                  geom = "area", fill = p[2], alpha = 0.5) +
    stat_function(fun = dlin_pool, n = 201, args = list(mean, sd, w, k3),
                  geom = "area", fill = p[1], alpha = 0.5) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "grey"))
}

p1 <- pool_plot(mean = c(-1, 1), sd = c(1, 1), w = c(0.5, 0.5))
p2 <- pool_plot(mean = c(-1, 1), sd = c(1, 1), w = c(0.75, 0.25))

cowplot::plot_grid(p1, p2, ncol = 2)
