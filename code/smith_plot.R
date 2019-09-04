library(ggplot2)

run <- readRDS("../output/mcmc_smith.Rds") # Import MCMC
truth <- readRDS("../output/truth_smith.Rds") # Import ground truth

midblue <- "#3D9BD0"
midgreen <- "#00855A"
midpink <- "#B3608E"

nsim <- dim(run$chain)[1] # Length of chain
thinfactor <- 1000 # Thinning for trace-plot display
thin <- run$chain[thinfactor*1:(nsim / thinfactor), ]

trace_hist <- function(num, latex, mycol = midblue) {
  p1 <- ggplot(thin, aes_string(y = paste0("V", num), x = 1:(nsim/thinfactor))) +
    geom_line(alpha = 0.8, col = "grey") +
    labs(x = "Iterations", y = latex) +
    scale_x_continuous(breaks = c(0, (nsim/(2*thinfactor)), (nsim/thinfactor)), labels = c(0, (nsim/2), nsim)) +
    theme(panel.background = element_blank(),
          axis.line = element_line(colour = "grey"),
          plot.margin = unit(c(5.5, 7, 5.5, 5.5), "points"))

  p2 <- ggplot(run$chain, aes_string(x = paste0("V", num))) +
    geom_histogram(alpha = 0.6, bins = 30, fill = mycol, col = mycol) +
    labs(x = latex, y = "Count") +
    geom_vline(xintercept = truth[num], colour = "black", linetype = "longdash") +
    theme(panel.background = element_blank(),
          axis.line = element_line(colour = "grey"))
  return(cowplot::plot_grid(p1, p2, ncol = 2))
}

trace_hist(1, "$\\kappa$")
trace_hist(2, "$\\upsilon$")

trace_hist(3, "$p_1^C$")
trace_hist(8, "$p_1^T$")

trace_hist(4, "$p_2^C$")
trace_hist(9, "$p_2^T$")

trace_hist(5, "$p_3^C$")
trace_hist(10, "$p_3^T$")

trace_hist(6, "$p_4^C$")
trace_hist(11, "$p_4^T$")

trace_hist(7, "$p_5^C$")
trace_hist(12, "$p_5^T$")
