# Sub-model 1
chains1 <- readRDS("results/chains1.Rds")

for (i in c(1, 2, 5, 6)) {
  assign(paste0("m1_plot", i), qplot(x = 1:100000, y = chains1[[paste0("b", i)]], geom = "line"))
}

cowplot::plot_grid(m1_plot1, m1_plot2, m1_plot5, m1_plot6) # These all look good

colMeans(chains1)

# Sub-model2
chains2 <- readRDS("results/chains2.Rds")

for (i in c(1, 3, 4, 6)) {
  assign(paste0("m2_plot", i), qplot(x = 1:100000, y = chains1[[paste0("b", i)]], geom = "line"))
}

cowplot::plot_grid(m2_plot1, m2_plot3, m2_plot4, m2_plot6) # TAlso look good

colMeans(chains2)