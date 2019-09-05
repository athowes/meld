library(ggplot2)

midblue <- "#3D9BD0"
midgreen <- "#00855A"
midpink <- "#B3608E"

smith <- as.numeric(readRDS("../output/means_smith.Rds"))
fixed <- as.numeric(readRDS("../output/means_fixed.Rds"))
meld <- as.numeric(readRDS("../output/means_meld.Rds"))

df <- data.frame(method = c(rep("meld", 100), rep("fixed", 100), rep("smith", 100)), kappa = c(meld, fixed, smith))

ggplot(df, aes(x = kappa)) +
  geom_boxplot(aes(fill = method), alpha = 0.5, position = "dodge") +
  geom_vline(xintercept = -0.25, colour = "black", linetype = "longdash") +
  scale_fill_manual(values = c(midgreen, midpink, midblue),
                    name = "",
                    labels = c("Fixed Effects", "Melded", "Random Effects")) +
  labs(x = "$\\kappa$", y = "Count") +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "grey"))

