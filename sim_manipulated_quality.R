
library(tidyverse)
theme_set(theme_bw())


N <- 10000      # Number of simulations
G <- 3          # Threshold number of test a technician will take
c <- 0.08       # cut off value where technician will repeat the test
mux <- 0.03     # mean of process
sdx <- 0.056    # irreducible error
sdz <- 0.008    # measurement error


Ygx <- function(x,c,s,G) {
  Z <- (c-x)/s
  p <- pnorm( Z )
  return(x - ( s*dnorm( Z )/p ) * ( 1-(1-p)^G))
}



# Plot of the expected process <- not in Journal, only for reference.
data.frame(y=sapply(rnorm(N,mean=mux,sd=sdx), function(x) Ygx(x,c,sdz,G))) %>%
  ggplot(aes(x=y,y=..density..)) +
  geom_histogram(alpha=0.7,bins=50) +
  geom_vline(xintercept=c,color='red',linewidth=0.8)



X <- rnorm(N,mean=mux,sd=sdx)
Z <- lapply(1:N,function(x) rnorm(G,mean=0,sd=sdz))
y <- numeric(N)
for ( i in 1:N) {
  idx <- which(Z[[i]]+X[i]<c)
  if(length(idx)!=0) {
    y[i] <- X[i] + Z[[i]][min(idx)]
  } else {
    y[i] <- X[i] + Z[[i]][G]
  }
}


data.frame(y=y) %>%
  ggplot(aes(x=y,y=..density..)) +
  geom_histogram(alpha=0.7,bins=80) +
  geom_vline(xintercept=c,color='red',linewidth=0.8) +
  theme(text = element_text(size = plot_scale_factor*10)) +
  xlim(-0.25,0.25)


ggsave(
  filename = "s4_quality_distribution_simulation.svg",
  plot = last_plot(),
  device = "svg",
  path = "./results/charts",
  scale = plot_scale_factor,
  width = 12,
  height = 6.75,
  units = "cm",#c("in", "cm", "mm", "px"),
  dpi = 1000
)


