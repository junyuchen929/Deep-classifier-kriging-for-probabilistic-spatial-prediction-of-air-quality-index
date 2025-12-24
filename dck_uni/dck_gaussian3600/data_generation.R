rm(list = ls())
library(geoR)
library(MASS)
library(fields)
setwd("/Users/beauniverse/Documents/GitHub/Deep-classifier-kriging-for-probabilistic-spatial-prediction-of-air-quality-index/dck_uni/dck_gaussian3600")
mainDir <- "."
subDir <- "synthetic_data_simulations/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

Dtgh1 <- function(nx, theta, g, h) {
  grid <- expand.grid(1:nx, 1:nx)
  xy <- matrix(runif(n * 2, -0.4, 0.4), n, 2)
  loc <- (grid - 0.5 + xy) / nx
  x <- loc[, 1]
  y <- loc[, 2]
  
  dm <- as.matrix(dist(loc))
  omega1 <- exp(-dm/theta)
  z <- t(chol(omega1))%*%rnorm(n)
  
  if(g==0){tgh = z*exp(h*z^2/2)}
  else {tgh = (exp(g*z)-1)*exp(h*z^2/2)/g}
  out = data.frame(x, y, tgh)
  return(list(simulation=out, Sigma22=omega1))
}

set.seed(16739)
theta <- 0.5 
sim_size <- 100

g_values <- c(0)
h_values <- c(0)
data_list  <- vector("list", length(g_values))

nx <- 60
n <- nx^2


for (c in 1:length(g_values)) {
  data_list[[c]]  <- vector("list", sim_size)
  
  for (i in 1:sim_size) {
    g <- g_values[c]
    h <- h_values[c]
    
    res <- Dtgh1(nx, theta = theta, g = g, h = h)
    dat <- res$simulation  
    
    data_list[[c]][[i]]  <- dat
  }
}

for (c in 1:length(g_values)) {
  for (i in 1:sim_size) {
    filename <- paste0("synthetic_data_simulations/Tgh_sim_g", g_values[c],
                       "_h", h_values[c], "_", i, ".csv")
    write.csv(data_list[[c]][[i]], file = filename, row.names = FALSE)
  }
}

mainDir <- "./synthetic_data_simulations/"
subDir <- "training_data/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
subDir <- "testing_data/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

mainDir <- "."
subDir <- "Results_DNN/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
