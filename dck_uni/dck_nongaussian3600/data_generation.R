library(parallel)
library(geoR)
library(MASS)
library(fields)

setwd("/Users/beauniverse/Documents/GitHub/Deep-classifier-kriging-for-probabilistic-spatial-prediction-of-air-quality-index/dck_uni/dck_nongaussian3600")
mainDir <- "."
subDir <- "synthetic_data_simulations/"
trainDir <- paste0(subDir, "training_data/")
testDir <- paste0(subDir, "testing_data/")
DNNDir <- paste0(subDir, "Results_DNN/")
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
dir.create(file.path(mainDir, trainDir), showWarnings = FALSE)
dir.create(file.path(mainDir, testDir), showWarnings = FALSE)
dir.create(file.path(mainDir, DNNDir), showWarnings = FALSE)

set.seed(12345567)
num_sim <- 100
g <- 0.8
h <- 0.5
x <- seq(0, 1, length.out = 80)
y <- seq(0, 1, length.out = 80)
d1 <- expand.grid(x = x, y = y)
sample_idx <- sample(1:6400, 3600)
d1 <- d1[sample_idx, ]
X <- d1$x
Y <- d1$y
m <- as.matrix(dist(data.frame(X = X, Y = Y)))

s11 <- 0.7; nu11 <- 0.3; alpha11 <- 0.05
matern_cov1 <- s11 * matern(m, sqrt(alpha11), nu11)
a <- 0.1; s <- 0.9; nu <- 0.5
cov <- s * matern(m, a, nu)

run_simulation <- function(sim) {
  x1 <- mvrnorm(1, rep(sample(seq(-2.5, 2.5, length.out = 100), 1), 3600), cov)
  x2 <- mvrnorm(1, rep(sample(seq(-2.5, 2.5, length.out = 100), 1), 3600), cov)
  x3 <- mvrnorm(1, rep(sample(seq(-2.5, 2.5, length.out = 100), 1), 3600), cov)
  x4 <- mvrnorm(1, rep(sample(seq(-2.5, 2.5, length.out = 100), 1), 3600), cov)
  x5 <- mvrnorm(1, rep(sample(seq(-2.5, 2.5, length.out = 100), 1), 3600), cov)
  
  mean_val <- x1^2 - x2^2 + x3^2 - x4^2 - x5^2 + 2*x1*x2 + 3*x2*x3 - 2*x3*x5 +
    10*x1*x4 + sin(x1)*x2*x3 + cos(x2)*x3*x5 + (x1*x2*x4*x5)
  
  simulation <- mvrnorm(1, rep(0, 3600), matern_cov1)
  tukey_var1 <- if (g == 0) simulation * exp(h * simulation^2 / 2) + mean_val else
    (exp(g * simulation) - 1) * exp(h * simulation^2 / 2) / g + mean_val
  
  df <- data.frame(
    x = X, y = Y,
    X1 = x1, X2 = x2, X3 = x3, X4 = x4, X5 = x5,
    var1 = tukey_var1
  )
  
  write.csv(
    df,
    file = paste0(subDir, "sim3600_g", g, "_h", h, "_", sim, ".csv"),
    row.names = FALSE
  )
  return(NULL)
}

# Detect number of cores and run in parallel
num_cores <- detectCores() - 1
mclapply(1:num_sim, run_simulation, mc.cores = num_cores)



