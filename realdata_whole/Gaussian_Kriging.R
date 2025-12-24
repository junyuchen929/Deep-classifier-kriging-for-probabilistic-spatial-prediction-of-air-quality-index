rm(list = ls())
source("LMC_functions_not_colocated.R")
library(sp)


data_train_part1 <- read.csv("CSV1_train.csv", header = TRUE)
data_train_part2_full <- read.csv("CSV2.csv", header = TRUE)
data_test <- read.csv("CSV1_test.csv", header = TRUE)


coords1 <- data_train_part1[, c("uj_lat", "uj_lon")]
hull_idx <- chull(coords1)
poly_coords <- coords1[c(hull_idx, hull_idx[1]), ]  
sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(poly_coords)), ID = "1")))
points2 <- SpatialPoints(data_train_part2_full[, c("uj_lat", "uj_lon")])
inside_idx <- !is.na(over(points2, sp_poly))
data_train_part2_inside <- data_train_part2_full[inside_idx, ]
cat("number of samples =", nrow(data_train_part2_inside), "\n")

set.seed(123)
sample_indices <- sample(1:nrow(data_train_part2_inside), 10000) 
data_train_part2 <- data_train_part2_inside[sample_indices, ]
#write.csv(data_train_part2, "CSV2_10000.csv", row.names = FALSE)

normalize_minmax <- function(x, a, b) (x - a) / (b - a)

uj_lat_all <- c(data_train_part1$uj_lat, data_train_part2$uj_lat)
uj_lon_all <- c(data_train_part1$uj_lon, data_train_part2$uj_lon)
min_lat <- min(uj_lat_all, na.rm = TRUE); max_lat <- max(uj_lat_all, na.rm = TRUE)
min_lon <- min(uj_lon_all, na.rm = TRUE); max_lon <- max(uj_lon_all, na.rm = TRUE)

data_train_part1$uj_lat <- normalize_minmax(data_train_part1$uj_lat, min_lat, max_lat)
data_train_part1$uj_lon <- normalize_minmax(data_train_part1$uj_lon, min_lon, max_lon)
data_train_part2$uj_lat <- normalize_minmax(data_train_part2$uj_lat, min_lat, max_lat)
data_train_part2$uj_lon <- normalize_minmax(data_train_part2$uj_lon, min_lon, max_lon)
data_test$uj_lat <- normalize_minmax(data_test$uj_lat, min_lat, max_lat)
data_test$uj_lon <- normalize_minmax(data_test$uj_lon, min_lon, max_lon)

data_train <- list(part1 = data_train_part1, part2 = data_train_part2)


start_time <- Sys.time()

init.ind <- c(1, 1, 1, 1, 1, 1, 0, 0)
indmat.estim <- optim_indmat_loglik(init.ind, data_train)

d_med = 0.44

a1_init <- d_med        
a2_init <- d_med        

nu1_init <- 1.5          
nu2_init <- 1.5
sigma1_init <- var(data_train$part1$Z1)                    
sigma2_init <- var(c(data_train$part1$Z2, data_train$part2$Z2))


b11_init <- 1; b12_init <- 0
b21_init <- 0; b22_init <- 1

nug1_init <- 0.05 * sigma1_init   
nug2_init <- 0.05 * sigma2_init

init.lmc <- c(a1_init, nu1_init, sigma1_init,
              a2_init, nu2_init, sigma2_init,
              b11_init, b12_init, b21_init, b22_init,
              nug1_init, nug2_init)


fit.Model.lmc <- optim_lmc.loglikelihood(par = init.lmc, data_train)



pred <- lmc.pred.summary(fit.Model.lmc$par, data_train, data_test, var_type = "obs")

end_time <- Sys.time()
comp_time <- as.numeric(difftime(end_time, start_time, units = "secs"))


saveRDS(pred, "pred_Gaussian.rds")
pred$pred
diag(pred$conditional_var)

cat("Computation time (in seconds):", comp_time, "\n")


