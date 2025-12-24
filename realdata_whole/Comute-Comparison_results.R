rm(list = ls())
source("LMC_functions_not_colocated.R")
source("density_computation_DNN.R")
library(VGAM)

pred_list_all <- readRDS("pred_Gaussian.rds")

C = 12
data_test = read.csv("projection_matrix_test.csv", header = TRUE)

data_test_dnn = read.csv("predictions_results.csv", header = TRUE)

data_train = read.csv("projection_matrix_train.csv", header = TRUE)

sigma <- mad(data_train$Z1, constant = 1)
n <- nrow(data_train)
alpha <- 1/3
h <- 2 * C * sigma * n^(-alpha) / 3

total_classes <- ncol(data_test_dnn) - 11
Z1_node <- matrix(NA, nrow = 1, ncol = total_classes)
Z2_node <- matrix(NA, nrow = 1, ncol = total_classes)
for (i in 1:total_classes) {
  df_sample <- data_train[data_train$class == i, ]
  Z1_node[, i] <- mean(df_sample$Z1)
  Z2_node[, i] <- mean(df_sample$Z2)
}

pred_obj = pred_list_all

n_test <- 97

library(doParallel)
library(foreach)

n <- n_test

pred1 <- pred_obj$pred[seq_len(n), 1]
pred2 <- pred_obj$pred[(n + seq_len(n)), 1]
v11 <- diag(pred_obj$conditional_var)[seq_len(n)]
v22 <- diag(pred_obj$conditional_var)[n + seq_len(n)]
v12 <- pred_obj$conditional_var[cbind(seq_len(n), n + seq_len(n))]
Z2_obs <- data_test$Z2

cl <- makeCluster(max(1L, parallel::detectCores() - 1L))
registerDoParallel(cl)


res <- foreach(i = seq_len(n), .combine = rbind,
               .packages = c()) %dopar% {
                 ro <- v12[i] / sqrt(v11[i] * v22[i])
                 cond_pred <- pred1[i] + ro * (v11[i] / v22[i]) * (Z2_obs[i] - pred2[i])
                 cond_var_scalar <- (1 - ro^2) * v11[i]
                 sd_i <- sqrt(max(cond_var_scalar, 0))
                 
                 g_med <- qnorm(0.5,  cond_pred, sd_i)
                 g_lb  <- qnorm(0.025, cond_pred, sd_i)
                 g_ub  <- qnorm(0.975, cond_pred, sd_i)
                 
                 d_med <- marginal.cond.quantile_z1(h, 0.5,   Z2_obs[i], i, data_test_dnn, Z1_node, Z2_node)
                 d_lb  <- marginal.cond.quantile_z1(h, 0.025, Z2_obs[i], i, data_test_dnn, Z1_node, Z2_node)
                 d_ub  <- marginal.cond.quantile_z1(h, 0.975, Z2_obs[i], i, data_test_dnn, Z1_node, Z2_node)
                 
                 c(g_med, g_lb, g_ub, d_med, d_lb, d_ub)
               }

stopCluster(cl); registerDoSEQ()

median_Gaussian <- res[, 1]; lb_Gaussian <- res[, 2]; ub_Gaussian <- res[, 3]
median_dnn      <- res[, 4]; lb_dnn      <- res[, 5]; ub_dnn      <- res[, 6]


mad_Gaussian <- mean(abs(data_test$Z1 - median_Gaussian))
mad_dnn <- mean(abs(data_test$Z1 - median_dnn))

picp_Gaussian <- mean(data_test$Z1 > lb_Gaussian & data_test$Z1 < ub_Gaussian)
picp_dnn <- mean(data_test$Z1 > lb_dnn & data_test$Z1 < ub_dnn)

length_Gaussian <- mean(ub_Gaussian - lb_Gaussian)
length_dnn <- mean(ub_dnn - lb_dnn)

cat("Gaussian results:\n")
cat("  MAD   =", mad_Gaussian, "\n")
cat("  PICP  =", picp_Gaussian, "\n")
cat("  Length=", length_Gaussian, "\n\n")

cat("DNN results:\n")
cat("  MAD   =", mad_dnn, "\n")
cat("  PICP  =", picp_dnn, "\n")
cat("  Length=", length_dnn, "\n")

