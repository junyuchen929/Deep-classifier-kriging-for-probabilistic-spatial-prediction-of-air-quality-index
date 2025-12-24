rm(list = ls())

source("LMC_functions.R")
source("density_computation_DNN.R")
library(VGAM)
library(foreach)
library(doParallel)

num_sim <- 100
C_grid <- c(12)
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

pred_list_all <- readRDS("pred_list_Gaussian.rds")
results_all <- list()

for (C in C_grid) {
  results <- foreach(sim_num = 1:num_sim, .combine = rbind,
                     .packages = c("VGAM")) %dopar% {
                       
                       file_path = paste0("synthetic_data_simulations/testing_data/2D_Gaussian_3600_projection_", sim_num, "test.csv")
                       data_test = read.csv(file_path, header = TRUE)
                       
                       file_path = paste0("Results_DNN/2D_Gaussian_3600_predictions_", sim_num, ".csv")
                       data_test_dnn = read.csv(file_path, header = TRUE)
                       
                       file_path = paste0("synthetic_data_simulations/training_data/2D_Gaussian_3600_projection_", sim_num, "train.csv")
                       data_train = read.csv(file_path, header = TRUE)
                       
                       sigma <- mad(data_train$var1, constant = 1)
                       n <- nrow(data_train)
                       alpha <- 1/3
                       h <- 2 * C * sigma * n^(-alpha) / 3
                       
                       total_classes <- ncol(data_test_dnn) - 11
                       Z1_node <- matrix(NA, nrow = 1, ncol = total_classes)
                       Z2_node <- matrix(NA, nrow = 1, ncol = total_classes)
                       for (i in 1:total_classes) {
                         df_sample <- data_train[data_train$class == i, ]
                         Z1_node[, i] <- mean(df_sample$var1)
                         Z2_node[, i] <- mean(df_sample$var2)
                       }
                       
                       pred_list_Gaussian <- lapply(pred_list_all, function(x) x$pred_value)
                       pred_obj <- pred_list_Gaussian[[sim_num]]
                       
                       n_test <- 360
                       median_Gaussian <- lb_Gaussian <- ub_Gaussian <- 
                         median_dnn <- lb_dnn <- ub_dnn <- rep(NA, n_test)
                       
                       for (index_loc in 1:n_test) {
                         preds <- pred_obj$pred[c(index_loc, n_test + index_loc), 1]
                         cond_var <- pred_obj$conditional_var[c(index_loc, n_test + index_loc), c(index_loc, n_test + index_loc)]
                         ro <- cond_var[1, 2] / sqrt(cond_var[1, 1] * cond_var[2, 2])
                         varx <- cond_var[1, 1]  # sigma1^2
                         vary <- cond_var[2, 2]  # sigma2^2
                         cond_pred <- preds[1] + ro * sqrt(varx / vary) * (data_test$var2[index_loc] - preds[2])
                         cond_var_scalar <- (1 - ro^2) * cond_var[1, 1]
                         
                         median_Gaussian[index_loc] <- qnorm(0.5, cond_pred, sqrt(cond_var_scalar))
                         lb_Gaussian[index_loc] <- qnorm(0.025, cond_pred, sqrt(cond_var_scalar))
                         ub_Gaussian[index_loc] <- qnorm(0.975, cond_pred, sqrt(cond_var_scalar))
                         
                         median_dnn[index_loc] <- marginal.cond.quantile_z1(h, 0.5, data_test$var2[index_loc], index_loc, data_test_dnn, Z1_node, Z2_node)
                         lb_dnn[index_loc] <- marginal.cond.quantile_z1(h, 0.025, data_test$var2[index_loc], index_loc, data_test_dnn, Z1_node, Z2_node)
                         ub_dnn[index_loc] <- marginal.cond.quantile_z1(h, 0.975, data_test$var2[index_loc], index_loc, data_test_dnn, Z1_node, Z2_node)
                       }
                       
                       mad_Gaussian <- mean(abs(data_test$var1 - median_Gaussian))
                       mad_dnn <- mean(abs(data_test$var1 - median_dnn))
                       
                       picp_Gaussian <- mean(data_test$var1 > lb_Gaussian & data_test$var1 < ub_Gaussian)
                       picp_dnn <- mean(data_test$var1 > lb_dnn & data_test$var1 < ub_dnn)
                       
                       length_Gaussian <- mean(ub_Gaussian - lb_Gaussian)
                       length_dnn <- mean(ub_dnn - lb_dnn)
                       
                       c(sim_num, mad_Gaussian, mad_dnn, picp_Gaussian, picp_dnn, length_Gaussian, length_dnn, h)
                     }
  
  results_all[[as.character(C)]] <- results
}


for (C_val in C_grid) {
  results <- results_all[[as.character(C_val)]]
  file_name <- paste0("C=", C_val, ".csv")
  write.csv(results, file_name, row.names = FALSE)
}


results <- read.csv("C=12.csv", header = TRUE)
colnames(results) <- c("sim_num", "mad_Gaussian", "mad_dnn",
                       "picp_Gaussian", "picp_dnn",
                       "length_Gaussian", "length_dnn", "h")
print(round(mean(results$mad_Gaussian), 3))
print(round(sd(results$mad_Gaussian)/10, 3))
print(round(mean(results$picp_Gaussian), 3))
print(round(sd(results$picp_Gaussian)/10, 3))
print(round(mean(results$length_Gaussian), 3))
print(round(sd(results$length_Gaussian)/10, 3))
print(round(mean(results$mad_dnn), 3))
print(round(sd(results$mad_dnn)/10, 3))
print(round(mean(results$picp_dnn), 3))
print(round(sd(results$picp_dnn)/10, 3))
print(round(mean(results$length_dnn), 3))
print(round(sd(results$length_dnn)/10, 3))
print(round(mean(results$h), 3))
