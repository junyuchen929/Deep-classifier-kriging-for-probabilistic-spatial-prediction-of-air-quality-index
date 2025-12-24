library(geoR)
library(doParallel)
library(foreach)
library(fields)

rm(list = ls())

g_val <- 0.8
h_val <- 0.5
num_sim <- 100

num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)

time_list_Gaussian <- vector("numeric", num_sim)
pred_list_Gaussian <- foreach(sim = 1:num_sim, .packages = c("fields")) %dopar% {
  start_time <- Sys.time()
  
  file_path = paste0("synthetic_data_simulations/training_data/sim1600_classified_g", g_val,
                     "_h", h_val, "_", sim, "train.csv")
  data_train = read.csv(file_path, header = TRUE)
  file_path = paste0("synthetic_data_simulations/testing_data/sim1600_classified_g", g_val,
                     "_h", h_val, "_", sim, "test.csv")
  data_test = read.csv(file_path, header = TRUE)
  
  Z_train <- as.matrix(data_train[, c("X1", "X2", "X3", "X4", "X5")])
  Z_test <- as.matrix(data_test[, c("X1", "X2", "X3", "X4", "X5")])
  
  krig_model <- Krig(
    x = as.matrix(data_train[, c("x", "y")]),
    Y = data_train$var1,
    Z = Z_train,
    cov.function = "stationary.cov",
    cov.args = list(Covariance = "Matern", nu = 0.3) 
  )
  
  predicted_values <- predict(krig_model,
                              x = as.matrix(data_test[, c("x", "y")]),
                              Z = Z_test)
  predicted_se <- predictSE(krig_model,
                            x = as.matrix(data_test[, c("x", "y")]),
                            Z = Z_test)
  
  end_time <- Sys.time()
  comp_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  list(pred_value = predicted_values, pred_SE = predicted_se, computation_time = comp_time)
}

time_list_Gaussian <- sapply(pred_list_Gaussian, function(x) x$computation_time)
pred_list_clean <- lapply(pred_list_Gaussian, function(x) {
  x$computation_time <- NULL
  return(x)
})
saveRDS(pred_list_clean, paste0("pred_list_Gaussian_exp_g", g_val, "_h", h_val,".rds"))
write.csv(time_list_Gaussian, paste0("Krig_time_list_g", g_val, "_h", h_val, ".csv"), row.names = FALSE)
cat("Average computation time per simulation (Original):", mean(time_list_Gaussian), "seconds\n")



time_list_Gaussian_Stand <- vector("numeric", num_sim)
pred_list_Gaussian_Stand <- foreach(sim = 1:num_sim, .packages = c("fields")) %dopar% {
  start_time <- Sys.time()
  
  file_path_train <- paste0("synthetic_data_simulations/training_data/sim1600_classified_g",
                            g_val, "_h", h_val, "_", sim, "train_Stand.csv")
  data_train <- read.csv(file_path_train, header = TRUE)
  
  file_path_test <- paste0("synthetic_data_simulations/testing_data/sim1600_classified_g",
                           g_val, "_h", h_val, "_", sim, "test_Stand.csv")
  data_test <- read.csv(file_path_test, header = TRUE)
  
  Z_train <- as.matrix(data_train[, c("X1", "X2", "X3", "X4", "X5")])
  Z_test <- as.matrix(data_test[, c("X1", "X2", "X3", "X4", "X5")])
  
  krig_model <- Krig(
    x = as.matrix(data_train[, c("x", "y")]),
    Y = data_train$var1,
    Z = Z_train,
    cov.function = "stationary.cov",
    cov.args = list(Covariance = "Matern", nu = 0.3) 
  )
  
  predicted_values <- predict(krig_model,
                              x = as.matrix(data_test[, c("x", "y")]),
                              Z = Z_test)
  predicted_se <- predictSE(krig_model,
                            x = as.matrix(data_test[, c("x", "y")]),
                            Z = Z_test)
  
  end_time <- Sys.time()
  comp_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  list(pred_value = predicted_values, pred_SE = predicted_se, computation_time = comp_time)
}

stopCluster(cl)

time_list_Gaussian_Stand <- sapply(pred_list_Gaussian_Stand, function(x) x$computation_time)
pred_list_clean_Stand <- lapply(pred_list_Gaussian_Stand, function(x) {
  x$computation_time <- NULL
  return(x)
})
saveRDS(pred_list_clean_Stand, paste0("pred_list_Gaussian_exp_g", g_val, "_h", h_val, "_Stand.rds"))
write.csv(time_list_Gaussian_Stand, paste0("Krig_time_list_g", g_val, "_h", h_val, "_Stand.csv"), row.names = FALSE)
cat("Average computation time per simulation (Standardized):", mean(time_list_Gaussian_Stand), "seconds\n")
