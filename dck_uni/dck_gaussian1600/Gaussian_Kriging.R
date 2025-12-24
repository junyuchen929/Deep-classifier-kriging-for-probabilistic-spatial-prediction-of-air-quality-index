library(geoR)
library(doParallel)
library(foreach)

rm(list = ls())

g_val <- 0
h_val <- 0
num_sim <- 100

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

time_list_Gaussian <- vector("numeric", num_sim)

pred_list_Gaussian <- foreach(sim = 1:num_sim, .packages = c("geoR"), .combine = "c") %dopar% {
  
  start_time <- Sys.time()  # start
  
  file_path <- paste0("synthetic_data_simulations/training_data/Tgh_nonGaussian_1600_classification_g", g_val,
                      "_h", h_val, "_", sim, "train.csv")
  data_train <- read.csv(file_path, header = TRUE)
  
  file_path <- paste0("synthetic_data_simulations/testing_data/Tgh_nonGaussian_1600_classification_g", g_val,
                      "_h", h_val, "_", sim, "test.csv")
  data_test <- read.csv(file_path, header = TRUE)
  
  geodata_train <- as.geodata(data_train, coords.col = c("x", "y"), data.col = "tgh")
  train_loc <- data_train[,1:2]
  test_loc <- data_test[,1:2]
  
  fit.Model <- likfit(geodata_train,
                      ini.cov.pars = c(4, 2), nugget = 0.1,
                      cov.model = "exponential",
                      fix.nugget = FALSE)
  
  pred <- krige.conv(geodata_train, locations = test_loc,
                     krige = krige.control(cov.model = "exponential", type.krige = "SK",
                                           beta = mean(data_train$tgh),
                                           cov.pars = c(fit.Model$sigmasq, fit.Model$phi),
                                           nugget = fit.Model$nugget))
  
  end_time <- Sys.time()  # end
  comp_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  list(list(pred_value = pred$predict,
            conditional_var = pred$krige.var,
            computation_time = comp_time))
}

stopCluster(cl)

# computation times
time_list_Gaussian <- sapply(pred_list_Gaussian, function(x) x$computation_time)

# pred results
pred_list_clean <- lapply(pred_list_Gaussian, function(x) {
  x$computation_time <- NULL
  return(x)
})
saveRDS(pred_list_clean, paste0("pred_list_Gaussian_exp_g", g_val, "_h", h_val, ".rds"))
write.csv(time_list_Gaussian, paste0("GK_time_list_g", g_val, "_h", h_val, ".csv"), row.names = FALSE)

# avg time
avg_time <- mean(time_list_Gaussian)
cat("Average computation time per simulation (in seconds):", avg_time, "\n")
