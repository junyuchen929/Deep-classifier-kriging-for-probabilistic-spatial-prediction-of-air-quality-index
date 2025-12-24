rm(list = ls())
library(doParallel)
library(foreach)
source("LMC_functions.R")

num_sim <- 100
num_cores <- detectCores()-1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

pred_list_Gaussian <- foreach(sim = 1:num_sim, .packages = c("fields", "MASS", "geoR"), .combine = "c") %dopar% {
  start_time <- Sys.time()

  file_path <- paste0("synthetic_data_simulations/training_data/2D_nonGaussian_1600_projection_", sim, "train.csv")
  data_train <- read.csv(file_path, header = TRUE)
  
  file_path <- paste0("synthetic_data_simulations/testing_data/2D_nonGaussian_1600_projection_", sim, "test.csv")
  data_test <- read.csv(file_path, header = TRUE)
  
  init.ind <- c(1, 1, 1, 1, 1, 1, 0, 0)
  indmat.estim <- optim_indmat_loglik(init.ind)
  
  init.lmc <- c(indmat.estim$par[1:6], 0.1,0.1,0.1,0.1, indmat.estim$par[7:8])
  fit.Model.lmc <- optim_lmc.loglikelihood(par = init.lmc)
  
  pred <- lmc.pred.summary(fit.Model.lmc$par, data_train, data_test)
  
  end_time <- Sys.time()
  comp_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  list(list(pred_value = pred, computation_time = comp_time))
}

stopCluster(cl)

time_list_Gaussian <- sapply(pred_list_Gaussian, function(x) x$computation_time)

pred_list_clean <- lapply(pred_list_Gaussian, function(x) {
  x$computation_time <- NULL
  return(x)
})

saveRDS(pred_list_clean, "pred_list_nonGaussian.rds")
write.csv(time_list_Gaussian, "LMC_time_list.csv", row.names = FALSE)

cat("Average computation time per simulation (in seconds):", mean(time_list_Gaussian), "\n")
