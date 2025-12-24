rm(list = ls())
library(VGAM)
library(ggplot2)
library(gridExtra)
library(foreach)
library(doParallel)

g_val <- 0.8
h_val <- 0.5

sim_list <- 1:100
m <- 160
total_classes <- 31

cond.cdf <- function(y, index, data, Z1_node, h) {
  CDF = 0
  for (j in 11:dim(data)[2]) {
    CDF = CDF + data[index, j] * pnorm((y - Z1_node[j - 10])/ h) 
  }
  return(CDF)
}

cdf.quantile <- function(p, index, data, Z1_node, h) {
  y = seq(-200, 200, length.out = 600) 
  
  cdf_vals = rep(NA, 600)
  
  for (i in 1:600) {
    cdf_vals[i] = cond.cdf(y[i], index, data, Z1_node, h)
  }
  if (p < min(cdf_vals, na.rm = TRUE)) {
    return(-250)
  } else if (p > max(cdf_vals, na.rm = TRUE)) {
    return(250)
  } else {
    return(approx(cdf_vals, y, p, rule = 2)$y)  
  }
}


pred_list_Gaussian <- readRDS(paste0("pred_list_Gaussian_exp_g", g_val, "_h", h_val,".rds"))

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

C_list <- c(12)
results_allC <- list()
for (C in C_list) {
  results_C <- foreach(sim = sim_list, .packages = c("VGAM","gridExtra"), .combine = "rbind") %dopar% {

    file_path_dnn <- paste0("synthetic_data_simulations/Results_DNN/sim1600_classified_g", 
                            g_val, "_h", h_val, "_", sim, ".csv")
    data_test <- read.csv(file_path_dnn, header = TRUE)
    
    file_path_train <- paste0("synthetic_data_simulations/training_data/sim1600_classified_g", 
                              g_val, "_h", h_val, "_", sim, "train.csv")
    data_train <- read.csv(file_path_train, header = TRUE)
    
    # choose h
    sigma <- mad(data_train$var1, constant = 1)   
    n <- nrow(data_train)   
    alpha = 1/3
    h <- C*sigma*n^(-alpha)/3
    
    tau <- seq(0.01, 0.99, length.out = total_classes-1)
    q <- quantile(data_train$var1, probs = tau)
    
    Z1_node <- rep(NA, total_classes)
    for(i in 1:total_classes) {
      df_sample <- data_train[data_train$threshold == i, ]
      Z1_node[i] <- mean(df_sample$var1)
    }
    
    median_Gaussian <- rep(NA, m)
    lb_Gaussian <- rep(NA, m)
    ub_Gaussian <- rep(NA, m)
    
    median_dnn <- rep(NA, m)
    lb_dnn <- rep(NA, m)
    ub_dnn <- rep(NA, m)
    
    for(index_loc in 1:m) {
      
      # Gaussian-based prediction and intervals
      preds <- pred_list_Gaussian[[sim]]$pred_value[index_loc]
      pred_SE <- pred_list_Gaussian[[sim]]$pred_SE[index_loc]
      
      median_Gaussian[index_loc] <- qnorm(0.5, preds, pred_SE)
      lb_Gaussian[index_loc] <- qnorm(0.025, preds, pred_SE)
      ub_Gaussian[index_loc] <- qnorm(0.975, preds, pred_SE)
      
      # DNN-based prediction and intervals
      median_dnn[index_loc] <- cdf.quantile(0.5, index_loc, data_test, Z1_node, h)
      lb_dnn[index_loc] <- cdf.quantile(0.025, index_loc, data_test, Z1_node, h)
      ub_dnn[index_loc] <- cdf.quantile(0.975, index_loc, data_test, Z1_node, h)
  
    }
    
    ## MAD
    mad_Gaussian <- mean(abs(data_test$var1 - median_Gaussian))
    mad_dnn <- mean(abs(data_test$var1 - median_dnn),na.rm = TRUE)
    
    ## PICP
    count_Gaussian <- sum(data_test$var1 > lb_Gaussian & data_test$var1 < ub_Gaussian)
    count_dnn <- sum(data_test$var1 > lb_dnn & data_test$var1 < ub_dnn,na.rm = TRUE)
   
    picp_Gaussian <- count_Gaussian / m
    picp_dnn <- count_dnn / m
    
    data.frame(
      simulation = sim,
      mad_Gaussian = mad_Gaussian,
      picp_Gaussian = picp_Gaussian,
      CI_Length_Gaussian = mean(ub_Gaussian - lb_Gaussian),
      
      mad_dnn = mad_dnn,
      picp_dnn = picp_dnn,
      CI_Length_dnn = mean(ub_dnn - lb_dnn),
      
      h = h,
      sigma = sigma
    )
  }
  results_allC[[as.character(C)]] <- results_C
}


results = results_allC[[as.character(12)]]
for (C_val in C_list) {
  results <- results_allC[[as.character(C_val)]]
  file_name <- paste0("C=", C_val, ".csv")
  write.csv(results, file_name, row.names = FALSE)
}

stopCluster(cl)
results <- read.csv("C=12.csv", header = TRUE)
mean(results$mad_Gaussian)
sd(results$mad_Gaussian) / 10
mean(results$mad_dnn)
sd(results$mad_dnn) / 10
mean(results$picp_Gaussian)
mean(results$picp_dnn)
sd(results$picp_Gaussian) / 10
sd(results$picp_dnn) / 10
mean(results$CI_Length_Gaussian)
sd(results$CI_Length_Gaussian) / 10
mean(results$CI_Length_dnn)
sd(results$CI_Length_dnn) / 10
mean(results$h)
