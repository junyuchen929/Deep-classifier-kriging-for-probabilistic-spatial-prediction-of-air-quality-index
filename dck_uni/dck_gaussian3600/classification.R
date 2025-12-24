rm(list = ls())
setwd("/Users/beauniverse/Documents/GitHub/Deep-classifier-kriging-for-probabilistic-spatial-prediction-of-air-quality-index/dck_uni/dck_gaussian3600")

g_val <- 0
h_val <- 0

num_sim <- 100
n_class <- 30  # change this later

tau <- seq(0.01, 0.99, length.out = n_class)
for(sim in 1:num_sim){

  file_path <- paste0("synthetic_data_simulations/Tgh_sim_g", g_val,
                      "_h", h_val, "_", sim, ".csv")
  data <- read.csv(file_path, sep = ",", header = TRUE)

  quantiles <- quantile(data$tgh, probs = tau)
  
  data$threshold <- cut(data$tgh, 
                        breaks = c(min(data$tgh)-1, quantiles, max(data$tgh)), 
                        labels = 1:(length(quantiles) + 1), 
                        right = TRUE)

  data$index <- as.numeric(row.names(data))

  out_file <- paste0("synthetic_data_simulations/Tgh_Gaussian_3600_classification_g", g_val,
                     "_h", h_val, "_", sim, ".csv")
  write.csv(x = data, file = out_file, row.names = FALSE)
}
