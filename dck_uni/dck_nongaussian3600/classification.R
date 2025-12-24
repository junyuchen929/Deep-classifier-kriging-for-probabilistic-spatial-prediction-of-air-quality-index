rm(list = ls())

# choose different g, h here
g <- 0.8
h <- 0.5

num_sim <- 100
n_class <- 30  # change this later

tau <- seq(0.01, 0.99, length.out = n_class)


for(sim in 1:num_sim){

  file_path <- paste0("synthetic_data_simulations/sim3600_g", g, "_h", h, "_", sim, ".csv")
  data <- read.csv(file_path, sep = ",", header = TRUE)

  quantiles <- quantile(data$var1, probs = tau)
  
  data$threshold <- cut(data$var1, 
                        breaks = c(min(data$var1)-1, quantiles, max(data$var1)), 
                        labels = 1:(length(quantiles) + 1), 
                        right = TRUE)

  data$index <- as.numeric(row.names(data))

  out_file <- paste0("synthetic_data_simulations/sim3600_classified_g", g, "_h", h, "_", sim, ".csv")
  write.csv(x = data, file = out_file, row.names = FALSE)
}

