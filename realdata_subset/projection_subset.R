rm(list = ls())
setwd("/Users/beauniverse/Documents/GitHub/Deep-classifier-kriging-for-probabilistic-spatial-prediction-of-air-quality-index/realdata_subset")


library(ddalpha)
library(ggplot2)
library(RColorBrewer)
library(quantreg)
library(calculus)
library(foreach)
library(doParallel)

# -------------------------
# Parameters
# -------------------------
n_model = 5
tau = seq(0.05, 0.95, length.out = n_model)
min_samples_per_class = 55

alpha  = 1.0
sigma0 = 1e-6
gamma  = 1.0
min_weight_part2 = 0.3

regions = c("CA", "NE", "SE")

for(reg in regions){
  
  cat("\n====================\n")
  cat("Processing Region:", reg, "\n")
  cat("====================\n")
  
  # Load CSV1: region-specific TRAINING ONLY
  data <- read.csv(paste0("CSV1_", reg, ".csv"), header = TRUE)
  if(!all(c("uj_lat","uj_lon","Z1","Z2") %in% names(data))){
    stop("CSV1 must contain columns: uj_lat, uj_lon, Z1, Z2")
  }
  
  # Fit quantile regressions Z1 ~ Z2
  a = numeric(n_model)
  b = numeric(n_model)
  for(j in 1:n_model){
    fit <- rq(Z1 ~ Z2, data = data, tau = tau[j])
    a[j] <- unname(coef(fit)[1])
    b[j] <- unname(coef(fit)[2])
  }
  
  # Part1 projection
  n1 <- nrow(data)
  proj_Z1 <- numeric(n1)
  proj_Z2 <- numeric(n1)
  threshold <- integer(n1)
  resid_selected <- numeric(n1)
  dist_frm_org <- numeric(n1)
  
  for(i in 1:n1){
    z1 <- data$Z1[i]
    z2 <- data$Z2[i]
    r_vec <- z1 - (a + b * z2)
    jstar <- which.min(abs(r_vec))
    threshold[i] <- jstar
    proj_Z2[i] <- z2
    proj_Z1[i] <- a[jstar] + b[jstar] * z2
    resid_selected[i] <- r_vec[jstar]
    dist_frm_org[i] = (proj_Z1[i] - a[jstar])^2 + (proj_Z2[i]^2)
  }
  
  data$proj_Z1 <- proj_Z1
  data$proj_Z2 <- proj_Z2
  data$threshold <- threshold
  data$resid <- resid_selected
  data$dist_frm_orig <- dist_frm_org
  
  # Robust scale (MAD-based)
  s_med <- median(abs(data$resid), na.rm = TRUE)
  if(s_med <= 0) s_med <- sd(data$resid, na.rm = TRUE)
  if(is.na(s_med) || s_med == 0) s_med <- 1e-6
  data$weight <- exp(-alpha * abs(data$resid) / s_med)
  data$sigma2 <- sigma0 + gamma * (data$resid^2)
  
  # Load CSV2: region-specific CMAQ ONLY
  data2 <- read.csv(paste0("CSV2_", reg, ".csv"), header = TRUE)
  if(!all(c("uj_lat","uj_lon","Z2") %in% names(data2))){
    stop("CSV2 must contain columns: uj_lat, uj_lon, Z2")
  }
  
  # Soft projection
  data2_proj <- do.call(rbind, lapply(1:nrow(data2), function(i){
    z2 <- data2$Z2[i]
    dx <- data$uj_lon - data2$uj_lon[i]
    dy <- data$uj_lat - data2$uj_lat[i]
    distances <- sqrt(dx^2 + dy^2)
    nn_idx <- order(distances)[1:5]
    qvals <- a + b * z2
    
    resid_mat <- sapply(1:n_model, function(j){
      data$Z1[nn_idx] - (a[j] + b[j] * data$Z2[nn_idx])
    })
    abs_resid_mat <- abs(resid_mat)
    min_position <- which(abs_resid_mat == min(abs_resid_mat),
                          arr.ind = TRUE)
    
    threshold_idx <- min_position[2]
    center_idx <- threshold_idx
    r_from_center <- abs(qvals - qvals[center_idx])
    
    p <- exp(-(r_from_center /
                 (median(r_from_center) + 1e-6)))
    p <- p / sum(p)
    
    z1_expected <- sum(p * qvals)
    
    data.frame(
      uj_lat = data2$uj_lat[i],
      uj_lon = data2$uj_lon[i],
      Z1 = z1_expected,
      Z2 = z2,
      proj_Z1 = z1_expected,
      proj_Z2 = z2,
      threshold = center_idx,
      dist_frm_orig = (z1_expected - a[center_idx])^2 + (z2^2)
    )
  }))
  
  # Combine & group into final classes
  combined_data <- rbind(
    data[,c("uj_lat","uj_lon","Z1","Z2",
            "proj_Z1","proj_Z2","threshold","dist_frm_orig")],
    data2_proj
  )
  
  missing_var1_idx <- is.na(combined_data$Z1)
  combined_data$Z1[missing_var1_idx] <-
    combined_data$proj_Z1[missing_var1_idx]
  
  data_list <- list()
  tot_class_number <- 0
  
  for(i in 1:n_model){
    df <- combined_data[combined_data$threshold==i,]
    if(nrow(df)==0) next
    df <- df[order(df$dist_frm_orig),]
    
    if(nrow(df)>=min_samples_per_class){
      if(nrow(df) %% min_samples_per_class >= 20){
        n_class <- nrow(df)%/%min_samples_per_class + 1
      } else {
        n_class <- nrow(df)%/%min_samples_per_class
      }
    } else n_class <- 1
    
    batch <- c()
    for(class_number in 1:(n_class - 1)){
      tot_class_number <- tot_class_number + 1
      batch <- c(batch,
                 rep(tot_class_number, min_samples_per_class))
    }
    li <- length(batch)
    tot_class_number <- tot_class_number + 1
    threshold2 <- c(
      batch,
      rep((tot_class_number),(nrow(df)-li))
    )
    
    df$threshold2 <- threshold2
    data_list[[length(data_list)+1]] <- df
    cat(" Quantile", i, "→ Classes:", n_class,
        "Rows:", nrow(df), "\n")
  }
  
  df_threshold2 <- do.call(rbind, data_list)
  df_threshold2 <- na.omit(df_threshold2)
  df_threshold2$class <- df_threshold2$threshold2
  
  output_file <- paste0("projection_matrix_", reg, ".csv")
  write.csv(df_threshold2, file=output_file, row.names=FALSE)
  
  cat("\n🎯 Saved:", output_file, "\n")
}

cat("\n===============================\n")
cat(" 🎯 All Regions Processed! 🎯\n")
cat("===============================\n")







