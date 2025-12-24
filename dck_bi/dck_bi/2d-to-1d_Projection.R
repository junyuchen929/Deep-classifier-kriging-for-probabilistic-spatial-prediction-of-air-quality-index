rm(list = ls())
library(ddalpha)
library(ggplot2)
library(RColorBrewer)
library(quantreg)
library(calculus)

# 2d-nonGaussian data generated from tukey g-h transformation over parsimonious matern
num_sim = 100
min_samples_per_class = 15
for(sim in 1:num_sim){
  file_path <- paste0("synthetic_data_simulations/2d_gaussian_1600_",sim,".csv")
  data <- read.csv(file_path, sep = ",", header = T)
  
  n_model = 5
  # n_splits_per_model = 5
  
  tau = seq(0.05,0.95,length.out = n_model)
  a = rep(NA,n_model)
  b = rep(NA,n_model)
  for(i in 1:n_model){
    fit <- rq(data$var1 ~ data$var2, tau = tau[i])
    a[i] = unname(fit$coefficients[1])
    b[i] = unname(fit$coefficients[2])
  }

  
  df_threshold_list = list()
  threshold = rep(NA,dim(data)[1])
  proj_var1 = rep(NA,dim(data)[1])
  proj_var2 = rep(NA,dim(data)[1])
  for(i in 1:dim(data)[1]){
    
    vec_distance = abs(data$var1[i] - (a + b*data$var2[i])) / sqrt(1+b^2)
    index = which.min(vec_distance)
    
    proj.var2 = (data$var2[i] + b[index]*(data$var1[i] - a[index])) / (1 + b[index]^2)
    proj.var1 = a[index] + b[index]*proj.var2
    
    proj_var1[i] = proj.var1
    proj_var2[i] = proj.var2
    threshold[i] = index
  }
  
  data$proj_var1 = proj_var1
  data$proj_var2 = proj_var2
  data$threshold = threshold
  rm(threshold,proj_var1,proj_var2)
  dist_frm_org = rep(NA,dim(data)[1])
  
  for(i in 1:dim(data)[1]){
    dist_frm_org[i] = (data$proj_var1[i] - a[data$threshold[i]])^2 + (data$proj_var2[i]^2)
  }
  
  data$dist_frm_orig = dist_frm_org
  
  data_list = list()
  tot_class_number = 0
  for(i in 1:n_model){
    df = data[data$threshold == i,]
    # threshold2 = rep(NA,dim(df)[1])
    min_dist = min(df$dist_frm_orig)
    max_dist = max(df$dist_frm_orig)
    df = df[order(df$dist_frm_orig),]
    # row.names(df) <- NULL
    if(dim(df)[1] >= min_samples_per_class){
      if(dim(df)[1]%%min_samples_per_class >= 10){
        n_class = dim(df)[1]%/%min_samples_per_class + 1
      }else{
        n_class = dim(df)[1]%/%min_samples_per_class
      }
      }
    else{n_class = 1}

    batch = c()
    for(class_number in 1:(n_class-1)){
      tot_class_number = tot_class_number + 1
      batch = c(batch,rep((tot_class_number),min_samples_per_class))
    }
    li = length(batch)
    tot_class_number = tot_class_number + 1
    threshold2 = c(batch,rep((tot_class_number),(dim(df)[1] - li)))
    
    df$threshold2 = threshold2
    # df$constant = rep(((i-1)*n_class),dim(df)[1])
    data_list[[i]] = df
    print(i)
  }
  
  df_threshold2 = do.call(rbind, data_list)
  df_threshold2 = na.omit(df_threshold2)
  
  df_threshold2$class = df_threshold2$threshold2
  df_threshold2$index <- as.numeric(row.names(df_threshold2))
  df_threshold2 = df_threshold2[order(df_threshold2$index), ]
  total_classes = length(unique(df_threshold2$class))
  df = df_threshold2
  Z1_node = matrix(NA, nrow = 1, ncol = total_classes )
  Z2_node = matrix(NA, nrow = 1, ncol = total_classes )
  
  for(i in 1:total_classes){
    df_sample = df[df$class == i,]
    Z1_node[, i] <- mean(df_sample$proj_var1)  
    Z2_node[, i] <- mean(df_sample$proj_var2)  
  }
  write.csv(x = df_threshold2,
            file = paste0("synthetic_data_simulations/2D_Gaussian_1600_projection_",sim,".csv"),
            row.names = F)
   
}  

