rm(list = ls())
setwd("/Users/beauniverse/Documents/GitHub/Deep-classifier-kriging-for-probabilistic-spatial-prediction-of-air-quality-index/realdata_whole")

library(ddalpha)
library(ggplot2)
library(RColorBrewer)
library(quantreg)
library(calculus)
library(foreach)
library(doParallel)

# Load CSV1.csv for quantile regression
data <- read.csv("CSV1.csv", sep = ",", header = T)

# Perform quantile regression using Z1 ~ Z2
n_model = 5
tau = seq(0.05,0.95,length.out = n_model)
a = rep(NA,n_model)
b = rep(NA,n_model)
for(i in 1:n_model){
  fit <- rq(data$Z1 ~ data$Z2, tau = tau[i])
  a[i] = unname(fit$coefficients[1])
  b[i] = unname(fit$coefficients[2])
}

# Load CSV2.csv
data2 <- read.csv("CSV2_10000.csv", sep = ",", header = T)
#data2 <- read.csv("CSV2.csv", sep = ",", header = T)

# Part 1: projection of data points in CSV1 onto quantile regression lines
df_threshold_list = list()
threshold = rep(NA,dim(data)[1])
proj_Z1 = rep(NA,dim(data)[1])
proj_Z2 = rep(NA,dim(data)[1])
for(i in 1:dim(data)[1]){
  
  vec_distance = abs(data$Z1[i] - b*data$Z2[i] - a)/sqrt(1+b^2)
  index = order(vec_distance)[1]
  
  proj.Z2 = (b[index]*data$Z1[i] + data$Z2[i] - a[index]*b[index])/(1 + b[index]^2)
  proj.Z1 = b[index]*proj.Z2 + a[index]
  
  proj_Z1[i] = proj.Z1
  proj_Z2[i] = proj.Z2
  threshold[i] = index
}

data$proj_Z1 = proj_Z1
data$proj_Z2 = proj_Z2
data$threshold = threshold
rm(threshold,proj_Z1,proj_Z2)

# Compute distance from original points to the projection
dist_frm_org = rep(NA,dim(data)[1])
for(i in 1:dim(data)[1]){
  dist_frm_org[i] = (data$proj_Z1[i] - a[data$threshold[i]])^2 + (data$proj_Z2[i]^2)
}
data$dist_frm_orig = dist_frm_org

# Part 2: Project CSV2 points onto all 5 quantile regression lines
# For each point in CSV2, compute 5 predicted Z1 values (1 per quantile line)
data2_expanded <- data.frame()

# Initialize parallel backend
cores <- detectCores() - 1
cl <- makeCluster(cores)
registerDoParallel(cl)

# Perform parallel projection
data2_expanded <- foreach(i = 1:nrow(data2), .combine = rbind, .packages = "base") %dopar% {
  rows <- vector("list", n_model)
  for(j in 1:n_model){
    z1_predicted <- b[j] * data2$Z2[i] + a[j]
    rows[[j]] <- data.frame(
      uj_lat = data2$uj_lat[i],
      uj_lon = data2$uj_lon[i],
      Z1 = z1_predicted,
      Z2 = data2$Z2[i],
      proj_Z1 = z1_predicted,
      proj_Z2 = data2$Z2[i],
      threshold = j,
      dist_frm_orig = (z1_predicted - a[j])^2 + (data2$Z2[i]^2)
    )
  }
  do.call(rbind, rows)
}

# Stop parallel backend
stopCluster(cl)

# Combine CSV1-projected data with CSV2-expanded projection
combined_data <- rbind(
  data[, c("uj_lat", "uj_lon", "Z1", "Z2", "proj_Z1", "proj_Z2", "threshold", "dist_frm_orig")],
  data2_expanded
)

# Group data into subclasses based on distance within each quantile regression line
data_list = list()
tot_class_number = 0
min_samples_per_class = 15

for(i in 1:n_model){
  df = combined_data[combined_data$threshold == i,]
  min_dist = min(df$dist_frm_orig)
  max_dist = max(df$dist_frm_orig)
  df = df[order(df$dist_frm_orig),]
  
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
  data_list[[i]] = df
  print(i)
}

# Merge all grouped data
df_threshold2 = do.call(rbind, data_list)
df_threshold2 = na.omit(df_threshold2)

# Rename class variable
df_threshold2$class = df_threshold2$threshold2
df_threshold2$index <- as.numeric(row.names(df_threshold2))
df_threshold2 = df_threshold2[order(df_threshold2$index), ]

# Save final projection matrix
write.csv(x = df_threshold2, file = "projection_matrix_total.csv",row.names = F)
