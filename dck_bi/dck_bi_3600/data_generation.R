rm(list = ls())
library(geoR)
library(MASS)
library(fields)
setwd("/Users/beauniverse/Documents/GitHub/Deep-classifier-kriging-for-probabilistic-spatial-prediction-of-air-quality-index/dck_bi/dck_bi_3600")
mainDir <- "."
subDir <- "synthetic_data_simulations/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

set.seed(12345567)
num_sim = 100

x <- seq(0,1, length.out = 60)
y <- seq(0,1, length.out = 60)

d1 <- expand.grid(x = x, y = y)
X = d1$x             
Y = d1$y
m <- as.matrix(dist(data.frame(X=X,Y=Y)))

R = 0.8
del1 = 0.7
del2 = 0.9
R_corr1 = 0.3
R_corr2 = 0.88

s11 = 0.89
s22 = 1.3

nu11 = 0.8 
nu22 = 0.8 
nu12 = (nu11 + nu22) /2 #+ (del1*(1-R_corr1))

alpha11 = 0.2
alpha22 = 0.4
alpha12 = (alpha11 + alpha22)/2 #+ (del2*(1-R_corr2))


s12 = sqrt(s11*s22)*(alpha11^(nu11/2))*(alpha22^(nu22/2))*gamma(nu12)*R / 
  ((alpha12^nu12)*sqrt(gamma(nu11)*gamma(nu22)))

constant = (sqrt(s11*s22)*R*(gamma(nu12)/sqrt(gamma(nu11)*gamma(nu22)))) /
  ((alpha12^nu12)/sqrt(alpha11^nu11 * alpha22^nu22))

matern_cov1 = s11*matern(m,sqrt(alpha11),nu11)
matern_cov2 <- constant*matern(m,sqrt(alpha12),nu12)
matern_cov4 = s22*matern(m,sqrt(alpha22),nu22)
full_matern_cov = rbind(cbind(matern_cov1,matern_cov2),cbind(t(matern_cov2),matern_cov4))

for(i in 1:num_sim){
  simulation = mvrnorm(1,rep(0,2*3600),full_matern_cov)
  
  var1 = simulation[1:3600]
  var2 = simulation[3601:(2*3600)]
  g = 0
  h = 0
  
  if (g == 0) {
    tukey_var1 = var1 * exp(h * (var1^2) / 2)
    tukey_var2 = var2 * exp(h * (var2^2) / 2)
  } else {
    tukey_var1 = ((exp(g * var1) - 1) / g) * exp(h * (var1^2) / 2)
    tukey_var2 = ((exp(g * var2) - 1) / g) * exp(h * (var2^2) / 2)
  }
  
  

  
  df = data.frame(x=X,y=Y,var1 = tukey_var1, var2 = tukey_var2)
  write.csv(x = df,
            file = paste0("synthetic_data_simulations/2d_gaussian_3600_",i,".csv"),
            row.names = F)
}

mainDir <- "./synthetic_data_simulations/"
subDir <- "training_data/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
subDir <- "testing_data/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

mainDir <- "."
subDir <- "Results_DNN/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)


