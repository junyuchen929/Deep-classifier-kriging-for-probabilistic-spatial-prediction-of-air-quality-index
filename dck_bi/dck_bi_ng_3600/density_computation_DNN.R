### bivariate conditional cdfs
library(VGAM)
cond.cdf_2d <- function(h, z, index, data, Z1_node, Z2_node){
  sum = 0
  z1 = z[1]
  z2 = z[2]
  
  for(i in 12:dim(data)[2]){
    z1_res = z1 - Z1_node[,(i-11)]
    z2_res = z2 - Z2_node[,(i-11)]
    p = data[index,i]
    sum = sum + (p*pbinorm(z1_res,z2_res,var1 = h,var2 = h))
    # sum = sum + p*pnorm(z1_res/h)*pnorm(z2_res/h)
  }
  return(sum)
}

marginal.z2 <- function(h, z2, index, data, Z2_node){
  sum = 0
  for(i in 12:dim(data)[2]){
    z2_res = (z2 - Z2_node[,(i-11)])
    p = data[index,i]
    sum = sum + (p*dnorm(z2_res/h))
  }
  return(sum)
}

marginal.cond.cdf_z1 <- function(h, z, index, data, Z1_node, Z2_node){
  sum = 0
  z1 = z[1]
  z2 = z[2]
  for(i in 12:dim(data)[2]){
    z1_res = (z1 - Z1_node[,(i-11)])
    z2_res = (z2 - Z2_node[,(i-11)])
    p = data[index,i]
    sum = sum + (p*pnorm(z1_res/h)*dnorm(z2_res/h))
  }
  sum = sum/marginal.z2(h,z2,index,data,Z2_node)
  return(sum)
}

marginal.cond.quantile_z1 <- function(h,p,z2_given,index,data,Z1_node,Z2_node){
  z1_val = seq(-25,25,length.out = 250)
  z1_cond.cdf = rep(NA,250)
  for(i in 1:250){
    z = c(z1_val[i],z2_given)
    z1_cond.cdf[i] = marginal.cond.cdf_z1(h,z,index,data, Z1_node, Z2_node)
  }
  # return(approx(z1_cond.cdf,z1_val,p)$y)
  
  valid_cdf <- z1_cond.cdf[!is.na(z1_cond.cdf)]
  if (p < min(valid_cdf)) {
    return(-25)
  } else if (p > max(valid_cdf)) {
    return(25)
  } else {
    return(approx(z1_cond.cdf, z1_val, p, rule = 2)$y)
  }
}
  
