library(geoR)
library(MASS)
library(fields)
mle_ind_mat<-function(p,data_train)
{
  un.grd.train <- data_train[,c(1,2,3,4)]
  dist.mat <- rdist(un.grd.train[,-c(3,4)])
  z=c(un.grd.train$var1-mean(un.grd.train$var1) , un.grd.train$var2-mean(un.grd.train$var2))
  a1<-p[1]
  nu1<-p[2]
  sigma1<-p[3]
  a2<-p[4]
  nu2<-p[5]
  sigma2<-p[6]
  nug1<-p[7]
  nug2<-p[8]
  if(sum(p[1:6]<0)!=0||nug1<0||nug2<0)
  {
    nloglikelihood<-10000000
    return(list(mlv=nloglikelihood,params=NULL))
  }
  else
  {
    C11<-sigma1*matern(dist.mat,a1,nu1)
    C22<-sigma2*matern(dist.mat,a2,nu2)
    
    COV12<-matrix(0,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
    
    NUG1<-diag(nug1,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
    NUG2<-diag(nug2,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
    
    ############## Inverting C11 ##########
    C<-rbind(cbind(C11+NUG1,COV12),cbind(t(COV12),C22+NUG2))
    epsilon <- 1e-6
    C <- C + diag(epsilon, nrow(C))
    Cchol<-chol(C)
    Cinv<-chol2inv(Cchol)
    logD<-determinant(C)$modulus
    
    nloglikelihood <-(0.5 * logD + 0.5 * t(z) %*% Cinv %*% z+0.5*length(z)*log(2*pi))
    if(abs(nloglikelihood) == Inf || is.nan(nloglikelihood)){ nloglikelihood <- 1e+08}
    return(list(mlv=nloglikelihood,a1=a1,a2=a2,nu1=nu1,nu2=nu2,
                sigma1=sigma1,sigma2=sigma2,full.cov=C, df_train=un.grd.train))
  }
  
}


mle_ind_mlv<-function(pars)
{
  return(mle_ind_mat(p=pars,data_train)$mlv)
}

optim_indmat_loglik <- function(par){
  optim(par=par,
        fn = mle_ind_mlv,
        hessian=FALSE, 
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=20))
}

lmc.loglikelihood_allcomponents<-function(p,data_train){
  
  un.grd.train <- data_train[,c(1,2,3,4)]
  dist.mat <- rdist(un.grd.train[,-c(3,4)])
  z=c(un.grd.train$var1-mean(un.grd.train$var1) , un.grd.train$var2-mean(un.grd.train$var2))
  theta<-p
  a1<-theta[1]
  nu1<-theta[2]
  sigma1<-theta[3]
  a2<-theta[4]
  nu2<-theta[5]
  sigma2<-theta[6]
  b11<-theta[7]
  b12<-theta[8]
  b21<-theta[9]
  b22<-theta[10]
  nug1<-theta[11]
  nug2<-theta[12]
  ######## Putting hard constraints on the parameters #############
  if( a1<=0 | nu1<=0 | sigma1<=0 | a2<=0 | nu2<=0 | sigma2<=0 |nug1<0|nug2<0)
  {
    return(list(mlv=Inf))
  }
  else
  {
    n<-nrow(dist.mat)
    C11<-sigma1*matern(dist.mat,a1,nu1)
    C22<-sigma2*matern(dist.mat,a2,nu2)
    
    COV11<-(b11^2)*C11+(b12^2)*C22
    COV22<-(b21^2)*C11+(b22^2)*C22
    COV12<-(b11*b21)*C11+(b12*b22)*C22
    
    NUG1<-diag(nug1,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
    NUG2<-diag(nug2,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
    
    
    
    ############## Inverting C11 ##########
    C<-rbind(cbind(COV11+NUG1,COV12),cbind(t(COV12),COV22+NUG2))
    epsilon <- 1e-6
    C <- C + diag(epsilon, nrow(C))
    
    
    Cchol<-chol(C)
    Cinv<-chol2inv(Cchol)
    logD<-determinant(C)$modulus
    
    nloglikelihood <-(0.5 * logD + 0.5 * t(z) %*% Cinv %*% z+0.5*length(z)*log(2*pi))
    if(abs(nloglikelihood) == Inf || is.nan(nloglikelihood)){ nloglikelihood <- 1e+08}
    return(list(mlv=nloglikelihood,full.cov=C, Cchol=Cchol))
  }
}

lmc.loglikelihood<-function(par)
{
  
  return(lmc.loglikelihood_allcomponents(p=par,data_train)$mlv)
}

optim_lmc.loglikelihood <- function(par){
  optim(par=par,
        fn = lmc.loglikelihood,
        hessian=FALSE, 
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=80))
}

lmc.pred.summary<-function(estim.par, train_data, test_data)
{
  un.grd.train <- train_data[,c(1,2,3,4)]
  test.data <- test_data[,c(1,2,3,4)]
  full.data.coords = rbind(as.matrix(test.data[,c(1,2)]),as.matrix(un.grd.train[,c(1,2)]))
  dist.mat<-rdist(full.data.coords)
  p<-estim.par
  print(dim(full.data.coords))
  
  theta<-p
  a1<-theta[1]
  nu1<-theta[2]
  sigma1<-theta[3]
  a2<-theta[4]
  nu2<-theta[5]
  sigma2<-theta[6]
  b11<-theta[7]
  b12<-theta[8]
  b21<-theta[9]
  b22<-theta[10]
  nug1<-theta[11]
  nug2<-theta[12]
  
  C11<-sigma1*matern(dist.mat,a1,nu1)
  C22<-sigma2*matern(dist.mat,a2,nu2)
  
  COV11<-(b11^2)*C11+(b12^2)*C22
  COV22<-(b21^2)*C11+(b22^2)*C22
  COV12<-(b11*b21)*C11+(b12*b22)*C22
  
  NUG1<-diag(nug1,nrow = nrow(COV11),ncol=ncol(COV11))
  NUG2<-diag(nug2,nrow = nrow(COV22),ncol=ncol(COV22))
  
  ############## Inverting C11 ##########
  C<-rbind(cbind(COV11+NUG1,COV12),cbind(t(COV12),COV22+NUG2))
  
  epsilon <- 1e-6
  C <- C + diag(epsilon, nrow(C))
  
  print(dim(C))
  test.index = c(1:dim(test.data)[1],
                 (dim(un.grd.train)[1]+dim(test.data)[1]+1):(dim(un.grd.train)[1]+
                                                               dim(test.data)[1]+dim(test.data)[1]))
  
  C.test<-C[test.index,test.index]
  C.train<-C[-test.index,-test.index]
  C.test.train<-C[test.index,-test.index]
  
  print(dim(C.train))
  print(dim(C.test.train))
  
  z=c(un.grd.train$var1-mean(un.grd.train$var1) , un.grd.train$var2-mean(un.grd.train$var2))
  
  
  
  prediction<-C.test.train%*%solve(C.train)%*%z #Conditional mean of a Multivariate Gaussian
  cond_var_MAT<- C.test - C.test.train%*%solve(C.train)%*%t(C.test.train)
  
  prediction <- prediction + c(rep(mean(un.grd.train$var1),dim(test.data)[1]),
                               rep(mean(un.grd.train$var2),dim(test.data)[1]))
  return(list(pred = prediction, conditional_var = cond_var_MAT))
  
  
}

inverse <- function(f, lower, upper){
  function(y){
    uniroot(function(x){f(x) - y}, lower = lower, upper = upper, tol=1e-4)[1]
  }
}


# Gaussian.cond.cdf <- function(sim_num, index_loc){
#   pred_list_Gaussian = readRDS("pred_list_Gaussian.rds")
#   preds = pred_list_Gaussian[[sim_num]]$pred[c(index_loc,(120+index_loc)),1]
#   cond_var = pred_list_Gaussian[[sim_num]]$conditional_var[c(index_loc,(120+index_loc)),
#                                                            c(index_loc,(120+index_loc))]
#   value_range = seq(-20,40,0.2)
#   proba = rep(NA,length(value_range))
#   for(i in 1:length(value_range)){
#     proba[i] = pbinorm(value_range[i],value_range[i], mean1= preds[1],mean2=preds[2],
#                        var1=cond_var[1,1],var2=cond_var[2,2],cov12=cond_var[1,2])
#   }
#   return(list(proba = proba,x_axis = value_range))
# }