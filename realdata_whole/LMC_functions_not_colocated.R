library(geoR)
library(MASS)
library(fields)


mle_ind_mat <- function(p, data_train) {

  part1 <- data_train$part1
  part2 <- data_train$part2
  
  all_coords <- rbind(
    part1[, c("uj_lat", "uj_lon")],
    part2[, c("uj_lat", "uj_lon")]
  )
  
  dist.mat <- rdist(all_coords)
  
  z1 <- part1$Z1 - mean(part1$Z1)
  mu_z2 <- mean(c(part1$Z2, part2$Z2))
  z2 <- c(part1$Z2 - mu_z2, part2$Z2 - mu_z2)
  
  a1 <- p[1]
  nu1 <- p[2]
  sigma1 <- p[3]
  a2 <- p[4]
  nu2 <- p[5]
  sigma2 <- p[6]
  nug1 <- p[7]
  nug2 <- p[8]
  
  if (sum(p[1:6] < 0) != 0 || nug1 < 0 || nug2 < 0) {
    nloglikelihood <- 10000000
    return(list(mlv = nloglikelihood, params = NULL))
  } else {
    n_part1 <- nrow(part1)
    n_part2 <- nrow(part2)
    n_total <- n_part1 + n_part2
    
    dist_part1 <- dist.mat[1:n_part1, 1:n_part1]
    dist_part2 <- dist.mat[(n_part1+1):n_total, (n_part1+1):n_total]
    dist_part1_part2 <- dist.mat[1:n_part1, (n_part1+1):n_total]
    
    C11_part1 <- sigma1 * matern(dist_part1, a1, nu1)
    C11_part2 <- sigma1 * matern(dist_part2, a1, nu1)

    C22_part1 <- sigma2 * matern(dist_part1, a2, nu2)
    C22_part2 <- sigma2 * matern(dist_part2, a2, nu2)
    C22_part1_part2 <- sigma2 * matern(dist_part1_part2, a2, nu2)
    

    C <- rbind(
      cbind(C11_part1 + diag(nug1, n_part1),
            matrix(0, n_part1, n_part1),
            matrix(0, n_part1, n_part2)),
      cbind(matrix(0, n_part1, n_part1),
            C22_part1 + diag(nug2, n_part1),
            C22_part1_part2),
      cbind(matrix(0, n_part2, n_part1),
            t(C22_part1_part2),
            C22_part2 + diag(nug2, n_part2))
    )
    
    C <- (C + t(C)) / 2
    
    epsilon <- 1e-8
    C <- C + diag(epsilon, nrow(C))

    eigenvals <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
    if (any(eigenvals <= 0)) {
      nloglikelihood <- 1e+08
      return(list(mlv = nloglikelihood, params = NULL))
    }
    
    Cchol <- chol(C)
    Cinv <- chol2inv(Cchol)
    logD <- determinant(C)$modulus
    
    z <- c(z1, z2)
    
    nloglikelihood <- (0.5 * logD + 0.5 * t(z) %*% Cinv %*% z + 0.5 * length(z) * log(2 * pi))
    if (abs(nloglikelihood) == Inf || is.nan(nloglikelihood)) {
      nloglikelihood <- 1e+08
    }
    
    return(list(mlv = nloglikelihood, a1 = a1, a2 = a2, nu1 = nu1, nu2 = nu2,
                sigma1 = sigma1, sigma2 = sigma2, full.cov = C, df_train = data_train))
  }
}

mle_ind_mlv <- function(pars, data_train) {
  mle_ind_mat(p = pars, data_train = data_train)$mlv
}

optim_indmat_loglik <- function(par, data_train) {
  optim(par = par,
        fn = mle_ind_mlv,
        data_train = data_train,
        hessian = FALSE, 
        control = list(trace = 6,
                       pgtol = 0,
                       parscale = rep(0.1, length(par)),
                       maxit = 20))
}

lmc.loglikelihood_allcomponents <- function(p, data_train) {
  part1 <- data_train$part1
  part2 <- data_train$part2
  
  all_coords <- rbind(
    part1[, c("uj_lat", "uj_lon")],
    part2[, c("uj_lat", "uj_lon")]
  )
  
  dist.mat <- rdist(all_coords)
  
  z1 <- part1$Z1 - mean(part1$Z1)
  mu_z2 <- mean(c(part1$Z2, part2$Z2))
  z2 <- c(part1$Z2 - mu_z2, part2$Z2 - mu_z2)
  
  theta <- p
  a1 <- theta[1]
  nu1 <- theta[2]
  sigma1 <- theta[3]
  a2 <- theta[4]
  nu2 <- theta[5]
  sigma2 <- theta[6]
  b11 <- theta[7]
  b12 <- theta[8]
  b21 <- theta[9]
  b22 <- theta[10]
  nug1 <- theta[11]
  nug2 <- theta[12]
  
  if (a1 <= 0 | nu1 <= 0 | sigma1 <= 0 | a2 <= 0 | nu2 <= 0 | sigma2 <= 0 | nug1 < 0 | nug2 < 0) {
    return(list(mlv = Inf))
  } else {
    n_part1 <- nrow(part1)
    n_part2 <- nrow(part2)
    n_total <- n_part1 + n_part2
    
    dist_part1 <- dist.mat[1:n_part1, 1:n_part1]
    dist_part2 <- dist.mat[(n_part1+1):n_total, (n_part1+1):n_total]
    dist_part1_part2 <- dist.mat[1:n_part1, (n_part1+1):n_total]

    C11_part1 <- sigma1 * matern(dist_part1, a1, nu1)
    C11_part2 <- sigma1 * matern(dist_part2, a1, nu1)
    C11_part1_part2 <- sigma1 * matern(dist_part1_part2, a1, nu1)
    
    C22_part1 <- sigma2 * matern(dist_part1, a2, nu2)
    C22_part2 <- sigma2 * matern(dist_part2, a2, nu2)
    C22_part1_part2 <- sigma2 * matern(dist_part1_part2, a2, nu2)
    
    COV11_part1 <- (b11^2) * C11_part1 + (b12^2) * C22_part1
    COV22_part1 <- (b21^2) * C11_part1 + (b22^2) * C22_part1
    COV22_part2 <- (b21^2) * C11_part2 + (b22^2) * C22_part2
    
    COV12_part1 <- (b11 * b21) * C11_part1 + (b12 * b22) * C22_part1
    COV12_part1_part2 <- (b11 * b21) * C11_part1_part2 + (b12 * b22) * C22_part1_part2
    
    C <- rbind(
      cbind(COV11_part1 + diag(nug1, n_part1),  COV12_part1,               COV12_part1_part2),
      cbind(t(COV12_part1),                     COV22_part1 + diag(nug2, n_part1),
            (b21^2) * C11_part1_part2 + (b22^2) * C22_part1_part2),
      cbind(t(COV12_part1_part2),               t((b21^2) * C11_part1_part2 + (b22^2) * C22_part1_part2),
            COV22_part2 + diag(nug2, n_part2))
    )

    C <- (C + t(C)) / 2
    
    epsilon <- 1e-8
    C <- C + diag(epsilon, nrow(C))
    
    eigenvals <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
    if (any(eigenvals <= 0)) {
      return(list(mlv = Inf))
    }
    
    Cchol <- chol(C)
    Cinv <- chol2inv(Cchol)
    logD <- determinant(C)$modulus
    
    z <- c(z1, z2)
    
    nloglikelihood <- (0.5 * logD + 0.5 * t(z) %*% Cinv %*% z + 0.5 * length(z) * log(2 * pi))
    if (abs(nloglikelihood) == Inf || is.nan(nloglikelihood)) {
      nloglikelihood <- 1e+08
    }
    
    return(list(mlv = nloglikelihood, full.cov = C, Cchol = Cchol))
  }
}

lmc.loglikelihood <- function(par, data_train) {
  lmc.loglikelihood_allcomponents(p = par, data_train = data_train)$mlv
}

optim_lmc.loglikelihood <- function(par, data_train) {
  optim(par = par,
        fn = lmc.loglikelihood,
        data_train = data_train,
        hessian = FALSE, 
        control = list(trace = 6,
                       pgtol = 0,
                       parscale = rep(0.1, length(par)),
                       maxit = 80))
}


lmc.pred.summary <- function(estim.par, train_data, test_data, var_type = c("latent","obs")) {
  var_type <- match.arg(var_type)
  
  stopifnot(length(estim.par) == 12)
  if (any(estim.par[1:6] <= 0) || any(estim.par[11:12] < 0))
    stop("All parameters must be positive (the nugget may be 0)")
  if (!all(c("part1","part2") %in% names(train_data)))
    stop("train_data must contain part1 and part2")
  if (!all(c("uj_lat","uj_lon","Z2","Z1") %in% names(test_data)))
    stop("test_data must contain uj_lat, uj_lon, Z2, and Z1")

  part1 <- train_data$part1   # 有 Z1,Z2
  part2 <- train_data$part2   # 仅 Z2
  test.data <- test_data[, c("uj_lat","uj_lon","Z2","Z1")]
  
  coords <- rbind(
    as.matrix(test.data[,c("uj_lat","uj_lon")]),
    as.matrix(part1[,c("uj_lat","uj_lon")]),
    as.matrix(part2[,c("uj_lat","uj_lon")])
  )
  D <- fields::rdist(coords)
  
  th <- estim.par
  a1 <- th[1]; nu1 <- th[2]; s1 <- th[3]
  a2 <- th[4]; nu2 <- th[5]; s2 <- th[6]
  b11 <- th[7]; b12 <- th[8]; b21 <- th[9]; b22 <- th[10]
  nug1 <- th[11]; nug2 <- th[12]
  
  n_test  <- nrow(test.data)
  n_p1    <- nrow(part1)
  n_p2    <- nrow(part2)
  
  D_tt    <- D[1:n_test, 1:n_test, drop=FALSE]
  D_t_p1  <- D[1:n_test, (n_test+1):(n_test+n_p1), drop=FALSE]
  D_t_p2  <- D[1:n_test, (n_test+n_p1+1):(n_test+n_p1+n_p2), drop=FALSE]
  D_p1_p1 <- D[(n_test+1):(n_test+n_p1), (n_test+1):(n_test+n_p1), drop=FALSE]
  D_p1_p2 <- D[(n_test+1):(n_test+n_p1), (n_test+n_p1+1):(n_test+n_p1+n_p2), drop=FALSE]
  D_p2_p2 <- D[(n_test+n_p1+1):(n_test+n_p1+n_p2),
               (n_test+n_p1+1):(n_test+n_p1+n_p2), drop=FALSE]
  
  K11_tt    <- s1 * matern(D_tt,    a1, nu1)
  K22_tt    <- s2 * matern(D_tt,    a2, nu2)
  K11_t_p1  <- s1 * matern(D_t_p1,  a1, nu1)
  K22_t_p1  <- s2 * matern(D_t_p1,  a2, nu2)
  K11_t_p2  <- s1 * matern(D_t_p2,  a1, nu1)
  K22_t_p2  <- s2 * matern(D_t_p2,  a2, nu2)
  K11_p1_p1 <- s1 * matern(D_p1_p1, a1, nu1)
  K22_p1_p1 <- s2 * matern(D_p1_p1, a2, nu2)
  K11_p1_p2 <- s1 * matern(D_p1_p2, a1, nu1)
  K22_p1_p2 <- s2 * matern(D_p1_p2, a2, nu2)
  K11_p2_p2 <- s1 * matern(D_p2_p2, a1, nu1)
  K22_p2_p2 <- s2 * matern(D_p2_p2, a2, nu2)
  
  # LMC 组合
  C11_tt    <- (b11^2)*K11_tt    + (b12^2)*K22_tt
  C22_tt    <- (b21^2)*K11_tt    + (b22^2)*K22_tt
  C12_tt    <- (b11*b21)*K11_tt  + (b12*b22)*K22_tt
  
  C11_t_p1  <- (b11^2)*K11_t_p1  + (b12^2)*K22_t_p1
  C22_t_p1  <- (b21^2)*K11_t_p1  + (b22^2)*K22_t_p1
  C12_t_p1  <- (b11*b21)*K11_t_p1+ (b12*b22)*K22_t_p1
  
  C11_t_p2  <- (b11^2)*K11_t_p2  + (b12^2)*K22_t_p2
  C22_t_p2  <- (b21^2)*K11_t_p2  + (b22^2)*K22_t_p2
  C12_t_p2  <- (b11*b21)*K11_t_p2+ (b12*b22)*K22_t_p2
  
  C11_p1_p1 <- (b11^2)*K11_p1_p1 + (b12^2)*K22_p1_p1
  C22_p1_p1 <- (b21^2)*K11_p1_p1 + (b22^2)*K22_p1_p1
  C12_p1_p1 <- (b11*b21)*K11_p1_p1+ (b12*b22)*K22_p1_p1
  
  C12_p1_p2 <- (b11*b21)*K11_p1_p2+ (b12*b22)*K22_p1_p2
  C22_p1_p2 <- (b21^2)*K11_p1_p2  + (b22^2)*K22_p1_p2
  C22_p2_p2 <- (b21^2)*K11_p2_p2  + (b22^2)*K22_p2_p2
  
  C <- rbind(
    cbind(C11_tt,                  C12_tt,                        C11_t_p1,                         C12_t_p1,                         C12_t_p2),
    cbind(t(C12_tt),               C22_tt,                        C12_t_p1,                         C22_t_p1,                         C22_t_p2),
    cbind(t(C11_t_p1),             t(C12_t_p1),                   C11_p1_p1,                        C12_p1_p1,                        C12_p1_p2),
    cbind(t(C12_t_p1),             t(C22_t_p1),                   t(C12_p1_p1),                     C22_p1_p1,                         C22_p1_p2),
    cbind(t(C12_t_p2),             t(C22_t_p2),                   t(C12_p1_p2),                     t(C22_p1_p2),                      C22_p2_p2)
  )
  C <- (C + t(C)) / 2
  
  test.idx  <- 1:(2*n_test)
  train.idx <- (2*n_test + 1):nrow(C)
  
  i1 <- (2*n_test + 1):(2*n_test + n_p1)                     # Z1_p1
  i2 <- (2*n_test + n_p1 + 1):(2*n_test + 2*n_p1)            # Z2_p1
  i3 <- (2*n_test + 2*n_p1 + 1):(2*n_test + 2*n_p1 + n_p2)   # Z2_p2
  C[i1,i1] <- C[i1,i1] + diag(nug1, n_p1)
  C[i2,i2] <- C[i2,i2] + diag(nug2, n_p1)
  C[i3,i3] <- C[i3,i3] + diag(nug2, n_p2)
  
  C <- (C + t(C)) / 2
  C <- C + diag(1e-8, nrow(C))
  
  C_tt <- C[test.idx,  test.idx,  drop=FALSE]
  C_xx <- C[train.idx, train.idx, drop=FALSE]
  C_tx <- C[test.idx,  train.idx, drop=FALSE]
  
  mu1 <- mean(part1$Z1)
  mu2 <- mean(c(part1$Z2, part2$Z2))
  z_train <- c(part1$Z1 - mu1,
               c(part1$Z2 - mu2, part2$Z2 - mu2))
  

  L <- chol(C_xx)
  alpha <- backsolve(L, forwardsolve(t(L), z_train))
  pred_centered <- C_tx %*% alpha
  pred <- pred_centered + c(rep(mu1, n_test), rep(mu2, n_test))
  
  # conditional variance
  V <- backsolve(L, forwardsolve(t(L), t(C_tx)))
  cond_var_latent <- C_tt - C_tx %*% V
  cond_var_latent <- (cond_var_latent + t(cond_var_latent)) / 2
  diag(cond_var_latent) <- pmax(diag(cond_var_latent), 0)
  
  if (var_type == "obs") {
    cond_var <- cond_var_latent
    idx1 <- 1:n_test
    idx2 <- (n_test+1):(2*n_test)
    diag(cond_var)[idx1] <- diag(cond_var)[idx1] + nug1
    diag(cond_var)[idx2] <- diag(cond_var)[idx2] + nug2
  } else {
    cond_var <- cond_var_latent
  }
  
  list(pred = pred, conditional_var = cond_var)
}




inverse <- function(f, lower, upper) {
  function(y) {
    uniroot(function(x) { f(x) - y }, lower = lower, upper = upper, tol = 1e-4)[1]
  }
}
