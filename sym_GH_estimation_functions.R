#' Modelo simétrico con lambda conocido y eta estimado con
#' loglik perfilada

library(tidyverse) # los graficos estan con ggplot y uso dplyr tambien 


nam <- function(p, q) {
  name_theta <- NULL
  for (i in 1:q) {
    name_a <- paste("a_", i, sep = "")
    name_theta <- c(name_theta, name_a)
  }
  
  for (k in 1:p) {
    for (l in 1:q) {
      name_b <- paste("b_", l, k, sep = "")
      name_theta <- c(name_theta, name_b)
    }
  }
  
  name_theta <- c(name_theta, "phi")
  
  for (i in 1:p) {
    name_mu <- paste("mu_", i, sep = "")
    name_theta <- c(name_theta, name_mu)
  }
  
  for (k in 1:p) {
    for (l in 1:p) {
      if (l >= k) {
        name_s <- paste("sigma_", l, k, sep = "")
        name_theta <- c(name_theta, name_s)
      }
    }
  }
  #name_theta <- c(name_theta, "lambda", "eta")
  name_theta
}

# simulacion --------------------------------------------------------------


# simulacion de datos bajo el modelo normal
sim_Ndata <- function(n, p, q, a, B, phi, mu_x, Sigma_x, testB = FALSE,vars.pos = c(1)) {
  r <- p + q
  # True x values, from a N(mu_x, Sigma_x)
  x <- matrix(, ncol = p, nrow = n)
  for (i in 1:n) {
    x[i, ] <- mvtnorm::rmvnorm(1, mu_x, Sigma_x, method = "chol")
  }
  
  # Simulate errors
  Sigma_e <- diag(phi, r, r)
  
  # Errors vector, they are independent and come from a multivariate normal
  e <- matrix(, ncol = r, nrow = n)
  for (i in 1:n) {
    e[i, ] <- mvtnorm::rmvnorm(1, rep(0, r), Sigma_e, method = "chol")
  }
  e1 <- as.matrix(e[, 1:p])
  e2 <- as.matrix(e[, -(1:p)])
  
  # Simulation of Y
  Y <- matrix(, ncol = q, nrow = n)
  for (i in 1:n) {
    Y[i, ] <- a + B %*% x[i, ] + e2[i, ]
  }
  
  # Simulation of Y under H0
  Y_h0 <- matrix(, ncol = q, nrow = n)
  B_h0 <- B
  B_h0[, vars.pos] <- 0
  
  for (i in 1:n) {
    Y_h0[i, ] <- a + B_h0 %*% x[i, ] + e2[i, ]
  }
  
  # X observed
  X <- x + e1
  
  # Simulation of Z
  alpha <- matrix(c(rep(0, p), a), ncol = 1)
  Lambda <- rbind(diag(1, p), B)
  Z <- matrix(, ncol = r, nrow = n)
  for (i in 1:n) {
    Z[i, ] <- alpha + Lambda %*% x[i, ] + e[i, ]
  }
  
  # eta
  eta <- alpha + Lambda %*% matrix(mu_x, ncol = 1)
  
  # psi
  psi <- Lambda %*% Sigma_x %*% t(Lambda) + diag(phi, r, r)
  
  theta_sim <- c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x))
  names(theta_sim) <- nam(p, q)
  
  # Resultados
  data = list(
    "X" = X,
    "Y" = Y,
    "Z" = Z,
    "a" = a,
    "B" = B,
    "phi" = phi,
    "mu_x" = mu_x,
    "Sigma_x" = Sigma_x,
    "eta" = eta,
    "psi" = psi,
    "theta" = theta_sim
    
  )
  if (testB == FALSE) {
    return(data)
  } else{
    data_h0 = list(
      "X" = X,
      "Y" = Y,
      "Y_h0" = Y_h0,
      "Z" = Z,
      "a" = a,
      "B" = B,
      "phi" = phi,
      "mu_x" = mu_x,
      "Sigma_x" = Sigma_x,
      "eta" = eta,
      "psi" = psi,
      "theta" = theta_sim
    )
    return(data_h0)
  }
}


sim_SGH_data <- function(n, p, q,
                         a, B, phi, mu_x, Sigma_x, lambda = -1/2, eta, 
                         type = "S-NIG",
                         test_a = FALSE) {
  r <- p + q
  x_e <- matrix(, ncol = (p+r), nrow = n)
  gamma_x = rep(0,p) #caso simetrico
  # S-NIG ----
  if(type == "S-NIG") { # NIG (lambda = -1/2)
    lambda <- -1/2
    NIGS_dist <- ghyp::NIG(alpha.bar = eta,
                           mu = c(mu_x,rep(0,r)),
                           sigma = as.matrix(Matrix::bdiag(Sigma_x, phi * diag(1,r))),
                           gamma= c(gamma_x,rep(0,r)))
    for (i in 1:n) {
      x_e[i, ] <- ghyp::rghyp(1, NIGS_dist)
    }
    theta <- c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x), lambda, eta)
  }
  # S-HYP ----
  if(type == "S-HYP"){ # HYP (lambda = 1)
    lambda <- 1
    hyp_dist <- ghyp::ghyp(lambda = lambda,
                           chi = eta,
                           psi = eta,
                           mu = c(mu_x, rep(0,r)),
                           sigma = as.matrix(Matrix::bdiag(Sigma_x, phi * diag(1,r))),
                           gamma= c(gamma_x,rep(0,r)))
    for (i in 1:n) {
      x_e[i, ] <- ghyp::rghyp(1, hyp_dist)
    }
    theta <- c(a, B, phi, mu_x,  matrixcalc::vech(Sigma_x), lambda, eta)
  }
  # S-GH ----
  if(type == "S-GH"){ # GH (lambda conocido)
    hyp_dist <- ghyp::ghyp(lambda = lambda,
                           chi = eta,
                           psi = eta,
                           mu = c(mu_x, rep(0,r)),
                           sigma = as.matrix(Matrix::bdiag(Sigma_x, phi * diag(1,r))),
                           gamma= c(gamma_x,rep(0,r)))
    for (i in 1:n) {
      x_e[i, ] <- ghyp::rghyp(1, hyp_dist)
    }
    theta <- c(a, B, phi, mu_x,  matrixcalc::vech(Sigma_x), lambda, eta)
  }
  
  alpha <- matrix(c(rep(0, p), a), ncol = 1)
  Lambda <- rbind(diag(1, p), B)
  #' Z_i = alpha + Lambda*x_i + e_i 
  Z <- matrix(, ncol = r, nrow = n)
  for (i in 1:n) {
    Z[i, ] <- alpha + (cbind(Lambda,diag(r)) %*% x_e[i,])  
    
  }
  X <- Z[,1:p]
  Y <- Z[,-c(1:p)]
  
  
  # Y bajo H_0
  Z_h0 <- matrix(, ncol = r, nrow = n)
  for (i in 1:n) {
    Z_h0[i, ] <- (cbind(Lambda,diag(r)) %*% x_e[i,])  
    
  }
  Y_h0 <- Z_h0[,-c(1:p)]
  
  # Results
  data = list("X" = X, "Y" = Y, "Z" = Z, "theta" = theta)
  
  if (test_a == FALSE){
    return(data)
  } else {
    data = list("X" = X, "Y" = Y, "Z" = Z, "theta" = theta, "Y_h0" = Y_h0)
    return(data)
  }
  
}

# Verosimilitud ----

# log verosimilitud del modelo normal 
logL <- function(theta, X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Z <- cbind(X,Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2
  
  # theta --> (a,B,phi,mu,Sigma)
  theta <- as.vector(theta)
  a <- matrix(theta[1:q], ncol = 1, nrow = q)
  B <- matrix(theta[(q + 1):(q + p * q)], nrow = q, ncol = p)
  phi <- theta[q + p * q + 1]
  mu_x <- theta[(q + p * q + 2):(q + p * q + 1 + p)]
  Sigma_x <- ks::invvech(theta[(q + p * q + 2 + p):d])
  Sigma_x <- 0.5 * (Sigma_x + t(Sigma_x))
  
  # eta and psi
  Lambda <- rbind(diag(1, p, p), B)
  alpha <- c(rep(0,p),a)
  eta <- c(alpha + Lambda %*% mu_x)
  psi <- Lambda %*% tcrossprod(Sigma_x, Lambda) + diag(phi, r, r)
  psi <- 0.5 * (psi + t(psi))
  
  # logL MRMME
  cte <- - r * 0.5 * log(2 * pi)
  det2 <- det(phi*diag(1,r,r)) * det(diag(1,p,p) + (1/phi) * crossprod(Lambda) %*% Sigma_x)
  log_det <- log(det2)
  delta_i <- stats::mahalanobis(Z,eta,psi)
  logL_i <- cte - log_det/2 - 1/2 * delta_i
  logL = sum(logL_i)
  logL
}

# log verosimilitud de cada observacion en el modelo normal
logLi_N <- function(theta, X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Z <- cbind(X,Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2
  
  # theta --> (a,B,phi,mu,Sigma)
  theta <- as.vector(theta)
  a <- matrix(theta[1:q], ncol = 1, nrow = q)
  B <- matrix(theta[(q + 1):(q + p * q)], nrow = q, ncol = p)
  phi <- theta[q + p * q + 1]
  mu_x <- theta[(q + p * q + 2):(q + p * q + 1 + p)]
  Sigma_x <- ks::invvech(theta[(q + p * q + 2 + p):d])
  Sigma_x <- 0.5 * (Sigma_x + t(Sigma_x))
  
  # eta and psi
  Lambda <- rbind(diag(1, p, p), B)
  alpha <- c(rep(0,p),a)
  eta <- c(alpha + Lambda %*% mu_x)
  psi <- Lambda %*% tcrossprod(Sigma_x, Lambda) + diag(phi, r, r)
  psi <- 0.5 * (psi + t(psi))
  
  # logL MRMME
  cte <- - r * 0.5 * log(2 * pi)
  det2 <- det(phi*diag(1,r,r)) * det(diag(1,p,p) + (1/phi) * crossprod(Lambda) %*% Sigma_x)
  log_det <- log(det2)
  delta_i <- stats::mahalanobis(Z,eta,psi)
  logL_i <- cte - log_det/2 - 1/2 * delta_i
  logL_i
}


lLz <- function(X, Y, 
                theta,
                type ="S-NIG") {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  Z <- cbind(X, Y)
  d <- (p + 1) * (2 + p + 2 * q) / 2
  
    # theta --> (a,B,phi,mu,Sigma)
  theta <- c(theta)
  a <- matrix(theta[1:q], ncol = 1, nrow = q)
  B <- matrix(theta[(q + 1):(q + p * q)], nrow = q, ncol = p)
  phi <- theta[q + p * q + 1]
  mu_x <- theta[(q + p * q + 2):(q + p * q + 1 + p)]
  Sigma_x <- ks::invvech(theta[(q + p * q + 2 + p):d])
  Sigma_x <- 0.5 * (Sigma_x + t(Sigma_x))
  lambda <- theta[(d + 1)]
  eta <- theta[(d + 2)]
  
  alpha <- c(rep(0, p), a)
  Lambda <- rbind(diag(1, p), B)
  Sigma_e <- diag(phi, r, r)
  gamma_x <- rep(0,p) #caso simetrico
  
  # S-NIG ----
  if (type == "S-NIG"){
      lambda <- -1/2
      mu_z <- c(alpha + Lambda %*% mu_x)
      Sigma_z <- Lambda %*% tcrossprod(Sigma_x, Lambda) + Sigma_e
      Sigma_z <- 0.5 * (Sigma_z + t(Sigma_z))
      gamma_z <- Lambda %*% gamma_x
      ghyp_params <- ghyp::ghyp(
        lambda = lambda,
        chi = eta,
        psi = eta,
        mu = mu_z,
        sigma = Sigma_z,
        gamma = gamma_z
      )
  }
  # S-HYP ----
  if (type == "S-HYP"){
      lambda <- 1
      mu_z <- c(alpha + Lambda %*% mu_x)
      Sigma_z <- Lambda %*% tcrossprod(Sigma_x, Lambda) + Sigma_e
      Sigma_z <- 0.5 * (Sigma_z + t(Sigma_z))
      gamma_z <- Lambda %*% gamma_x
      ghyp_params <- ghyp::ghyp(
        lambda = lambda,
        chi = eta,
        psi = eta,
        mu = mu_z,
        sigma = Sigma_z,
        gamma = gamma_z
      )
  }
  # S-GH ----
  if (type == "S-GH"){
      lambda <- lambda
      mu_z <- c(alpha + Lambda %*% mu_x)
      Sigma_z <- Lambda %*% tcrossprod(Sigma_x, Lambda) + Sigma_e
      Sigma_z <- 0.5 * (Sigma_z + t(Sigma_z))
      gamma_z <- Lambda %*% gamma_x
      ghyp_params <-  ghyp::ghyp(
        lambda = lambda,
        chi = eta,
        psi = eta,
        mu = mu_z,
        sigma = Sigma_z,
        gamma = gamma_z
      )
  }
  logL_i <- 0
  for (i in 1:n) {
    logL_i[i] <- ghyp::dghyp(Z[i, ], object = ghyp_params, logvalue = TRUE)
  }
  return(list(logL = sum(logL_i), logL_i = logL_i))
}

#' Para garantizar robustez el parametro eta se estima por verosimilitud 
#' perfilada y lambda se asume conocido

# Ajuste ---------------------------------------------------------------------

# esta funcion ajuesta el modelo para un eta dado
fit_SGH_MRMME <-
  function(X,
           Y,
           # a,
           # B,
           # phi,
           # mu_x,
           # Sigma_x,
           lambda = -1/2,
           eta, # este es fijo y se escoje de una grilla
           type = "S-NIG",
           StopCrit = "logL",
           MaxIter = 500,
           precision = 1e-10) {
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    n <- dim(Y)[1]
    q <- dim(Y)[2]
    p <- dim(X)[2]
    r <- p + q
    Z <- cbind(X, Y)
    d <- (p + 1) * (2 + p + 2 * q) / 2
    
    fit_normal <- MRMME::mrmme(Y=Y, X=X) #ajuste normal para los valores iniciales
    a <- fit_normal$a
    B <- fit_normal$B
    phi <- fit_normal$phi
    mu_x <- fit_normal$mu.x
    Sigma_x <- fit_normal$Sigma.x

    
    # mod_lm <- stats::lm(Y ~ X)
    # a <- matrix(stats::coef(mod_lm)[1,],ncol=1) # lm
    # B <- t(stats::coef(mod_lm)[-1,]) # lm
    # phi <- 1
    # mu_x <- colMeans(X) # observed
    # Sigma_x <- stats::cov(X) # observed
    
    
    
    if (type == "S-NIG") {
      # S-NIG ----
      lambda <- -1/2
      gamma_x <- rep(0, p)
      eta <- eta
      theta <- c(a,
                B,
                phi,
                mu_x,
                matrixcalc::vech(Sigma_x),
                lambda,
                eta)
      lk0 <-
          lLz(
            X,
            Y,
            theta,
            type = type
          )$logL
        loglik <- NULL
        criterio <- 1
        crits <- NULL
        count <- 0
        
        while (criterio > precision) {
          sum_m1 <- 0
          sum_m2 <- matrix(0, nrow = p, ncol = 1)
          sum_m3 <- matrix(0, nrow = p, ncol = p)
          mi2 <- matrix(0, nrow = p, ncol = 1)
          mi3 <- matrix(0, nrow = p, ncol = p)
          s1 <- rep(0, q) #\sum_{i=1}^{n}m_{i1}\bY_i
          s2 <- matrix(0, nrow = q, ncol = p)   #\sum_{i=1}^{n}\bY_i \bm_{i2}^{\top} 
          s3 <- 0 # \sum_{i=1}^{n} \Big\{ m_{i1}\Big[\bX_i^{\top} \bX_i + 
                  # (\bY_i-\hat{\ba})^{\top}(\bY_i-\hat{\ba})\Big] - \bX_i^{\top}\bm_{i2} 
                  #- (\bY_i - \hat{\ba})^{\top}\hat{\bB} \bm_{i2} - \bm_{i2}^{\top}\bX_i 
                  #- \bm_{i2}^{\top}\hat{\bB}^{\top}(\bY-\hat{\ba}) + \text{tr}(\bm_{i3}) 
                  # + \text{tr}(\bm_{i3}\hat{\bB}^{\top}\hat{\bB})\Big\}
          s4 <- matrix(0, nrow = p, ncol = p) # \sum_{i=1}^{n}\Big( \bm_{i3} - 
                                              # \bm_{i2}\hat{\bmu}_x^{\top} - 
                                              # \hat{\bmu}_x\bm_{i2}^{\top} + 
                                              # m_{i1}\hat{\bmu}_x\hat{\bmu}_x^{\top}\Big)
          count <- count + 1
          #print(count)
          alpha <- c(rep(0, p), a)
          Lambda <- rbind(diag(1, p), B)
          mu_z <- c(alpha + Lambda %*% mu_x)
          Sigma_z <- Lambda %*% Sigma_x %*% t(Lambda) + diag(phi, r, r)
          Sigma_z <- 0.5 * (Sigma_z + t(Sigma_z))
          Sigma_z_inv <- solve(Sigma_z, tol = 1e-30)
          Delta <- solve(Sigma_x, tol = 1e-30) + ((1 / phi) * crossprod(Lambda))
          Delta_inv <- solve(Delta, tol = 1e-30)
          #A <- diag(p) - (1 / phi) * Delta_inv %*% crossprod(Lambda)
          lambda_vz <- lambda - r / 2
          for (i in 1:n) {
            delta_i <-  t(Z[i, ] - mu_z) %*% Sigma_z_inv %*% (Z[i, ] - mu_z)
            m1_i <-
              ghyp::Egig(
                lambda = lambda_vz,
                chi = delta_i + eta,
                psi = eta,
                func = c("1/x")
              )
            b_i <-  (1 / phi) * Delta_inv %*% t(Lambda) %*% (Z[i, ] - mu_z)
            m2_i <- m1_i * (mu_x + b_i)
            m3_i <- Delta_inv + m1_i * tcrossprod(mu_x + b_i)
            m3_i <- 0.5 * (m3_i + t(m3_i)) # porque es simetrica la matriz
            sum_m1 <- sum_m1 + m1_i
            sum_m2 <- sum_m2 + m2_i
            sum_m3 <- sum_m3 + m3_i
            s1 <- s1 + (m1_i * Y[i, ])  
            s2 <- s2 + (Y[i, ] %*% t(m2_i))
            s3 <- s3 + c(
                  m1_i * (crossprod(X[i, ]) + crossprod(Y[i, ] - a)) -
                  crossprod(X[i, ], m2_i) -
                  t(Y[i, ] - a) %*% B %*% m2_i -
                  crossprod(m2_i, X[i, ]) -
                  t(m2_i) %*% t(B) %*% (Y[i, ] - a) +
                  sum(diag(m3_i)) +
                  sum(diag(m3_i %*% t(B) %*% B))
              )
            s4 <- s4 + m3_i -
              tcrossprod(m2_i, mu_x) -
              tcrossprod(mu_x, m2_i) +
              m1_i * tcrossprod(mu_x)
          }
          # parameters update
          a <- (s1 / sum_m1) - (B %*% (sum_m2 / sum_m1))
          B <- (s2 - (s1 / sum_m1) %*% t(sum_m2)) %*% solve((sum_m3 - (sum_m2 / sum_m1) %*% t(sum_m2)), tol = 1e-30)
          phi <- (1 / (n * r)) * s3
          mu_x <- sum_m2 / c(sum_m1)
          Sigma_x <- (1 / n) * s4
          
          theta <-
            c(a,
              B,
              phi,
              mu_x,
              matrixcalc::vech(Sigma_x),
              lambda,
              eta)
          
          #criterio
          lk <- lLz(
            X,
            Y,
            theta,
            type = type
          )$logL
          loglik[count] <- lk
          
          if (StopCrit == "logL") {
            lk1 <- lk
            criterio <- abs(lk1 / lk0 - 1)
            lk0 <- lk1
            if (count == MaxIter) {
              criterio <- 1e-11
            }
            # print(paste("criterio de convergencia simple:", criterio))
            # print(paste("loglikelihood", lk))
            # print(c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x), eta))
          }
          if (StopCrit == "Aitken") {
            if (count < 2) {
              lk1 <- lk
              criterio <- abs(lk1 / lk0 - 1)
            }
            if (count == 2) {
              lk2 <- lk
              c1 <- (lk2 - lk1) / (lk1 - lk0)
              l_A2 <- lk1 + ((lk2 - lk1) / (1 - c1))
              criterio <- abs(lk2 / lk1 - 1)
            }
            if (count > 2) {
              lk3 <- lk
              c <- (lk3 - lk2) / (lk2 - lk1)
              l_A3 <- lk2 + ((lk3 - lk2) / (1 - c))
              criterio <- abs(l_A3 - l_A2)
              l_A2 <- l_A3
              lk1 <- lk2
              lk2 <- lk3
            }
            if (count == MaxIter) {
              criterio <- 1e-11
            }
            # print(paste("criterio de convergencia Aitken:", criterio))
            # print(paste("loglikelihood", lk))
            # print(c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x), eta))
            
          }
          crits[count] <- criterio
          
          # theta <-
          #   c(a,
          #     B,
          #     phi,
          #     mu_x,
          #     matrixcalc::vech(Sigma_x),
          #     lambda,
          #     eta)
        }
    }
    if (type == "S-HYP") {
      # S-HYP ----
      lambda <- 1
      gamma_x <- rep(0, p)
      eta <- eta
      theta <-c(a,
                 B,
                 phi,
                 mu_x,
                 matrixcalc::vech(Sigma_x),
                 lambda,
                 eta)
      lk0 <-
        lLz(
          X,
          Y,
          theta,
          type = type
        )$logL
      loglik <- NULL
      criterio <- 1
      crits <- NULL
      count <- 0
      while (criterio > precision) {
        sum_m1 <- 0
        sum_m2 <- matrix(0, nrow = p, ncol = 1)
        sum_m3 <- matrix(0, nrow = p, ncol = p)
        mi2 <- matrix(0, nrow = p, ncol = 1)
        mi3 <- matrix(0, nrow = p, ncol = p)
        s1 <- rep(0, q) #\sum_{i=1}^{n}m_{i1}\bY_i
        s2 <- matrix(0, nrow = q, ncol = p)   #\sum_{i=1}^{n}\bY_i \bm_{i2}^{\top} 
        s3 <- 0 # \sum_{i=1}^{n} \Big\{ m_{i1}\Big[\bX_i^{\top} \bX_i + 
        # (\bY_i-\hat{\ba})^{\top}(\bY_i-\hat{\ba})\Big] - \bX_i^{\top}\bm_{i2} 
        #- (\bY_i - \hat{\ba})^{\top}\hat{\bB} \bm_{i2} - \bm_{i2}^{\top}\bX_i 
        #- \bm_{i2}^{\top}\hat{\bB}^{\top}(\bY-\hat{\ba}) + \text{tr}(\bm_{i3}) 
        # + \text{tr}(\bm_{i3}\hat{\bB}^{\top}\hat{\bB})\Big\}
        s4 <- matrix(0, nrow = p, ncol = p) # \sum_{i=1}^{n}\Big( \bm_{i3} - 
        # \bm_{i2}\hat{\bmu}_x^{\top} - 
        # \hat{\bmu}_x\bm_{i2}^{\top} + 
        # m_{i1}\hat{\bmu}_x\hat{\bmu}_x^{\top}\Big)
        count <- count + 1
        print(count)
        alpha <- c(rep(0, p), a)
        Lambda <- rbind(diag(1, p), B)
        mu_z <- c(alpha + Lambda %*% mu_x)
        Sigma_z <- Lambda %*% Sigma_x %*% t(Lambda) + diag(phi, r, r)
        Sigma_z <- 0.5 * (Sigma_z + t(Sigma_z))
        Sigma_z_inv <- solve(Sigma_z, tol = 1e-30)
        Delta <- solve(Sigma_x, tol = 1e-30) + ((1 / phi) * crossprod(Lambda))
        Delta_inv <- solve(Delta, tol = 1e-30)
        A <- diag(p) - (1 / phi) * Delta_inv %*% crossprod(Lambda)
        lambda_vz <- lambda - r / 2
        for (i in 1:n) {
          delta_i <-  t(Z[i, ] - mu_z) %*% Sigma_z_inv %*% (Z[i, ] - mu_z)
          m1_i <-
            ghyp::Egig(
              lambda = lambda_vz,
              chi = delta_i + eta,
              psi = eta,
              func = c("1/x")
            )
          b_i <-  (1 / phi) * Delta_inv %*% t(Lambda) %*% (Z[i, ] - mu_z)
          m2_i <- m1_i * (mu_x + b_i)
          m3_i <- Delta_inv + m1_i * tcrossprod(mu_x + b_i)
          m3_i <- 0.5 * (m3_i + t(m3_i)) # porque es simetrica la matriz
          sum_m1 <- sum_m1 + m1_i
          sum_m2 <- sum_m2 + m2_i
          sum_m3 <- sum_m3 + m3_i
          s1 <- s1 + (m1_i * Y[i, ])  
          s2 <- s2 + (Y[i, ] %*% t(m2_i))
          s3 <- s3 + c(
            m1_i * (crossprod(X[i, ]) + crossprod(Y[i, ] - a)) -
              crossprod(X[i, ], m2_i) -
              t(Y[i, ] - a) %*% B %*% m2_i -
              crossprod(m2_i, X[i, ]) -
              t(m2_i) %*% t(B) %*% (Y[i, ] - a) +
              sum(diag(m3_i)) +
              sum(diag(m3_i %*% t(B) %*% B))
          )
          s4 <- s4 + m3_i -
            tcrossprod(m2_i, mu_x) -
            tcrossprod(mu_x, m2_i) +
            m1_i * tcrossprod(mu_x)
        }
        # parameters update
        B <- (s2 - (s1 / sum_m1) %*% t(sum_m2)) %*% solve((sum_m3 - (sum_m2 / sum_m1) %*% t(sum_m2)), tol = 1e-30)
        a <- (s1 / sum_m1) - (B %*% (sum_m2 / sum_m1))
        phi <- (1 / (n * r)) * s3
        mu_x <- sum_m2 / c(sum_m1)
        Sigma_x <- (1 / n) * s4
        
        
        theta <-
          c(a,
            B,
            phi,
            mu_x,
            matrixcalc::vech(Sigma_x),
            lambda,
            eta)
        
        #criterio
        lk <- lLz(
          X,
          Y,
          theta,
          type = type
        )$logL
        loglik[count] <- lk
        
        if (StopCrit == "logL") {
          lk1 <- lk
          criterio <- abs(lk1 / lk0 - 1)
          lk0 <- lk1
          if (count == MaxIter) {
            criterio <- 1e-11
          }
          # print(paste("criterio de convergencia simple:", criterio))
          # print(paste("loglikelihood", lk))
          # print(c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x), eta))
        }
        if (StopCrit == "Aitken") {
          if (count < 2) {
            lk1 <- lk
            criterio <- abs(lk1 / lk0 - 1)
          }
          if (count == 2) {
            lk2 <- lk
            c1 <- (lk2 - lk1) / (lk1 - lk0)
            l_A2 <- lk1 + ((lk2 - lk1) / (1 - c1))
            criterio <- abs(lk2 / lk1 - 1)
          }
          if (count > 2) {
            lk3 <- lk
            c <- (lk3 - lk2) / (lk2 - lk1)
            l_A3 <- lk2 + ((lk3 - lk2) / (1 - c))
            criterio <- abs(l_A3 - l_A2)
            l_A2 <- l_A3
            lk1 <- lk2
            lk2 <- lk3
          }
          if (count == MaxIter) {
            criterio <- 1e-11
          }
          # print(paste("criterio de convergencia Aitken:", criterio))
          # print(paste("loglikelihood", lk))
          # print(c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x), eta))
          
        }
        crits[count] <- criterio
        # 
        # theta <-
        #   c(a,
        #     B,
        #     phi,
        #     mu_x,
        #     matrixcalc::vech(Sigma_x),
        #     lambda,
        #     eta)
      }
    }
    if (type == "S-GH") {
      # S-GH ----
      lambda <- lambda
      gamma_x <- rep(0, p)
      eta <- eta
      theta <-c(a,
                 B,
                 phi,
                 mu_x,
                 matrixcalc::vech(Sigma_x),
                 lambda,
                 eta)
      lk0 <-
        lLz(
          X,
          Y,
          theta,
          type = type,
        )$logL
      loglik <- NULL
      criterio <- 1
      crits <- NULL
      count <- 0
      while (criterio > precision) {
        sum_m1 <- 0
        sum_m2 <- matrix(0, nrow = p, ncol = 1)
        sum_m3 <- matrix(0, nrow = p, ncol = p)
        mi2 <- matrix(0, nrow = p, ncol = 1)
        mi3 <- matrix(0, nrow = p, ncol = p)
        s1 <- rep(0, q) #\sum_{i=1}^{n}m_{i1}\bY_i
        s2 <- matrix(0, nrow = q, ncol = p)   #\sum_{i=1}^{n}\bY_i \bm_{i2}^{\top} 
        s3 <- 0 # \sum_{i=1}^{n} \Big\{ m_{i1}\Big[\bX_i^{\top} \bX_i + 
        # (\bY_i-\hat{\ba})^{\top}(\bY_i-\hat{\ba})\Big] - \bX_i^{\top}\bm_{i2} 
        #- (\bY_i - \hat{\ba})^{\top}\hat{\bB} \bm_{i2} - \bm_{i2}^{\top}\bX_i 
        #- \bm_{i2}^{\top}\hat{\bB}^{\top}(\bY-\hat{\ba}) + \text{tr}(\bm_{i3}) 
        # + \text{tr}(\bm_{i3}\hat{\bB}^{\top}\hat{\bB})\Big\}
        s4 <- matrix(0, nrow = p, ncol = p) # \sum_{i=1}^{n}\Big( \bm_{i3} - 
        # \bm_{i2}\hat{\bmu}_x^{\top} - 
        # \hat{\bmu}_x\bm_{i2}^{\top} + 
        # m_{i1}\hat{\bmu}_x\hat{\bmu}_x^{\top}\Big)
        count <- count + 1
        print(count)
        alpha <- c(rep(0, p), a)
        Lambda <- rbind(diag(1, p), B)
        mu_z <- c(alpha + Lambda %*% mu_x)
        Sigma_z <- Lambda %*% Sigma_x %*% t(Lambda) + diag(phi, r, r)
        Sigma_z <- 0.5 * (Sigma_z + t(Sigma_z))
        Sigma_z_inv <- solve(Sigma_z, tol = 1e-30)
        Delta <- solve(Sigma_x, tol = 1e-30) + ((1 / phi) * crossprod(Lambda))
        Delta_inv <- solve(Delta, tol = 1e-30)
        A <- diag(p) - (1 / phi) * Delta_inv %*% crossprod(Lambda)
        lambda_vz <- lambda - r / 2
        for (i in 1:n) {
          delta_i <-  t(Z[i, ] - mu_z) %*% Sigma_z_inv %*% (Z[i, ] - mu_z)
          m1_i <-
            ghyp::Egig(
              lambda = lambda_vz,
              chi = delta_i + eta,
              psi = eta,
              func = c("1/x")
            )
          b_i <-  (1 / phi) * Delta_inv %*% t(Lambda) %*% (Z[i, ] - mu_z)
          m2_i <- m1_i * (mu_x + b_i)
          m3_i <- Delta_inv + m1_i * tcrossprod(mu_x + b_i)
          m3_i <- 0.5 * (m3_i + t(m3_i)) # porque es simetrica la matriz
          sum_m1 <- sum_m1 + m1_i
          sum_m2 <- sum_m2 + m2_i
          sum_m3 <- sum_m3 + m3_i
          s1 <- s1 + (m1_i * Y[i, ])  
          s2 <- s2 + (Y[i, ] %*% t(m2_i))
          s3 <- s3 + c(
            m1_i * (crossprod(X[i, ]) + crossprod(Y[i, ] - a)) -
              crossprod(X[i, ], m2_i) -
              t(Y[i, ] - a) %*% B %*% m2_i -
              crossprod(m2_i, X[i, ]) -
              t(m2_i) %*% t(B) %*% (Y[i, ] - a) +
              sum(diag(m3_i)) +
              sum(diag(m3_i %*% t(B) %*% B))
          )
          s4 <- s4 + m3_i -
            tcrossprod(m2_i, mu_x) -
            tcrossprod(mu_x, m2_i) +
            m1_i * tcrossprod(mu_x)
        }
        # parameters update
        B <- (s2 - (s1 / sum_m1) %*% t(sum_m2)) %*% solve((sum_m3 - (sum_m2 / sum_m1) %*% t(sum_m2)), tol = 1e-30)
        a <- (s1 / sum_m1) - (B %*% (sum_m2 / sum_m1))
        phi <- (1 / (n * r)) * s3
        mu_x <- sum_m2 / c(sum_m1)
        Sigma_x <- (1 / n) * s4
        
        theta <-
          c(a,
            B,
            phi,
            mu_x,
            matrixcalc::vech(Sigma_x),
            lambda,
            eta)
        
        #criterio
        lk <- lLz(
          X,
          Y,
          theta,
          type = type
        )$logL
        loglik[count] <- lk
        
        if (StopCrit == "logL") {
          lk1 <- lk
          criterio <- abs(lk1 / lk0 - 1)
          lk0 <- lk1
          if (count == MaxIter) {
            criterio <- 1e-11
          }
          print(paste("criterio de convergencia simple:", criterio))
          print(paste("loglikelihood", lk))
          print(c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x), eta))
        }
        if (StopCrit == "Aitken") {
          if (count < 2) {
            lk1 <- lk
            criterio <- abs(lk1 / lk0 - 1)
          }
          if (count == 2) {
            lk2 <- lk
            c1 <- (lk2 - lk1) / (lk1 - lk0)
            l_A2 <- lk1 + ((lk2 - lk1) / (1 - c1))
            criterio <- abs(lk2 / lk1 - 1)
          }
          if (count > 2) {
            lk3 <- lk
            c <- (lk3 - lk2) / (lk2 - lk1)
            l_A3 <- lk2 + ((lk3 - lk2) / (1 - c))
            criterio <- abs(l_A3 - l_A2)
            l_A2 <- l_A3
            lk1 <- lk2
            lk2 <- lk3
          }
          if (count == MaxIter) {
            criterio <- 1e-11
          }
          print(paste("criterio de convergencia Aitken:", criterio))
          print(paste("loglikelihood", lk))
          print(c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x), eta))
          
        }
        crits[count] <- criterio
        
        # theta <-
        #   c(a,
        #     B,
        #     phi,
        #     mu_x,
        #     matrixcalc::vech(Sigma_x),
        #     lambda,
        #     eta)
      }
    }

    
      names(theta) <- c(nam(p,q),"lambda", "eta")
      resultados <- list(
        a = a,
        B = B,
        phi = phi,
        mu.x = mu_x,
        Sigma.x = Sigma_x,
        lambda = lambda,
        eta = eta,
        theta = theta,
        logLz = lk,
        logs = loglik, # despues se puede quitar
        iter = count,
        crits = crits, # despues se puede quitar
        aic = 2 * d - 2 * lk
      )
      resultados

  }

# en esta función se ajusta el modelo para una grilla de etas usando
# fit_SGH_MRMME() y proporciona eta estimado por verosimilitud profile.

fit_SGH_MRMME_profile <- function(X,Y,etas=seq(0.01,5,length = 20), type = "S-NIG"){
  # ajuste del modelo para cada eta
  logliks <- data.frame()
  for (eta in etas){
    fit <- fit_SGH_MRMME(X = X, Y = Y,
                         eta = eta, type = type)
    temp_df <- data.frame(eta = eta, logL = fit$logLz)
  logliks <- rbind(logliks, temp_df)
  }

  # grafico de la verosimilitud profile
  plot_prof <- ggplot() + 
  geom_line(data = logliks ,aes(x = eta, y = logL)) +
  geom_point(data = logliks ,aes(x = eta, y = logL)) +
  labs(x = expression(eta)) +
  theme_bw()

  # eta optimo
  hat_eta <- logliks %>% arrange(desc(logL)) %>% slice(1)

  # ajuste con el eta optimo
  fit_gh <- fit_SGH_MRMME(X = X, Y = Y,
                          eta = hat_eta$eta, 
                          type = type)
  lLz_i <- lLz(X, Y, 
               fit_gh$theta,
               type = type)$logL_i
  

  resultados <- list(
    a = fit_gh$a,
    B = fit_gh$B,
    phi = fit_gh$phi,
    mu.x = fit_gh$mu.x,
    Sigma.x = fit_gh$Sigma.x,
    lambda = fit_gh$lambda,
    eta = fit_gh$eta,
    theta = fit_gh$theta,
    logLz = fit_gh$logLz,
    logLz_i = lLz_i,
    iter = fit_gh$iter,
    plot_profile = plot_prof,
    logliks = logliks,
    aic = fit_gh$aic
  )
  resultados
}



# Assessment --------------------------------------------------------------

# funcion para evaluar el ajuste del modelo

model_assess <- function(X,Y,logLz_i, 
                         theta,
                         n_sample, type){
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Z <- cbind(X, Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2
  
  # theta --> (a,B,phi,mu,Sigma)
  theta <- c(theta)
  a <- matrix(theta[1:q], ncol = 1, nrow = q)
  B <- matrix(theta[(q + 1):(q + p * q)], nrow = q, ncol = p)
  phi <- theta[q + p * q + 1]
  mu.x <- theta[(q + p * q + 2):(q + p * q + 1 + p)]
  Sigma.x <- ks::invvech(theta[(q + p * q + 2 + p):d])
  Sigma.x <- 0.5 * (Sigma.x + t(Sigma.x))
  lambda <- theta[(d + 1)]
  eta <- theta[(d + 2)]
  
  
  assess <- data.frame(ll = logLz_i, id = 1)
  for(i in 2:n_sample){
    data <- sim_SGH_data(n = length(logLz_i), p = p, q = q,
                         a =  a, B = B, phi = phi, 
                         mu_x = mu.x, Sigma_x = Sigma.x, lambda = lambda, eta = eta,
                         type = type)
    logLz_i_fit <- lLz(data$X, data$Y, 
                       theta,
                       type = type)$logL_i
    assess <- rbind(assess, data.frame(ll = logLz_i_fit, id = i))
  }
  
  
  # grafico con los límites de las ecdf
  
  # calcula las ecdf:
  ecdf <- assess  %>% 
    group_by(id) %>% 
    mutate(ecdf = ecdf(ll)(unique(ll))) %>% 
    mutate(ll2 = round(ll,1)) %>% ungroup() 
  
  # calcula los lim inferior y superior de las ecdf calculadas con 
  # los datos simulados
  
  lims_ecdf <- ecdf %>% 
    filter(! id %in% c(1)) %>% # se quita el primero porque es la ecdf del modelo
    group_by(ll2) %>% 
    summarise(ymin = min(ecdf),
              ymax = max(ecdf)) %>% 
    ungroup()
  
  # plot1_SGH <- ggplot() +
  #   geom_ribbon(data = lims_ecdf, aes(x = ll2, ymax = ymax, ymin = ymin), 
  #               linewidth = 0.1, fill = "blue", alpha = 0.5) +
  #   geom_line(data = ecdf %>% filter(id == 1), aes(x = ll, y = ecdf),
  #             linewidth = 0.6, color = "darkblue") + 
  #   geom_point(data = ecdf %>% filter(id == 1), aes(x = ll, y = ecdf),
  #              color = "darkblue") +
  #   theme_bw() +
  #   xlab("Log density")
  
  # plot con las lineas de las ecdf
  
  plot2_SGH <- ggplot(assess, aes(x = ll, fill = as.factor(id))) +
    stat_ecdf(color = "darkgray") + theme_bw() +
    theme(legend.position = "none") + 
    geom_line(data = ecdf %>% filter(id == 1), aes(x = ll, y = ecdf),
              linewidth = 0.6, color = "darkblue") + 
    geom_point(data = ecdf %>% filter(id == 1), aes(x = ll, y = ecdf),
               color = "darkblue", size = 0.8) +
    xlab("Log density")
  
  #plot1_SGH | plot2_SGH
  plot2_SGH
  
}

# model assessment modelo normal 

model_assess_N <- function(X, Y, 
                           theta,
                           n_sample){
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Z <- cbind(X, Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2
  
  logLz_i = logLi_N(theta = theta, X = X, Y=Y)
  
  # theta --> (a,B,phi,mu,Sigma)
  theta <- c(theta)
  a <- matrix(theta[1:q], ncol = 1, nrow = q)
  B <- matrix(theta[(q + 1):(q + p * q)], nrow = q, ncol = p)
  phi <- theta[q + p * q + 1]
  mu.x <- theta[(q + p * q + 2):(q + p * q + 1 + p)]
  Sigma.x <- ks::invvech(theta[(q + p * q + 2 + p):d])
  Sigma.x <- 0.5 * (Sigma.x + t(Sigma.x))
  #lambda <- theta[(d + 1)]
  #eta <- theta[(d + 2)]
  
  assess_N <- data.frame(ll = logLz_i, id = 1)
  for(i in 2:n_sample){
    data <- sim_Ndata(n = length(logLz_i), p = p, q = q,
                      a =  a, B = B, phi = phi, 
                      mu_x = mu.x, Sigma_x = Sigma.x)
    # logLz_i_fit <- lLz(data$X, data$Y, 
    #                    a, B, phi, mu.x, Sigma.x, 
    #                    lambda = lambda, eta = eta, type = type)$logL_i
    logLz_i_fit <- logLi_N(theta = theta, X = data$X, Y = data$Y)
    assess_N <- rbind(assess_N, data.frame(ll = logLz_i_fit, id = i))
  }
  
  
  # # grafico con los límites de las ecdf
  ecdf_N <- assess_N  %>%
    group_by(id) %>%
    mutate(ecdf = ecdf(ll)(unique(ll))) %>%
    mutate(ll2 = round(ll,1)) %>% ungroup()
  # 
  # lims_ecdf_N <- ecdf_N %>% 
  #   filter(! id %in% c(1)) %>% # se quita el primero porque es la ecdf del modelo
  #   group_by(ll2) %>% 
  #   summarise(ymin = min(ecdf),
  #             ymax = max(ecdf)) %>% 
  #   ungroup()
  
  # plot1_N <- ggplot() +
  #   geom_ribbon(data = lims_ecdf_N, aes(x = ll2, ymax = ymax, ymin = ymin), 
  #               linewidth = 0.1, fill = "blue", alpha = 0.5) +
  #   geom_line(data = ecdf_N %>% filter(id == 1), aes(x = ll, y = ecdf),
  #             linewidth = 0.6, color = "darkblue") + 
  #   geom_point(data = ecdf_N %>% filter(id == 1), aes(x = ll, y = ecdf),
  #              color = "darkblue") +
  #   theme_bw() +
  #   xlab("Log density")
  # plot con las lineas de las ecdf
  
  plot2_N <- ggplot(assess_N, aes(x = ll, fill = as.factor(id))) +
    stat_ecdf(geom = "step", color = "darkgray") + theme_bw() +
    theme(legend.position = "none") + 
    geom_line(data = ecdf_N %>% filter(id == 1), aes(x = ll, y = ecdf),
              linewidth = 0.6, color = "darkblue") + 
    geom_point(data = ecdf_N %>% filter(id == 1), aes(x = ll, y = ecdf),
               color = "darkblue", size = 0.8) +
    xlab("Log density")
  # plot1_N | plot2_N
  plot2_N
}




# Score ----
score_matrix_SGH <- function(theta, X, Y, total = FALSE) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Z <- cbind(X, Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2
  
  # theta --> (a,B,phi,mu,Sigma)
  theta <- c(theta)
  a <- matrix(theta[1:q], ncol = 1, nrow = q)
  B <- matrix(theta[(q + 1):(q + p * q)], nrow = q, ncol = p)
  phi <- theta[q + p * q + 1]
  mu_x <- theta[(q + p * q + 2):(q + p * q + 1 + p)]
  Sigma_x <- ks::invvech(theta[(q + p * q + 2 + p):d])
  Sigma_x <- 0.5 * (Sigma_x + t(Sigma_x))
  lambda <- theta[(d + 1)]
  eta <- theta[(d + 2)]
  
  alpha <- matrix(c(rep(0, p), a), ncol = 1)
  Lambda <- rbind(diag(1, p), B)
  Sigma_e <- diag(phi, r, r)
  mu_z <- c(alpha + Lambda %*% mu_x)
  Sigma_z <- Lambda %*% tcrossprod(Sigma_x, Lambda) + Sigma_e
  Sigma_z <- 0.5 * (Sigma_z + t(Sigma_z))
  tau <- t(t(Z) - mu_z)
  delta <- diag(tau %*% solve(Sigma_z, tol = 1e-30) %*% t(tau))
  
  #' (v_i | z_i) ind GIG(lambda_vz, chi_vz_{i}, psi_vz)
  lambda_vz <- lambda - (r / 2)
  chi_vz <- delta + eta
  psi_vz <- eta
  omega <- sqrt(chi_vz * psi_vz)
  
  # weights ----
  #' mi1, mi2, mi3 ----
  m1 <-
    ghyp::Egig(
      lambda = lambda_vz,
      chi = chi_vz,
      psi = psi_vz,
      func = c("1/x")
    )
  
  Delta <- solve(Sigma_x, tol = 1e-30) + ((1 / phi) * crossprod(Lambda))
  Delta <- 0.5 * (Delta + t(Delta))
  invDelta <- solve(Delta, tol = 1e-30)
  b <- t((1 / phi) * invDelta %*% t(Lambda) %*% t(tau))
  
  m2 <- (m1 * t(mu_x + t(b)))
  
  m3 <- list(NULL)
  mu_xb <- t(mu_x + t(b))
  for (i in 1:n) {
    aux1m3 <- m1[i] * tcrossprod(mu_xb[i,])
    m3[[i]] <- invDelta + aux1m3
    m3[[i]] <- 0.5 * (m3[[i]] + t(m3[[i]]))
  }
  
  score_matrix <- matrix(ncol = d, nrow = n)
  for (i in 1:n) {
    # a ----
    score_a <- c( (1 / phi) * ((Y[i,] - a) * m1[i] - B %*% m2[i,]))
    # B ----
    score_B <- c((1 / phi) * (tcrossprod((Y[i,] - a),m2[i,]) - (B %*% m3[[i]])))
    
    # phi ----
    score_phi <-
      c(
        -r / (2 * phi) + 1 / (2 * phi ^ 2) * (
          (crossprod(X[i,]) + crossprod(Y[i,] - a)) * m1[i] -
            crossprod(X[i,], m2[i,]) -
            t(Y[i,] - a) %*% B %*% m2[i,] -
            crossprod(m2[i,], X[i,]) -
            t(m2[i,]) %*% t(B) %*% (Y[i,] - a) +
            sum(diag(m3[[i]])) +
            sum(diag(m3[[i]] %*% t(B) %*% B))
        )
      )
    # mu_x ----
    score_mu_x <-
      c((solve(Sigma_x) %*% m2[i,]) - (solve(Sigma_x) %*% matrix(m1[i] * mu_x, ncol = 1)))
    
    # Sigma_x ----
    score_Sigma_x <-
      c(matrixcalc::vech(-1 / 2 * solve(Sigma_x) + 1 / 2 * (
        solve(Sigma_x) %*%
          (m3[[i]] - tcrossprod(m2[i,], mu_x) -
           tcrossprod(mu_x, m2[i,]) +
           m1[i]* tcrossprod(mu_x)
          ) %*% 
        solve(Sigma_x)
      )))
      
    score_matrix[i,] <-
      c(score_a,
        score_B,
        score_phi,
        score_mu_x,
        score_Sigma_x)
  }
  if (total == TRUE) {
    score_t <- colSums(score_matrix)
    return(score_t)
  } else{
    return(score_matrix)
  }
}

# funcion de score en el modelos normal para evaluar theta tilde en test de hipotesis

score_matrix_N <- function(theta, X, Y, total = FALSE) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Z <- cbind(X, Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2
  
  # theta --> (a,B,phi,mu,Sigma)
  theta <- c(theta)
  a <- matrix(theta[1:q], ncol = 1, nrow = q)
  B <- matrix(theta[(q + 1):(q + p * q)], nrow = q, ncol = p)
  phi <- theta[q + p * q + 1]
  mu_x <- theta[(q + p * q + 2):(q + p * q + 1 + p)]
  Sigma_x <- ks::invvech(theta[(q + p * q + 2 + p):d])
  Sigma_x <- 0.5 * (Sigma_x + t(Sigma_x))
  
  # eta and psi
  Lambda <- rbind(diag(1, p, p), B)
  alpha <- c(rep(0, p), a)
  eta <- c(alpha + Lambda %*% mu_x)
  psi <- Lambda %*% tcrossprod(Sigma_x, Lambda) + diag(phi, r, r)
  psi <- 0.5 * (psi + t(psi))
  psi_inv <- Rfast::spdinv(psi)
  
  # Score function
  score <- function(psi_j, partial_delta_j) {
    score =  -0.5 * (matrixcalc::matrix.trace(psi_inv %*% psi_j) + partial_delta_j)
    score
  }
  
  # jth - partial derivative
  partial_delta_j <- function(psi_j_inv, eta_j) {
    pdj <- stats::mahalanobis(Z, eta, psi_j_inv, inverted = TRUE) -
      2 * (t(t(Z) - eta) %*% psi_inv %*% eta_j)
    pdj
  }
  
  # score a----
  score_a <- matrix(, nrow = n, ncol = q)
  index <- NULL
  for (i in 1:q) {
    psi_j_a <- matrix(0, nrow = r, ncol = r)
    psi_j_inv_a <- -psi_inv %*% psi_j_a %*% psi_inv
    eta_j_a <- matrix(c(rep(0, p), diag(q)[, i]), ncol = 1)
    partial_delta_j_a <- partial_delta_j(psi_j_inv_a, eta_j_a)
    score_a[, i] <- score(psi_j_a, partial_delta_j_a)
  }
  
  # score vec(B)----
  # E_mat: matrix with 1 in k,l position
  E_mat <- list()
  for (l in 1:p) {
    E_mat[[l]] <- list()
    for (k in 1:q) {
      E_aux =  matrix(0, nrow = q, ncol = p)
      E_aux[k, l] =  1
      E_mat[[l]][[k]] <- E_aux
    }
  }
  E_mat <- unlist(E_mat, recursive = FALSE)
  
  eta_j_b <- list()
  for (k in 1:(p * q)) {
    aux <- E_mat[[k]] %*% mu_x
    eta_j_b[[k]] <- rbind(matrix(rep(0, p), ncol = 1), aux)
  }
  
  
  psi_j_b <- list()
  for (k in 1:(p * q)) {
    psi_j_b[[k]] <- list()
    row1 = cbind(matrix(0, nrow = p, ncol = p), tcrossprod(Sigma_x, E_mat[[k]]))
    row2 = cbind(E_mat[[k]] %*% Sigma_x,
                 B %*% tcrossprod(Sigma_x, E_mat[[k]]) +
                   E_mat[[k]] %*% Sigma_x %*% t(B))
    psi_j_b[[k]] <- rbind(row1, row2)
  }
  
  psi_j_inv_b <- list()
  for (k in 1:(p * q)) {
    psi_j_inv_b[[k]] <- -psi_inv %*% psi_j_b[[k]] %*% psi_inv
  }
  
  partial_delta_j_b <- list()
  for (k in 1:(p * q)) {
    partial_delta_j_b[[k]] <-
      partial_delta_j(psi_j_inv_b[[k]], eta_j_b[[k]])
  }
  
  score_b <- matrix(, nrow = n, ncol = (q * p))
  for (k in 1:(q * p)) {
    score_b[, k] <- score(psi_j_b[[k]], partial_delta_j_b[[k]])
  }
  
  # score phi----
  eta_j_phi <- matrix(rep(0, r), ncol = 1)
  psi_j_phi <- diag(r)
  psi_j_inv_phi <- -psi_inv %*% psi_j_phi %*% psi_inv
  partial_delta_j_phi <- partial_delta_j(psi_j_inv_phi, eta_j_phi)
  score_phi <- score(psi_j_phi, partial_delta_j_phi)
  
  
  # score mu----
  score_mu <- matrix(, nrow = n, ncol = p)
  for (i in 1:p) {
    eta_j_mu <- matrix(c(diag(p)[, i], c(B %*% diag(p)[, i])), ncol = 1)
    psi_j_mu <- matrix(0, nrow = r, ncol = r)
    psi_j_inv_mu <- -psi_inv %*% psi_j_mu %*% psi_inv
    partial_delta_j_mu <- partial_delta_j(psi_j_inv_mu, eta_j_mu)
    score_mu[, i] <- score(psi_j_mu, partial_delta_j_mu)
  }
  
  # score vech(Sigma)----
  E_mat_Sigma <- list()
  for (l in 1:p) {
    E_mat_Sigma[[l]] <- list()
    for (k in 1:p) {
      if (k >= l) {
        E_aux =  matrix(0, nrow = p, ncol = p)
        E_aux[k, l] =  1
        E_aux[l, k] =  1
        E_mat_Sigma[[l]][[k]] <- E_aux
      }
    }
  }
  
  eta_j_Sigma <- list()
  for (l in 1:p) {
    eta_j_Sigma[[l]] <- list()
    for (k in 1:p) {
      if (k >= l) {
        eta_j_Sigma[[l]][[k]] <- matrix(rep(0, r), ncol = 1)
      }
    }
  }
  
  psi_j_Sigma <- list()
  for (l in 1:p) {
    psi_j_Sigma[[l]] <- list()
    for (k in 1:p) {
      if (k >= l) {
        psi_j_Sigma[[l]][[k]] <-
          Lambda %*% E_mat_Sigma[[l]][[k]] %*% t(Lambda)
      }
    }
  }
  
  psi_j_inv_Sigma <- list()
  for (l in 1:p) {
    psi_j_inv_Sigma[[l]] <- list()
    for (k in 1:p) {
      if (k >= l) {
        psi_j_inv_Sigma[[l]][[k]] <-
          -psi_inv %*% psi_j_Sigma[[l]][[k]] %*% psi_inv
      }
    }
  }
  
  partial_delta_j_Sigma <- list()
  for (l in 1:p) {
    partial_delta_j_Sigma[[l]] <- list()
    for (k in 1:p) {
      if (k >= l) {
        partial_delta_j_Sigma[[l]][[k]] <-
          partial_delta_j(psi_j_inv_Sigma[[l]][[k]], eta_j_Sigma[[l]][[k]])
      }
    }
  }
  
  score_Sigma <- list()
  for (l in 1:p) {
    score_Sigma[[l]] <- list()
    for (k in 1:p) {
      if (k >= l) {
        score_Sigma[[l]][[k]] <-
          score(psi_j_Sigma[[l]][[k]], partial_delta_j_Sigma[[l]][[k]])
      }
    }
  }
  score_Sigma <-
    matrix(unlist(score_Sigma),
           ncol = p * (p + 1) / 2,
           byrow = FALSE)
  
  # score matrix----
  score_theta <-
    cbind(score_a, score_b, score_phi, score_mu, score_Sigma)
  colnames(score_theta) <- nam(p, q)
  
  if (total == TRUE) {
    score_t <- colSums(score_theta)
    return(score_t)
  } else{
    return(score_theta)
  }
}

# Standar errors ----

se_SGH <- function(theta, X, Y){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2
  score_mat <- score_matrix_SGH(theta, X, Y)
  fim <- crossprod(score_mat) # empirical fisher information matrix 
  fim <- 0.5 * (fim + t(fim))
  colnames(fim) <- nam(p, q)
  rownames(fim) <- nam(p, q)
  D <-  Rfast::spdinv(fim)
  D <- 0.5 * (D + t(D))
  se <- sqrt(diag(D))
  colnames(D) <- nam(p, q)
  rownames(D) <- nam(p, q)
  names(se) <- nam(p, q)
  res <- list("covmat" = D,
              "se" = se)
  res
}  

fim_N <- function(theta, X, Y, type = "emp") {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2
  
  # expected ----
  if (type == "exp") {
    # theta --> (a,B,phi,mu,Sigma)
    theta <- as.vector(theta)
    a <- matrix(theta[1:q], ncol = 1, nrow = q)
    B <- matrix(theta[(q + 1):(q + p * q)], nrow = q, ncol = p)
    phi <- theta[q + p * q + 1]
    mu_x <- theta[(q + p * q + 2):(q + p * q + 1 + p)]
    Sigma_x <- ks::invvech(theta[(q + p * q + 2 + p):d])
    Sigma_x <- 0.5 * (Sigma_x + t(Sigma_x))
    
    # eta and psi
    Lambda <- rbind(diag(1, p, p), B)
    alpha <- c(rep(0, p), a)
    eta <- c(alpha + Lambda %*% mu_x)
    psi <- Lambda %*% tcrossprod(Sigma_x, Lambda) + diag(phi, r, r)
    psi <- 0.5 * (psi + t(psi))
    psi_inv <-  Rfast::spdinv(psi)
    
    # aux's
    V <-
      rbind(matrix(0, nrow = p, ncol = q), diag(nrow = q, ncol = q))
    A <- diag(nrow = r, ncol = r)
    
    # block a
    cov_a <- n * crossprod(V, psi_inv) %*% V
    
    # block B
    E_mat <- list()
    for (j in 1:(q * p)) {
      E_mat[[j]] <- list()
      for (l in 1:p) {
        E_mat[[j]][[l]] <- list()
        for (k in 1:q) {
          E_aux =  matrix(0, nrow = q, ncol = p)
          E_aux[k, l] =  1
          E_mat[[j]][[l]][[k]] <-
            rbind(diag(0, nrow = p, ncol = p), E_aux)
        }
      }
    }
    
    E_mat_aux <- unlist(E_mat[[1]], recursive = FALSE)
    
    cov_B <- matrix(0, nrow = q * p, ncol = q * p)
    for (i in 1:(q * p)) {
      izq <- E_mat_aux[[i]]
      fim_aux <- NULL
      for (k in 1:p) {
        for (l in 1:q) {
          der <- E_mat[[i]][[k]][[l]]###
          p1 <- crossprod(izq, psi_inv) %*% der
          s1 <- crossprod(mu_x, p1) %*% mu_x
          p2 <- Lambda %*% tcrossprod(Sigma_x, izq)
          p3 <- izq %*% tcrossprod(Sigma_x, Lambda)
          p4 <- Lambda %*% tcrossprod(Sigma_x, der)
          s2 <- psi_inv %*% (p2 + p3) %*% psi_inv %*% p4
          aux <- n * (s1 + sum(diag(s2)))
          fim_aux <- c(fim_aux, aux)
        }
      }
      cov_B[, i] <- fim_aux
    }
    
    # block phi
    cov_phi <- (n / 2) * sum(diag(psi_inv %*% A %*% psi_inv %*% A))
    
    # block mu
    cov_mu <- n * crossprod(Lambda, psi_inv) %*% Lambda
    
    # block Sigma
    E_mat_s <- list()
    for (j in 1:(p * p)) {
      E_mat_s[[j]] <- list()
      for (l in 1:p) {
        E_mat_s[[j]][[l]] <- list()
        for (k in 1:p) {
          if (k >= l) {
            E_aux =  matrix(0, nrow = p, ncol = p)
            E_aux[k, l] =  1
            E_aux[l, k] =  1
            E_mat_s[[j]][[l]][[k]] <- E_aux
          }
        }
      }
    }
    
    E_mat_aux_s <- unlist(E_mat_s[[1]], recursive = FALSE)
    E_mat_s_aux <- E_mat_aux_s[lengths(E_mat_aux_s) != 0]
    
    p_aux <- length(matrixcalc::vech(Sigma_x))
    cov_S <- matrix(0, nrow = p_aux, ncol = p_aux)
    for (i in 1:p_aux) {
      izq <- E_mat_s_aux[[i]]
      fim_aux <- NULL
      for (k in 1:p) {
        for (l in 1:p) {
          if (k >= l) {
            der <- E_mat_s[[i]][[l]][[k]]
            p1 <- crossprod(Lambda, psi_inv) %*% Lambda
            s1 <- izq %*% p1 %*% der %*% p1
            aux <- (n / 2) * sum(diag(s1))
            fim_aux <- c(fim_aux, aux)
          }
        }
      }
      cov_S[i,] <- fim_aux
    }
    
    # block a with b
    cov_a_B <- matrix(0, nrow = q, ncol = (q * p))
    for (i in 1:(q * p)) {
      cov_a_B[, i] <-
        n * crossprod(V, psi_inv) %*% E_mat_aux[[i]] %*% mu_x
    }
    
    # block a with phi
    cov_a_phi <- matrix(0, nrow = q, ncol = 1)
    
    # block a with mu
    cov_a_mu <- n * crossprod(V, psi_inv) %*% Lambda
    
    # block a with Sigma
    cov_a_S <- matrix(0, nrow = q, ncol = p_aux)
    
    # block B with phi
    cov_B_phi <- matrix(0, nrow = (q * p), ncol = 1)
    for (i in 1:(q * p)) {
      p1 <- Sigma_x %*% crossprod(E_mat_aux[[i]], psi_inv)
      s1 <- A %*% psi_inv %*% Lambda %*% p1
      cov_B_phi[i,] <- n * sum(diag(s1))
    }
    
    # block B with mu
    cov_B_mu <- matrix(0, nrow = (q * p), ncol = p)
    for (i in 1:(q * p)) {
      p1 <- crossprod(Lambda, psi_inv)
      p2 <- E_mat_aux[[i]] %*% mu_x
      cov_B_mu[i,] <- n * p1 %*% p2
    }
    
    # block B with Sigma
    cov_B_S <- matrix(0, nrow = (q * p), ncol = p_aux)
    for (i in 1:(q * p)) {
      cov_aux <- NULL
      for (j in 1:p_aux) {
        p1 <- E_mat_s_aux[[j]] %*% crossprod(Lambda, psi_inv)
        p2 <-
          Lambda %*% Sigma_x %*% crossprod(E_mat_aux[[i]], psi_inv) %*% Lambda
        aux <- n * sum(diag(p1 %*% p2))
        cov_aux <- c(cov_aux, aux)
      }
      cov_B_S[i, ] <- cov_aux
    }
    
    # block phi with mu
    cov_phi_mu <- matrix(0, nrow = 1, ncol = p)
    
    # block phi with Sigma
    cov_phi_S <- NULL
    for (j in 1:p_aux) {
      p1 <- psi_inv %*% Lambda
      p2 <- E_mat_s_aux[[j]] %*% crossprod(Lambda, psi_inv) %*% A
      cov_phi_S[j] <- (n / 2) * sum(diag(p1 %*% p2))
    }
    cov_phi_S <- matrix(cov_phi_S, nrow = 1)
    
    # block mu with Sigma
    cov_mu_S <- matrix(0, nrow = p, ncol = p_aux)
    
    # join together by cols
    a <- cbind(cov_a, cov_a_B, cov_a_phi, cov_a_mu, cov_a_S)
    B <- cbind(t(cov_a_B), cov_B, cov_B_phi, cov_B_mu, cov_B_S)
    phi <-
      c(t(cov_a_phi), t(cov_B_phi), cov_phi, cov_phi_mu, cov_phi_S)
    mu <-
      cbind(t(cov_a_mu), t(cov_B_mu), t(cov_phi_mu), cov_mu, cov_mu_S)
    S <-
      cbind(t(cov_a_S), t(cov_B_S), t(cov_phi_S), t(cov_mu_S), cov_S)
    
    # join together by rows
    fim_CM <- rbind(a, B, phi, mu, S)
    fim_CM <- 0.5 * (fim_CM + t(fim_CM)) # this is the FIM total
    colnames(fim_CM) <- nam(p, q)
    rownames(fim_CM) <- nam(p, q)
    fim_CM # this is the FIM total
    
  } else {
    # observed ----
    if (type == "obs") {
      fim_obs <- -numDeriv::jacobian(
        func = score_matrix_N,
        x = theta,
        X = X,
        Y = Y,
        total = TRUE
      )
      colnames(fim_obs) <- nam(p, q)
      rownames(fim_obs) <- nam(p, q)
      fim_obs <-  0.5 * (fim_obs + t(fim_obs))
      fim_obs # This is the FIM total
      
    } else {
      # empirical ----
      if (type == "emp") {
        score_mat <- score_matrix_N(theta, X, Y)
        fim_emp <- crossprod(score_mat)
        fim_emp <- 0.5 * (fim_emp + t(fim_emp))
        colnames(fim_emp) <- nam(p, q)
        rownames(fim_emp) <- nam(p, q)
        fim_emp # This is the FIM total
      }
    }
  }
}

se_N <- function(theta, X, Y, type = "emp") {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Z <- cbind(X, Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2
  type <- type
  
  fim <- fim_N(theta, X, Y, type)
  D <-  Rfast::spdinv(fim)
  D <- 0.5 * (D + t(D))
  se <- sqrt(diag(D))
  colnames(D) <- nam(p, q)
  rownames(D) <- nam(p, q)
  names(se) <- nam(p, q)
  res <- list("covmat" = D,
              "se" = se)
  res
}


# Hypothesis testing ------------------------------------------------------

# esta funcion ajusta el modelo para un eta dado y para el modelo bajo H0
fit_SGH_MRMME_R <-
  function(X,
           Y,
           # a,
           # B,
           # phi,
           # mu_x,
           # Sigma_x,
           lambda = -1/2,
           eta, # este es fijo y se escoje de una grilla
           type = "S-NIG",
           StopCrit = "logL",
           MaxIter = 500,
           precision = 1e-10) {
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    n <- dim(Y)[1]
    q <- dim(Y)[2]
    p <- dim(X)[2]
    r <- p + q
    Z <- cbind(X, Y)
    d <- (p + 1) * (2 + p + 2 * q) / 2
    
    
    
    fit_normal <- MRMME::mrmme(Y=Y, X=X) #ajuste normal para los valores iniciales
    # a <- fit_normal$a
    a <- rep(0,q) # H0
    B <- fit_normal$B
    phi <- fit_normal$phi
    mu_x <- fit_normal$mu.x
    Sigma_x <- fit_normal$Sigma.x
    
    
    # mod_lm <- stats::lm(Y ~ X)
    # a <- matrix(stats::coef(mod_lm)[1,],ncol=1) # lm
    # B <- t(stats::coef(mod_lm)[-1,]) # lm
    # phi <- 1
    # mu_x <- colMeans(X) # observed
    # Sigma_x <- stats::cov(X) # observed
    
    
    
    if (type == "S-NIG") {
      # S-NIG ----
      lambda <- -1/2
      gamma_x <- rep(0, p)
      eta <- eta
      theta <-c(a,
                 B,
                 phi,
                 mu_x,
                 matrixcalc::vech(Sigma_x),
                 lambda,
                 eta)
      lk0 <-
        lLz(
          X,
          Y,
          theta,
          type = type
        )$logL
      loglik <- NULL
      criterio <- 1
      crits <- NULL
      count <- 0
      while (criterio > precision) {
        sum_m1 <- 0
        sum_m2 <- matrix(0, nrow = p, ncol = 1)
        sum_m3 <- matrix(0, nrow = p, ncol = p)
        mi2 <- matrix(0, nrow = p, ncol = 1)
        mi3 <- matrix(0, nrow = p, ncol = p)
        s1 <- rep(0, q) #\sum_{i=1}^{n}m_{i1}\bY_i
        s2 <- matrix(0, nrow = q, ncol = p)   #\sum_{i=1}^{n}\bY_i \bm_{i2}^{\top} 
        s3 <- 0 # \sum_{i=1}^{n} \Big\{ m_{i1}\Big[\bX_i^{\top} \bX_i + 
        # (\bY_i-\hat{\ba})^{\top}(\bY_i-\hat{\ba})\Big] - \bX_i^{\top}\bm_{i2} 
        #- (\bY_i - \hat{\ba})^{\top}\hat{\bB} \bm_{i2} - \bm_{i2}^{\top}\bX_i 
        #- \bm_{i2}^{\top}\hat{\bB}^{\top}(\bY-\hat{\ba}) + \text{tr}(\bm_{i3}) 
        # + \text{tr}(\bm_{i3}\hat{\bB}^{\top}\hat{\bB})\Big\}
        s4 <- matrix(0, nrow = p, ncol = p) # \sum_{i=1}^{n}\Big( \bm_{i3} - 
        # \bm_{i2}\hat{\bmu}_x^{\top} - 
        # \hat{\bmu}_x\bm_{i2}^{\top} + 
        # m_{i1}\hat{\bmu}_x\hat{\bmu}_x^{\top}\Big)
        count <- count + 1
        #print(count)
        alpha <- c(rep(0, p), a)
        Lambda <- rbind(diag(1, p), B)
        mu_z <- c(alpha + Lambda %*% mu_x)
        Sigma_z <- Lambda %*% Sigma_x %*% t(Lambda) + diag(phi, r, r)
        Sigma_z <- 0.5 * (Sigma_z + t(Sigma_z))
        Sigma_z_inv <- solve(Sigma_z, tol = 1e-30)
        Delta <- solve(Sigma_x, tol = 1e-30) + ((1 / phi) * crossprod(Lambda))
        Delta_inv <- solve(Delta, tol = 1e-30)
        #A <- diag(p) - (1 / phi) * Delta_inv %*% crossprod(Lambda)
        lambda_vz <- lambda - r / 2
        for (i in 1:n) {
          delta_i <-  t(Z[i, ] - mu_z) %*% Sigma_z_inv %*% (Z[i, ] - mu_z)
          m1_i <-
            ghyp::Egig(
              lambda = lambda_vz,
              chi = delta_i + eta,
              psi = eta,
              func = c("1/x")
            )
          b_i <-  (1 / phi) * Delta_inv %*% t(Lambda) %*% (Z[i, ] - mu_z)
          m2_i <- m1_i * (mu_x + b_i)
          m3_i <- Delta_inv + m1_i * tcrossprod(mu_x + b_i)
          m3_i <- 0.5 * (m3_i + t(m3_i)) # porque es simetrica la matriz
          sum_m1 <- sum_m1 + m1_i
          sum_m2 <- sum_m2 + m2_i
          sum_m3 <- sum_m3 + m3_i
          s1 <- s1 + (m1_i * Y[i, ])  
          s2 <- s2 + (Y[i, ] %*% t(m2_i))
          s3 <- s3 + c(
            m1_i * (crossprod(X[i, ]) + crossprod(Y[i, ] - a)) -
              crossprod(X[i, ], m2_i) -
              t(Y[i, ] - a) %*% B %*% m2_i -
              crossprod(m2_i, X[i, ]) -
              t(m2_i) %*% t(B) %*% (Y[i, ] - a) +
              sum(diag(m3_i)) +
              sum(diag(m3_i %*% t(B) %*% B))
          )
          s4 <- s4 + m3_i -
            tcrossprod(m2_i, mu_x) -
            tcrossprod(mu_x, m2_i) +
            m1_i * tcrossprod(mu_x)
        }
        # parameters update
        # a <- (s1 / sum_m1) - (B %*% (sum_m2 / sum_m1))
        a <- rep(0,q) #H0
        B <- s2  %*% solve(sum_m3, tol = 1e-30) # este estimador cambia por la hipotesis H0
        phi <- (1 / (n * r)) * s3
        mu_x <- sum_m2 / c(sum_m1)
        Sigma_x <- (1 / n) * s4
        
        theta <-
          c(a,
            B,
            phi,
            mu_x,
            matrixcalc::vech(Sigma_x),
            lambda,
            eta)
        
        #criterio
        lk <- lLz(
          X,
          Y,
          theta,
          type = type
        )$logL
        loglik[count] <- lk
        
        if (StopCrit == "logL") {
          lk1 <- lk
          criterio <- abs(lk1 / lk0 - 1)
          lk0 <- lk1
          if (count == MaxIter) {
            criterio <- 1e-11
          }
          # print(paste("criterio de convergencia simple:", criterio))
          # print(paste("loglikelihood", lk))
          # print(c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x), eta))
        }
        if (StopCrit == "Aitken") {
          if (count < 2) {
            lk1 <- lk
            criterio <- abs(lk1 / lk0 - 1)
          }
          if (count == 2) {
            lk2 <- lk
            c1 <- (lk2 - lk1) / (lk1 - lk0)
            l_A2 <- lk1 + ((lk2 - lk1) / (1 - c1))
            criterio <- abs(lk2 / lk1 - 1)
          }
          if (count > 2) {
            lk3 <- lk
            c <- (lk3 - lk2) / (lk2 - lk1)
            l_A3 <- lk2 + ((lk3 - lk2) / (1 - c))
            criterio <- abs(l_A3 - l_A2)
            l_A2 <- l_A3
            lk1 <- lk2
            lk2 <- lk3
          }
          if (count == MaxIter) {
            criterio <- 1e-11
          }
          print(paste("criterio de convergencia Aitken:", criterio))
          print(paste("loglikelihood", lk))
          print(c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x), eta))
          
        }
        crits[count] <- criterio
        
        # theta <-
        #   c(a,
        #     B,
        #     phi,
        #     mu_x,
        #     matrixcalc::vech(Sigma_x),
        #     lambda,
        #     eta)
      }
    }
    if (type == "S-HYP") {
      # S-HYP ----
      lambda <- 1
      gamma_x <- rep(0, p)
      eta <- eta
      theta <-c(a,
                 B,
                 phi,
                 mu_x,
                 matrixcalc::vech(Sigma_x),
                 lambda,
                 eta)
      lk0 <-
        lLz(
          X,
          Y,
          theta,
          type = type
        )$logL
      loglik <- NULL
      criterio <- 1
      crits <- NULL
      count <- 0
      while (criterio > precision) {
        sum_m1 <- 0
        sum_m2 <- matrix(0, nrow = p, ncol = 1)
        sum_m3 <- matrix(0, nrow = p, ncol = p)
        mi2 <- matrix(0, nrow = p, ncol = 1)
        mi3 <- matrix(0, nrow = p, ncol = p)
        s1 <- rep(0, q) #\sum_{i=1}^{n}m_{i1}\bY_i
        s2 <- matrix(0, nrow = q, ncol = p)   #\sum_{i=1}^{n}\bY_i \bm_{i2}^{\top} 
        s3 <- 0 # \sum_{i=1}^{n} \Big\{ m_{i1}\Big[\bX_i^{\top} \bX_i + 
        # (\bY_i-\hat{\ba})^{\top}(\bY_i-\hat{\ba})\Big] - \bX_i^{\top}\bm_{i2} 
        #- (\bY_i - \hat{\ba})^{\top}\hat{\bB} \bm_{i2} - \bm_{i2}^{\top}\bX_i 
        #- \bm_{i2}^{\top}\hat{\bB}^{\top}(\bY-\hat{\ba}) + \text{tr}(\bm_{i3}) 
        # + \text{tr}(\bm_{i3}\hat{\bB}^{\top}\hat{\bB})\Big\}
        s4 <- matrix(0, nrow = p, ncol = p) # \sum_{i=1}^{n}\Big( \bm_{i3} - 
        # \bm_{i2}\hat{\bmu}_x^{\top} - 
        # \hat{\bmu}_x\bm_{i2}^{\top} + 
        # m_{i1}\hat{\bmu}_x\hat{\bmu}_x^{\top}\Big)
        count <- count + 1
        print(count)
        alpha <- c(rep(0, p), a)
        Lambda <- rbind(diag(1, p), B)
        mu_z <- c(alpha + Lambda %*% mu_x)
        Sigma_z <- Lambda %*% Sigma_x %*% t(Lambda) + diag(phi, r, r)
        Sigma_z <- 0.5 * (Sigma_z + t(Sigma_z))
        Sigma_z_inv <- solve(Sigma_z, tol = 1e-30)
        Delta <- solve(Sigma_x, tol = 1e-30) + ((1 / phi) * crossprod(Lambda))
        Delta_inv <- solve(Delta, tol = 1e-30)
        A <- diag(p) - (1 / phi) * Delta_inv %*% crossprod(Lambda)
        lambda_vz <- lambda - r / 2
        for (i in 1:n) {
          delta_i <-  t(Z[i, ] - mu_z) %*% Sigma_z_inv %*% (Z[i, ] - mu_z)
          m1_i <-
            ghyp::Egig(
              lambda = lambda_vz,
              chi = delta_i + eta,
              psi = eta,
              func = c("1/x")
            )
          b_i <-  (1 / phi) * Delta_inv %*% t(Lambda) %*% (Z[i, ] - mu_z)
          m2_i <- m1_i * (mu_x + b_i)
          m3_i <- Delta_inv + m1_i * tcrossprod(mu_x + b_i)
          m3_i <- 0.5 * (m3_i + t(m3_i)) # porque es simetrica la matriz
          sum_m1 <- sum_m1 + m1_i
          sum_m2 <- sum_m2 + m2_i
          sum_m3 <- sum_m3 + m3_i
          s1 <- s1 + (m1_i * Y[i, ])  
          s2 <- s2 + (Y[i, ] %*% t(m2_i))
          s3 <- s3 + c(
            m1_i * (crossprod(X[i, ]) + crossprod(Y[i, ] - a)) -
              crossprod(X[i, ], m2_i) -
              t(Y[i, ] - a) %*% B %*% m2_i -
              crossprod(m2_i, X[i, ]) -
              t(m2_i) %*% t(B) %*% (Y[i, ] - a) +
              sum(diag(m3_i)) +
              sum(diag(m3_i %*% t(B) %*% B))
          )
          s4 <- s4 + m3_i -
            tcrossprod(m2_i, mu_x) -
            tcrossprod(mu_x, m2_i) +
            m1_i * tcrossprod(mu_x)
        }
        # parameters update
        # a <- (s1 / sum_m1) - (B %*% (sum_m2 / sum_m1))
        a <- rep(0,q) #H0
        B <- s2  %*% solve(sum_m3, tol = 1e-30) # este estimador cambia por la hipotesis H0
        phi <- (1 / (n * r)) * s3
        mu_x <- sum_m2 / c(sum_m1)
        Sigma_x <- (1 / n) * s4
        
        theta <-
          c(a,
            B,
            phi,
            mu_x,
            matrixcalc::vech(Sigma_x),
            lambda,
            eta)
        
        #criterio
        lk <- lLz(
          X,
          Y,
          theta,
          type = type
        )$logL
        loglik[count] <- lk
        
        if (StopCrit == "logL") {
          lk1 <- lk
          criterio <- abs(lk1 / lk0 - 1)
          lk0 <- lk1
          if (count == MaxIter) {
            criterio <- 1e-11
          }
          print(paste("criterio de convergencia simple:", criterio))
          print(paste("loglikelihood", lk))
          print(c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x), eta))
        }
        if (StopCrit == "Aitken") {
          if (count < 2) {
            lk1 <- lk
            criterio <- abs(lk1 / lk0 - 1)
          }
          if (count == 2) {
            lk2 <- lk
            c1 <- (lk2 - lk1) / (lk1 - lk0)
            l_A2 <- lk1 + ((lk2 - lk1) / (1 - c1))
            criterio <- abs(lk2 / lk1 - 1)
          }
          if (count > 2) {
            lk3 <- lk
            c <- (lk3 - lk2) / (lk2 - lk1)
            l_A3 <- lk2 + ((lk3 - lk2) / (1 - c))
            criterio <- abs(l_A3 - l_A2)
            l_A2 <- l_A3
            lk1 <- lk2
            lk2 <- lk3
          }
          if (count == MaxIter) {
            criterio <- 1e-11
          }
          print(paste("criterio de convergencia Aitken:", criterio))
          print(paste("loglikelihood", lk))
          print(c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x), eta))
          
        }
        crits[count] <- criterio
        
        # theta <-
        #   c(a,
        #     B,
        #     phi,
        #     mu_x,
        #     matrixcalc::vech(Sigma_x),
        #     lambda,
        #     eta)
      }
    }
    if (type == "S-GH") {
      # S-GH ----
      lambda <- lambda
      gamma_x <- rep(0, p)
      eta <- eta
      theta <-c(a,
                 B,
                 phi,
                 mu_x,
                 matrixcalc::vech(Sigma_x),
                 lambda,
                 eta)
      lk0 <-
        lLz(
          X,
          Y,
          theta,
          type = type,
        )$logL
      loglik <- NULL
      criterio <- 1
      crits <- NULL
      count <- 0
      while (criterio > precision) {
        sum_m1 <- 0
        sum_m2 <- matrix(0, nrow = p, ncol = 1)
        sum_m3 <- matrix(0, nrow = p, ncol = p)
        mi2 <- matrix(0, nrow = p, ncol = 1)
        mi3 <- matrix(0, nrow = p, ncol = p)
        s1 <- rep(0, q) #\sum_{i=1}^{n}m_{i1}\bY_i
        s2 <- matrix(0, nrow = q, ncol = p)   #\sum_{i=1}^{n}\bY_i \bm_{i2}^{\top} 
        s3 <- 0 # \sum_{i=1}^{n} \Big\{ m_{i1}\Big[\bX_i^{\top} \bX_i + 
        # (\bY_i-\hat{\ba})^{\top}(\bY_i-\hat{\ba})\Big] - \bX_i^{\top}\bm_{i2} 
        #- (\bY_i - \hat{\ba})^{\top}\hat{\bB} \bm_{i2} - \bm_{i2}^{\top}\bX_i 
        #- \bm_{i2}^{\top}\hat{\bB}^{\top}(\bY-\hat{\ba}) + \text{tr}(\bm_{i3}) 
        # + \text{tr}(\bm_{i3}\hat{\bB}^{\top}\hat{\bB})\Big\}
        s4 <- matrix(0, nrow = p, ncol = p) # \sum_{i=1}^{n}\Big( \bm_{i3} - 
        # \bm_{i2}\hat{\bmu}_x^{\top} - 
        # \hat{\bmu}_x\bm_{i2}^{\top} + 
        # m_{i1}\hat{\bmu}_x\hat{\bmu}_x^{\top}\Big)
        count <- count + 1
        print(count)
        alpha <- c(rep(0, p), a)
        Lambda <- rbind(diag(1, p), B)
        mu_z <- c(alpha + Lambda %*% mu_x)
        Sigma_z <- Lambda %*% Sigma_x %*% t(Lambda) + diag(phi, r, r)
        Sigma_z <- 0.5 * (Sigma_z + t(Sigma_z))
        Sigma_z_inv <- solve(Sigma_z, tol = 1e-30)
        Delta <- solve(Sigma_x, tol = 1e-30) + ((1 / phi) * crossprod(Lambda))
        Delta_inv <- solve(Delta, tol = 1e-30)
        A <- diag(p) - (1 / phi) * Delta_inv %*% crossprod(Lambda)
        lambda_vz <- lambda - r / 2
        for (i in 1:n) {
          delta_i <-  t(Z[i, ] - mu_z) %*% Sigma_z_inv %*% (Z[i, ] - mu_z)
          m1_i <-
            ghyp::Egig(
              lambda = lambda_vz,
              chi = delta_i + eta,
              psi = eta,
              func = c("1/x")
            )
          b_i <-  (1 / phi) * Delta_inv %*% t(Lambda) %*% (Z[i, ] - mu_z)
          m2_i <- m1_i * (mu_x + b_i)
          m3_i <- Delta_inv + m1_i * tcrossprod(mu_x + b_i)
          m3_i <- 0.5 * (m3_i + t(m3_i)) # porque es simetrica la matriz
          sum_m1 <- sum_m1 + m1_i
          sum_m2 <- sum_m2 + m2_i
          sum_m3 <- sum_m3 + m3_i
          s1 <- s1 + (m1_i * Y[i, ])  
          s2 <- s2 + (Y[i, ] %*% t(m2_i))
          s3 <- s3 + c(
            m1_i * (crossprod(X[i, ]) + crossprod(Y[i, ] - a)) -
              crossprod(X[i, ], m2_i) -
              t(Y[i, ] - a) %*% B %*% m2_i -
              crossprod(m2_i, X[i, ]) -
              t(m2_i) %*% t(B) %*% (Y[i, ] - a) +
              sum(diag(m3_i)) +
              sum(diag(m3_i %*% t(B) %*% B))
          )
          s4 <- s4 + m3_i -
            tcrossprod(m2_i, mu_x) -
            tcrossprod(mu_x, m2_i) +
            m1_i * tcrossprod(mu_x)
        }
        # parameters update
        
        #a <- (s1 / sum_m1) - (B %*% (sum_m2 / sum_m1))
        a <- rep(0,q) #H0
        B <- s2  %*% solve(sum_m3, tol = 1e-30) # este estimador cambia por la hipotesis H0
        phi <- (1 / (n * r)) * s3
        mu_x <- sum_m2 / c(sum_m1)
        Sigma_x <- (1 / n) * s4
        
        theta <-
          c(a,
            B,
            phi,
            mu_x,
            matrixcalc::vech(Sigma_x),
            lambda,
            eta)
        
        #criterio
        lk <- lLz(
          X,
          Y,
          theta,
          type = type
        )$logL
        loglik[count] <- lk
        
        if (StopCrit == "logL") {
          lk1 <- lk
          criterio <- abs(lk1 / lk0 - 1)
          lk0 <- lk1
          if (count == MaxIter) {
            criterio <- 1e-11
          }
          print(paste("criterio de convergencia simple:", criterio))
          print(paste("loglikelihood", lk))
          print(c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x), eta))
        }
        if (StopCrit == "Aitken") {
          if (count < 2) {
            lk1 <- lk
            criterio <- abs(lk1 / lk0 - 1)
          }
          if (count == 2) {
            lk2 <- lk
            c1 <- (lk2 - lk1) / (lk1 - lk0)
            l_A2 <- lk1 + ((lk2 - lk1) / (1 - c1))
            criterio <- abs(lk2 / lk1 - 1)
          }
          if (count > 2) {
            lk3 <- lk
            c <- (lk3 - lk2) / (lk2 - lk1)
            l_A3 <- lk2 + ((lk3 - lk2) / (1 - c))
            criterio <- abs(l_A3 - l_A2)
            l_A2 <- l_A3
            lk1 <- lk2
            lk2 <- lk3
          }
          if (count == MaxIter) {
            criterio <- 1e-11
          }
          print(paste("criterio de convergencia Aitken:", criterio))
          print(paste("loglikelihood", lk))
          print(c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x), eta))
          
        }
        crits[count] <- criterio
        # 
        # theta <-
        #   c(a,
        #     B,
        #     phi,
        #     mu_x,
        #     matrixcalc::vech(Sigma_x),
        #     lambda,
        #     eta)
      }
    }
    
    names(theta) <- c(nam(p,q),"lambda", "eta")
    resultados <- list(
      a = a,
      B = B,
      phi = phi,
      mu.x = mu_x,
      Sigma.x = Sigma_x,
      lambda = lambda,
      eta = eta,
      theta = theta,
      logLz = lk,
      logs = loglik, # despues se puede quitar
      iter = count,
      crits = crits, # despues se puede quitar
      aic = 2 * d - 2 * lk
      
    )
    resultados
    
  }

# en esta función se ajusta el modelo para una grilla de etas usando
# fit_SGH_MRMME() y proporciona eta estimado por verosimilitud profile.

fit_SGH_MRMME_profile_R <- function(X,Y,etas=seq(0.01,5,length = 20), type = "S-NIG"){
  # ajuste del modelo para cada eta
  logliks <- data.frame()
  for (eta in etas){
    fit <- fit_SGH_MRMME_R(X = X, Y = Y,
                           eta = eta, type = type)
    temp_df <- data.frame(eta = eta, logL = fit$logLz)
    logliks <- rbind(logliks, temp_df)
  }
  
  # grafico de la verosimiliu profile
  plot_prof <- ggplot() + 
    geom_line(data = logliks ,aes(x = eta, y = logL)) +
    theme_bw()
  
  # eta optimo
  hat_eta <- logliks %>% arrange(desc(logL)) %>% slice(1)
  
  # ajuste con el eta optimo
  fit_gh <- fit_SGH_MRMME_R(X = X, Y = Y,
                          eta = hat_eta$eta, 
                          type = type)
  lLz_i <- lLz(X, Y, 
                 fit_gh$theta, 
                 type = type)$logL_i
  
  
  resultados <- list(
    a = fit_gh$a,
    B = fit_gh$B,
    phi = fit_gh$phi,
    mu.x = fit_gh$mu.x,
    Sigma.x = fit_gh$Sigma.x,
    lambda = fit_gh$lambda,
    eta = fit_gh$eta,
    theta = fit_gh$theta,
    logLz = fit_gh$logLz,
    logLz_i = lLz_i,
    iter = fit_gh$iter,
    plot_profile = plot_prof,
    logliks = logliks,
    aic = fit_gh$aic
  )
  resultados
}

# modelo restringido normal MRMME con algoritmo EM (bajo H0: a=0)

fit_EM_N_R <- function(X, Y, crit = 1e-10) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2
  Z <- cbind(X, Y)
  
  # Initial Values
  mod_lm <- stats::lm(Y ~ X)
  a <- rep(0,q) # h0
  B <- matrix(stats::coef(mod_lm)[-1,], ncol=p) # lm
  phi <- 1
  mu_x <- colMeans(X) # observed
  Sigma_x <- stats::cov(X) # observed
  
  # a <- rep(0, p)
  # B2 <- matrix(0, nrow = q, ncol = p2)
  # phi <- 1
  # mu_x <- rep(0, p)
  # Sigma_x <- diag(1, p)
  
  #log-likelihood
  theta_est <- c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x))
  logL_k <- logL(theta_est, X, Y)
  
  # Iterations
  s <- 0
  dif <- 1
  
  while (dif > crit) {
    s <- s + 1
    
    # E - step
    Lambda <- rbind(diag(1, p, p), B)
    alpha <- c(rep(0, p), a)
    
    # psi:  Z_i variance
    psi <- Lambda %*% Sigma_x %*% t(Lambda) + diag(phi, r, r)
    psi_inv <- Rfast::spdinv(psi)
    
    
    # m y M (x_i|Z_i,theta) ~ N(m_i,M)
    aux <- t(Z) - c(alpha + Lambda %*% mu_x)
    m <- t(mu_x + tcrossprod(Sigma_x, Lambda) %*% psi_inv %*% aux)
    M <-  Sigma_x - tcrossprod(Sigma_x, Lambda) %*% psi_inv %*% Lambda %*% Sigma_x
   
    # M - step
    y_bar <- colMeans(Y)
    m_bar <- colMeans(m)
    
    # S_ym
    S_ym <- 1 / n * tcrossprod(t(Y) - y_bar, t(m) - m_bar)
    # S_mm
    S_mm <- 1 / n * tcrossprod(t(m) - m_bar)
    
    # S_ym_r
    S_ym_r <- (1 / n) * crossprod(Y, m)
    # S_mm
    S_mm_r <- (1 / n) * crossprod(m)

    # theta_hat:
    # Sigma
    Sigma_x <- S_mm + M
    Sigma_x <- 0.5 * (Sigma_x + t(Sigma_x))
    # mu
    mu_x <- m_bar
    # B
    B <- S_ym_r %*% Rfast::spdinv(S_mm_r + M) # este estimador cambia en el modelo reducido
    # a
    a <- matrix(rep(0,q), ncol=p, nrow=q)
    row.names(a) <- colnames(Y)
    colnames(a) <- "intercept"
    # phi_hat
    #  S(theta_hat)
    #  use ||A||^2 = tr(A^t A)
    S_theta_v <- NULL
    for (i in 1:n) {
      S_theta_v[i] <- sum(diag(crossprod(X[i,] - m[i,]))) +
        sum(diag(crossprod(Y[i,] - a - B %*% m[i,]))) +
        sum(diag(M + t(B) %*% B %*% M))
    }
    S_theta <- sum(S_theta_v)
    phi <- c("phi" = S_theta / (n * r))
    
    #log-likelihood
    theta_est <- c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x))
    logL_k1 <- logL(theta_est, X, Y)
    dif <- abs((logL_k1 / logL_k) - 1)
    logL_k <- logL_k1
  }
  names(theta_est) <- nam(p, q)
  AIC <- 2 * d - 2 * logL(theta_est, X = X, Y = Y)
  BIC <- d * log(n) - 2 * logL(theta_est, X = X, Y = Y)
  
  res <- list(
    "a" = a,
    "B" = B,
    "phi" = phi,
    "mu_x" = mu_x,
    "Sigma_x" = Sigma_x,
    "coef" = cbind(a, B),
    "theta" = theta_est,
    "d" = d,
    "iter" = s,
    "AIC" = AIC,
    "BIC" = BIC,
    "logL" = logL_k
  )
  res
}

# Influence diagnostics ---------------------------------------------------


influence_SGH <- function(theta, X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Z <- cbind(X, Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2
  p_s <- p * (p + 1) / 2
  
  # theta --> (a,B,phi,mu,Sigma)
  theta <- c(theta)
  a <- matrix(theta[1:q], ncol = 1, nrow = q)
  B <- matrix(theta[(q + 1):(q + p * q)], nrow = q, ncol = p)
  phi <- theta[q + p * q + 1]
  mu_x <- theta[(q + p * q + 2):(q + p * q + 1 + p)]
  Sigma_x <- ks::invvech(theta[(q + p * q + 2 + p):d])
  Sigma_x <- 0.5 * (Sigma_x + t(Sigma_x))
  lambda <- theta[(d + 1)]
  eta <- theta[(d + 2)]
  
  alpha <- c(rep(0, p), a)
  Lambda <- rbind(diag(1, p), B)
  
  mu_z <- c(alpha + Lambda %*% mu_x)
  Sigma_z <- Lambda %*% Sigma_x %*% t(Lambda) + diag(phi, r, r)
  Sigma_z <- 0.5 * (Sigma_z + t(Sigma_z))
  Sigma_z_inv <- solve(Sigma_z, tol = 1e-30)
  
  Sigma_x_inv <- Rfast::spdinv(Sigma_x)
  
  Delta <- solve(Sigma_x, tol = 1e-30) + ((1 / phi) * crossprod(Lambda))
  Delta_inv <- solve(Delta, tol = 1e-30)
  lambda_vz <- lambda - r / 2
  
  if (p == 1) {
    D_p <- matrix(1)
  } else {
    D_p <- matrixcalc::duplication.matrix(p)
  }
  
  # Qpp ---------------------------------------------------------------------
  # second derivative of Q
  
  Q_a <- matrix(0, nrow = q, ncol = q)
  Q_B <- matrix(0, nrow = p * q, ncol = p * q)
  Q_phi <- matrix(0, nrow = 1, ncol = 1)
  Q_mu_x <- matrix(0, nrow = p, ncol = p)
  Q_Sigma_x <-  matrix(0, nrow = p_s, ncol = p_s)
  Q_aB <- matrix(0, nrow = q, ncol = p * q)
  Q_a_phi <- matrix(0, nrow = q, ncol = 1)
  Q_Bphi <- matrix(0, nrow = p * q, ncol = 1)
  Q_mu_x_Sigma_x <- matrix(0, nrow = p, ncol = p_s)
  
  for (i in 1:n) {
    delta_i <-  t(Z[i, ] - mu_z) %*% Sigma_z_inv %*% (Z[i, ] - mu_z)
    m1_i <- ghyp::Egig(
      lambda = lambda_vz,
      chi = delta_i + eta,
      psi = eta,
      func = c("1/x")
    )
    b_i <-  (1 / phi) * Delta_inv %*% t(Lambda) %*% (Z[i, ] - mu_z)
    m2_i <- m1_i * (mu_x + b_i)
    m3_i <- Delta_inv + m1_i * tcrossprod(mu_x + b_i)
    m3_i <- 0.5 * (m3_i + t(m3_i))
    
    # Q_a ----
    Q_a_i <- -(m1_i / phi) * diag(q)
    Q_a <- Q_a + Q_a_i
    
    # Q_B ----
    Q_B_i <-  -(1 / phi) * kronecker(t(m2_i), diag(q))
    Q_B <- Q_B + Q_B_i
    
    # Q_phi ----
    S_i <- c(
      m1_i * (crossprod(X[i, ]) + crossprod(Y[i, ] - a)) -
        crossprod(X[i, ], m2_i) -
        t(Y[i, ] - a) %*% B %*% m2_i -
        crossprod(m2_i, X[i, ]) -
        t(m2_i) %*% t(B) %*% (Y[i, ] - a) +
        sum(diag(m3_i)) +
        sum(diag(m3_i %*% t(B) %*% B))
    )
    Q_phi_i <- -1 / (2 * phi ^ 2) * ((2 / phi) * S_i - r)
    Q_phi <- Q_phi + Q_phi_i
    
    # Q_mu_x ----
    Q_mu_x_i <- -m1_i %*% Sigma_x_inv
    Q_mu_x <- Q_mu_x + Q_mu_x_i
    
    # Q_Sigma_x ----
    A_i <- m3_i - tcrossprod(m2_i, mu_x) - tcrossprod(mu_x, m2_i) +  m1_i * tcrossprod(mu_x)
    aux_Sig1 <- 1 / 2 * Sigma_x_inv - Sigma_x_inv %*% A_i %*% Sigma_x_inv
    Q_Sigma_x_i <- t(D_p) %*% kronecker(Sigma_x_inv, aux_Sig1) %*% D_p
    Q_Sigma_x <- Q_Sigma_x + Q_Sigma_x_i
    
    # Q_aB ----
    Q_aB_i <- -1 / phi * kronecker(m2_i, diag(q))
    Q_aB <- Q_aB + Q_aB_i
    
    # Q_a_phi ----
    Q_a_phi_i <- -1 / (phi ^ 2) * (m1_i * (Y[i, ] - a) - B %*% m2_i)
    Q_a_phi <- Q_a_phi + Q_a_phi_i
    
    # Q_Bphi----
    Q_Bphi_i <-  -1 / (phi ^ 2) * (kronecker(m2_i, (Y[i, ] - a)) - matrix(B %*%
                                                                            m3_i, ncol = 1))
    Q_Bphi <- Q_Bphi + Q_Bphi_i
    
    # Q_mu_x_Sigma_x -----
    aux1 <- t(m2_i - m1_i * mu_x) %*% Sigma_x_inv
    Q_mu_x_Sigma_x_i <-  -kronecker(aux1, Sigma_x_inv) %*% D_p
    Q_mu_x_Sigma_x <- Q_mu_x_Sigma_x + Q_mu_x_Sigma_x_i
  }
  
  Qpp <- pracma::blkdiag(Q_a, Q_B, Q_phi, Q_mu_x, Q_Sigma_x)
  Qpp[1:q, (q + 1):(q * (1 + p))] <- Q_aB
  Qpp[(q + 1):(q * (1 + p)), 1:q] <-  t(Q_aB)
  
  Qpp[1:q, (q * (1 + p) + 1)] <- Q_a_phi
  Qpp[(q * (1 + p) + 1) , 1:q] <-  t(Q_a_phi)
  
  Qpp[(q + 1):(q * (p + 1)), (q * (1 + p) + 1)] <- Q_Bphi
  Qpp[(q * (1 + p) + 1), (q + 1):(q * (p + 1))] <-  t(Q_Bphi)
  
  Qpp[(q * (1 + p) + 2):(q * (1 + p) + 1 + p) , (q * (1 + p) + 1 + p + 1):d] <- Q_mu_x_Sigma_x
  Qpp[(q * (1 + p) + 1 + p + 1):d, (q * (1 + p) + 2):(q * (1 + p) + 1 + p)] <-  t(Q_mu_x_Sigma_x)
  
  
  # Esquema 1 ----
  
  Delta_1 <- matrix(0, nrow = d, ncol = r)
  for (j in 1:r) {
    # e_j, c_qj
    e_j <- matrix(diag(r)[j, ], ncol = 1)
    c_qj <- cbind(matrix(0, nrow = q, ncol = p), diag(1, q, q)) %*% e_j
    c_qj
    Delta_j <- matrix(0, nrow = d, ncol = 1)
    
    for (i in 1:n) {
      delta_i <-  t(Z[i, ] - mu_z) %*% Sigma_z_inv %*% (Z[i, ] - mu_z)
      m1_i <- ghyp::Egig(
        lambda = lambda_vz,
        chi = delta_i + eta,
        psi = eta,
        func = c("1/x")
      )
      b_i <-  (1 / phi) * Delta_inv %*% t(Lambda) %*% (Z[i, ] - mu_z)
      m2_i <- m1_i * (mu_x + b_i)
      m3_i <- Delta_inv + m1_i * tcrossprod(mu_x + b_i)
      m3_i <- 0.5 * (m3_i + t(m3_i))
      
      A_theta_i <- m1_i * tcrossprod(Z[i, ] - alpha) - Lambda %*% m2_i %*% t(Z[i, ] -
                                                                               alpha) -
        (Z[i, ] - alpha) %*% t(m2_i) %*% t(Lambda)  + Lambda %*% m3_i %*% t(Lambda)
      
      # Q_ia_omega ----
      aux_a <- c(t(m1_i * Z[i, ] - alpha - Lambda %*% m2_i) %*% e_j)
      Q_ia_omega  <- 1 / phi *  aux_a * c_qj
      
      # Q_iB_omega ----
      auxB <- c(t(Z[i, ] - alpha - Lambda %*% m2_i) %*% e_j)
      auxB2 <- m3_i %*% t(Lambda) %*% e_j
      Q_iB_omega  <- 1 / phi * (auxB * kronecker(m2_i, c_qj)  - kronecker(auxB2, c_qj))
      
      # # Q_iphi_omega -----
      Q_iphi_omega  <- 1 / (2 * phi ^ 2) * sum(diag(t(e_j) %*% A_theta_i %*% e_j))
      Q_imu_omega  <- matrix(0, nrow = p, ncol = 1)
      Q_iSigma_omega  <- matrix(0, nrow = p_s, ncol = 1)
      
      aux_delta_j <- rbind(Q_ia_omega,
                           Q_iB_omega,
                           Q_iphi_omega,
                           Q_imu_omega,
                           Q_iSigma_omega)
      Delta_j <- Delta_j + aux_delta_j
    }
    Delta_1[, j] <- Delta_j
  }
  rownames(Delta_1) <- nam(p, q)
  
  
  # Esquema 2 ---------------------------------------------------------------
  
  Delta_2 <- matrix(0, nrow = d, ncol = n)
  for (i in 1:n) {
    b_i <-  (1 / phi) * Delta_inv %*% t(Lambda) %*% (Z[i, ] - mu_z)
    m2_i <- m1_i * (mu_x + b_i)
    m3_i <- Delta_inv + m1_i * tcrossprod(mu_x + b_i)
    m3_i <- 0.5 * (m3_i + t(m3_i))
    
    A_theta_i <- m1_i * tcrossprod(Z[i, ] - alpha) - Lambda %*% m2_i %*% t(Z[i, ] -
                                                                             alpha) -
      (Z[i, ] - alpha) %*% t(m2_i) %*% t(Lambda)  + Lambda %*% m3_i %*% t(Lambda)
    
    # Q_ia_omega ----
    aux_a2 <- Y[i, ] - a - B %*% m2_i
    Q_ia_omega  <- 1 / phi * aux_a2
    
    # Q_iB_omega ----
    aux_B_2 <- kronecker(m2_i, (Y[i, ] - a))
    Q_iB_omega  <- 1 / phi * (aux_B_2 - c(B %*% m3_i))
    
    # Q_iphi_omega ----
    Q_iphi_omega  <- (1 / (2 * phi ^ 2)) * sum(diag(A_theta_i))
    
    Q_imu_omega  <- matrix(0, nrow = p, ncol = 1)
    Q_iSigma_omega  <- matrix(0, nrow = p_s, ncol = 1)
    
    Delta_i <-
      c(Q_ia_omega,
        Q_iB_omega,
        Q_iphi_omega,
        Q_imu_omega,
        Q_iSigma_omega)
    Delta_2[, i] <- Delta_i
  }
  rownames(Delta_2) <- nam(p, q)
  
  # B_inf1 -------------------------------------------------------------------
  F_inf_1 <- t(Delta_1) %*% Rfast::spdinv(-Qpp) %*% Delta_1
  F_inf_1 <-  2 * abs(diag(F_inf_1)) / matrixcalc::matrix.trace(2 * F_inf_1)
  cut_point_1 <- mean(F_inf_1) + 2 * stats::sd(F_inf_1)
  names(F_inf_1) <-  colnames(Z)
  
  # B_inf2 -------------------------------------------------------------------
  F_inf_2 <- t(Delta_2) %*% Rfast::spdinv(-Qpp) %*% Delta_2
  F_inf_2 <-  2 * abs(diag(F_inf_2)) / matrixcalc::matrix.trace(2 * F_inf_2)
  cut_point_2 <- mean(F_inf_2) + 2 * stats::sd(F_inf_2)
  
  inf_1 <- which(round(F_inf_1, 2) > cut_point_1)
  inf_2 <- which(round(F_inf_2, 2) > cut_point_2)
  
  # Plot scheme 1----
  # p1 <-
  #   ggplot(data.frame(F_inf_1, ind = names(F_inf_1)), aes(label = ind)) +
  #   #scale_y_continuous(limits = c(0,1)) +
  #   geom_hline(aes(yintercept = cut_point_1), col = "red") +
  #   geom_point(aes(x = ind, y = F_inf_1)) +
  #   labs(x = "", y = "F values in Scheme I") +
  #   geom_text(aes(
  #     x = ind,
  #     y = F_inf_1,
  #     label = ifelse(F_inf_1 >= cut_point_1, ind, "")
  #   ),
  #   hjust = 0,
  #   vjust = -1) +
  #   theme_bw()
  # 
  
  # Plot scheme 2----
  
  # #p2 <-
  #   ##ggplot(data.frame(F_inf_2, ind = 1:length(F_inf_2)), aes(label = ind)) +
  #   ggplot(data.frame(F_inf_2, ind = 1:length(F_inf_2)), aes(label = ind)) +
  #   #scale_y_continuous(limits = c(0,1)) +
  #   geom_hline(aes(yintercept = cut_point_2), col = "red") +
  #   geom_point(aes(x = ind, y = F_inf_2)) +
  #   labs(x = "", y = "F values in Scheme II") +
  #   geom_text(
  #     aes(
  #       x = ind,
  #       y = F_inf_2,
  #       label = ifelse(round(F_inf_2, 2) > cut_point_2, ind, "")
  #     ),
  #     hjust = -1,
  #     vjust = 1,
  #     size = 2
  #   ) +
  #   theme_bw() +
  #   theme(panel.grid.minor.x = element_blank())
  # 
  # library(patchwork)
  # p2t <- p1 | p2
  # p2t <- p2t  +
  #   plot_annotation(
  #     tag_levels = "a",
  #     tag_prefix = '(',
  #     tag_sep = '',
  #     tag_suffix = ')'
  #   ) &
  #   theme(plot.tag.position = "top",
  #         plot.tag = element_text(
  #           face = "bold",
  #           size = 8,
  #           vjust = 5
  #         ))
  # 
  res = list(
    "G1" = F_inf_1,
    "G2" = F_inf_2,
    "cut_point_G1" = cut_point_1,
    "cut_point_G2" = cut_point_2,
    "inf_scheme_1" =  inf_1,
    "inf_scheme_2" =  inf_2#,
    #"plot1" = p1,
    #"plot2" = p2,
    #"plot" = p2t
  )
  
  return(res)
  
}


# *modelo normal----
influence_N <- function(theta, X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Z <- cbind(X, Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2
  p_s <- p * (p + 1) / 2
  
  # theta --> (a,B,phi,mu,Sigma)
  theta <- as.vector(theta)
  a <- matrix(theta[1:q], ncol = 1, nrow = q)
  B <- matrix(theta[(q + 1):(q + p * q)], nrow = q, ncol = p)
  phi <- theta[q + p * q + 1]
  mu_x <- theta[(q + p * q + 2):(q + p * q + 1 + p)]
  Sigma_x <- ks::invvech(theta[(q + p * q + 2 + p):d])
  Sigma_x <- 0.5 * (Sigma_x + t(Sigma_x))
  
  # eta and psi
  Lambda <- rbind(diag(1, p, p), B)
  alpha <- c(rep(0, p), a)
  eta <- c(alpha + Lambda %*% mu_x)
  psi <- Lambda %*% tcrossprod(Sigma_x, Lambda) + diag(phi, r, r)
  psi <- 0.5 * (psi + t(psi))
  psi_inv <- Rfast::spdinv(psi)
  
  # m_i y M
  aux_m <- t(Z) - c(alpha + Lambda %*% mu_x)
  m <- t(mu_x + tcrossprod(Sigma_x, Lambda) %*% psi_inv %*% aux_m)
  M <-
    Sigma_x - tcrossprod(Sigma_x, Lambda) %*% psi_inv %*% Lambda %*% Sigma_x
  m_bar <- colMeans(m)
  I_q <- diag(1, q, q)
  
  if (p == 1) {
    D_p <- matrix(1)
  } else {
    D_p <- matrixcalc::duplication.matrix(p)
  }
  
  S_mm <- 1 / n * tcrossprod(t(m) - m_bar)
  Sigma_x_inv <- Rfast::spdinv(Sigma_x)
  
  # Qpp ---------------------------------------------------------------------
  # second derivative of Q
  Q_a <- -(n / phi) * diag(1, q, q)
  Q_aB <- -(n / phi) * kronecker(t(m_bar), I_q)
  aux_B1 <-  kronecker(S_mm + tcrossprod(m_bar), I_q)
  aux_B2 <- kronecker(M, I_q)
  Q_B <- -(n / phi) * (aux_B1 + aux_B2)
  Q_phi <- matrix(-(n * r) / (2 * phi ^ 2), nrow = 1, ncol = 1)
  Q_mu <- -n * Sigma_x_inv
  aux_Sig1 <-
    kronecker(Sigma_x_inv %*% M %*% Sigma_x_inv, Sigma_x_inv)
  aux_Sig2 <- kronecker(Sigma_x_inv, Sigma_x_inv)
  Q_Sigma <- n * t(D_p) %*% (aux_Sig1 - (1 / 2 * aux_Sig2)) %*% D_p
  Qpp <- pracma::blkdiag(Q_a, Q_B, Q_phi, Q_mu, Q_Sigma)
  Qpp[1:q, (q + 1):(q * (1 + p))] <- Q_aB
  Qpp[(q + 1):(q * (1 + p)),1:q] <-  t(Q_aB) 
  
  # Esquema 1 ---------------------------------------------------------------
  Delta_1 <- matrix(0, nrow = d, ncol = r)
  for (j in 1:r) {
    # e_j, c_qj
    e_j <- diag(r)[j,]
    c_qj <-
      cbind(matrix(0, nrow = q, ncol = p), diag(1, q, q)) %*% e_j
    Delta_j <- matrix(0, nrow = d, ncol = 1)
    for (i in 1:n) {
      aux_delta <- c(t(Z[i,] - alpha - Lambda %*% m[i,]) %*% e_j)
      kron_prod <- kronecker(m[i,], c_qj)
      Q_ia_omega  <- 1 / phi * aux_delta * c_qj
      Q_iB_omega  <- 1 / phi * aux_delta * kron_prod
      Q_iphi_omega  <- 1 / (2 * phi ^ 2) * aux_delta ^ 2
      Q_imu_omega  <- matrix(0, nrow = p, ncol = 1)
      Q_iSigma_omega  <- matrix(0, nrow = p_s, ncol = 1)
      aux_delta_j <-
        rbind(Q_ia_omega,
              Q_iB_omega,
              Q_iphi_omega,
              Q_imu_omega,
              Q_iSigma_omega)
      Delta_j <- Delta_j + aux_delta_j
    }
    Delta_1[, j] <- Delta_j
  }
  rownames(Delta_1) <- nam(p, q)
  
  
  # Esquema 2 ---------------------------------------------------------------
  
  Delta_2 <- matrix(0, nrow = d, ncol = n)
  for (i in 1:n) {
    aux_Y <- t(Y[i,] - a - B %*% m[i,])
    Q_ia_omega  <- 1 / phi * aux_Y
    Q_iB_omega  <- 1 / phi * (kronecker(m[i, ], aux_Y) - c(B %*% M))
    # use ||A||^2 = tr(A^t A)
    aux_Z <- Z[i,] - alpha - Lambda %*% m[i,]
    S_i <-
      sum(diag(crossprod(aux_Z))) + sum(diag(M + t(B) %*% B %*% M))
    Q_iphi_omega  <- 1 / (2 * phi ^ 2) * S_i
    Q_imu_omega  <- matrix(0, nrow = p, ncol = 1)
    Q_iSigma_omega  <- matrix(0, nrow = p_s, ncol = 1)
    Delta_i <-
      c(Q_ia_omega,
        Q_iB_omega,
        Q_iphi_omega,
        Q_imu_omega,
        Q_iSigma_omega)
    Delta_2[, i] <- Delta_i
  }
  rownames(Delta_2) <- nam(p, q)
  
  # B_inf1 -------------------------------------------------------------------
  F_inf_1 <- t(Delta_1) %*% Rfast::spdinv(-Qpp) %*% Delta_1
  B_inf_1 <-
    2 * abs(diag(F_inf_1)) / matrixcalc::matrix.trace(2 * F_inf_1)
  cut_point_1 <- mean(B_inf_1) + 2 * stats::sd(B_inf_1)
  names(B_inf_1) <-  colnames(Z)
  # B_inf2 -------------------------------------------------------------------
  F_inf_2 <- t(Delta_2) %*% Rfast::spdinv(-Qpp) %*% Delta_2
  B_inf_2 <-     2 * abs(diag(F_inf_2)) / matrixcalc::matrix.trace(2 * F_inf_2)
  cut_point_2 <- mean(B_inf_2) + 2 * stats::sd(B_inf_2)
  
  inf_1 <- which(round(B_inf_1,2) > cut_point_1)
  inf_2 <- which(round(B_inf_2,2) > cut_point_2)
  
  # Plot scheme 1----
  p1 <-
    ggplot(data.frame(B_inf_1, ind = 1:length(B_inf_1)), aes(label = ind)) +
    #scale_y_continuous(limits = c(0,1)) +
    geom_hline(aes(yintercept = cut_point_1), col = "red") +
    geom_point(aes(x = ind, y = B_inf_1)) +
    labs(x = "", y = "F values in Scheme I") +
    geom_text(aes(
      x = ind,
      y = B_inf_1,
      label = ifelse(B_inf_1 >= cut_point_1, ind, "")
    ),
    hjust = 0,
    vjust = -1) +
    theme_bw()
  
  
  # Plot scheme 2----
  
  p2 <-
    ggplot(data.frame(B_inf_2, ind = 1:length(B_inf_2)), aes(label = ind)) +
    #scale_y_continuous(limits = c(0,1)) +
    geom_hline(aes(yintercept = cut_point_2), col = "red") +
    geom_point(aes(x = ind, y = B_inf_2)) +
    labs(x = "", y = "F values in Scheme II") +
    geom_text(
      aes(
        x = ind,
        y = B_inf_2,
        label = ifelse(round(B_inf_2,2) > cut_point_2, ind, "")
      ),
      hjust = -1,
      vjust = 1,
      size = 2
    ) +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank())
  
  library(patchwork)
  p2t <- p1 | p2
  p2t <- p2t  +
    plot_annotation(
      tag_levels = "a",
      tag_prefix = '(',
      tag_sep = '',
      tag_suffix = ')'
    ) &
    theme(plot.tag.position = "top",
          plot.tag = element_text(
            face = "bold",
            size = 8,
            vjust = 5
          ))
  
  res = list(
    "G1" = B_inf_1,
    "G2" = B_inf_2,
    "cut_point_G1" = cut_point_1,
    "cut_point_G2" = cut_point_2,
    "inf_scheme_1" =  inf_1,
    "inf_scheme_2" =  inf_2,
    "plot1" = p1,
    "plot2" = p2,
    "plot" = p2t
  )
  
  return(res)
  
}
