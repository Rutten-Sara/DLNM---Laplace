DLNM_Laplace <- function(model,
                         crossbasis,
                         ID,
                         covar.ri = "ind",
                         data,
                         vx,
                         vl){
  # Required packages
  required_packages <- c("Matrix", "spdep")
  for (package in required_packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      message(paste("Package", package, "is required but not installed."))
    }
  }
  
  if(!missing(data)){
    mod <- stats::model.frame(model, data = data)
    y <- as.numeric(stats::model.extract(mod, "response"))[!is.na(crossbasis[,1])] # response
    Z.linear <- Matrix::sparse.model.matrix(mod, data = data)[!is.na(crossbasis[,1]), ,drop = FALSE] # linear covariates
    p <- ncol(Z.linear)
  } else{
    mod <- stats::model.frame(model)
    y <- as.numeric(stats::model.extract(mod, "response"))[!is.na(crossbasis[,1])] # response
    Z.linear <- Matrix::sparse.model.matrix(mod)[!is.na(crossbasis[,1]), ,drop = FALSE] # linear covariates
    p <- ncol(Z.linear)
    }
  
  Wcb <- Matrix::Matrix(na.omit(crossbasis), sparse = TRUE) # cross basis
  zeta <- 1e-05 # precision for linear effect coefficient
  # Hyperparameters for penalty parameters
  a <- b <- 10^(-5)
  nu <- 3
  logpv.rand <- function(v) 0
  
  
  # Difference order of the penalty
  penorder <- 2
  
  # Penalty for exposure
  Dx <- Matrix::Diagonal(vx+1, x = NULL)
  for (k in 1:penorder) Dx <- Matrix::diff(Dx)
  Px <- Matrix::t(Dx) %*% Dx
  Px <- Px[-1,-1]
  Px <- Px + Matrix::Diagonal(vx, x = 1e-12)
  
  # Penalty for lag variable
  Dl <- Matrix::Diagonal(vl, x = NULL)
  for (k in 1:penorder) Dl <- Matrix::diff(Dl)
  Pl <- Matrix::t(Dl) %*% Dl
  Pl <- Pl + Matrix::Diagonal(vl, x = 1e-12)
  
  # Additional penalty for delay (force lag effect going to zero)
  Dl_add <- diag((0:(vl-1))^2)
  #Dl_add <- diag(rep(0:1,c(6,4)))
  Pl_add <- Dl_add + diag(1e-12, vl)
  
  Pv <- function(v) exp(v[1])*(kronecker(Px,Matrix::Diagonal(n = vl, x = NULL))) +
    exp(v[2])*(kronecker(Matrix::Diagonal(n = vx, x = NULL),Pl)) +
    exp(v[3])*(kronecker(Matrix::Diagonal(n = vx, x = NULL),Pl_add))
  
  
  if(!missing(ID)){
    atau <- btau <- 10^(-5)
    a.rho <- b.rho <- 0.5
    if(!missing(data)){
      Z.rand <- Matrix::sparse.model.matrix(~ as.factor(ID) + 0, data)[!is.na(crossbasis[,1]),] # model matrix for random effects
    } else {
      Z.rand <- Matrix::sparse.model.matrix(~ as.factor(ID) + 0)[!is.na(crossbasis[,1]),] # model matrix for random effects
    }
    q.rand <- dim(Z.rand)[2]
    v.rand <- 1

    
    if(covar.ri == "ind"){
      
      Gv <- function(v) Matrix::Diagonal(n = q.rand, x = exp(v[4]))
      logpv.rand <- function(v) 0.5 * (nu + q.rand) * v[4] - 
        (0.5*nu + a)*log(b + 0.5*nu*exp(v[4]))
      
    } else if (covar.ri == "ICAR"){
      neig.map <- spdep::poly2nb(map,row.names = map$MSOA11CD)
      Rn <- matrix(0, nrow = q.rand, ncol =q.rand)
      
      for (s in 1:q.rand) {
        # Diagonal elements (N_s)
        Rn[s, s] <- length(neig.map[[s]])
        
        # Off-diagonal elements (-1 for neighbors)
        for (u in neig.map[[s]]) {
          Rn[s, u] <- -1
        }
      }
      Rn <- Matrix::Matrix(Rn, sparse = TRUE)
      
      Gv <- function(v) exp(v[4])*Rn
      logpv.rand <- function(v) 0.5 * (nu + q.rand) * v[4] - 
        (0.5*nu + a)*log(b + 0.5*nu*exp(v[4]))
      
    } else if (covar.ri == "Convolution") {
      neig.map <- spdep::poly2nb(map,row.names = map$MSOA11CD)
      Rn <- matrix(0, nrow = q.rand, ncol =q.rand)
      
      for (s in 1:q.rand) {
        # Diagonal elements (N_s)
        Rn[s, s] <- length(neig.map[[s]])
        
        # Off-diagonal elements (-1 for neighbors)
        for (u in neig.map[[s]]) {
          Rn[s, u] <- -1
        }
      }
      Rn <- Matrix::Matrix(Rn, sparse = TRUE)
  
      Gv <- function(v) Matrix::bdiag(exp(v[4])*Rn,
            Matrix::Diagonal(n = q.rand, x = exp(v[5])))
      logpv.rand <- function(v) sum(0.5 * (nu + q.rand) * v[4:5]) - 
        sum((0.5*nu + a)*log(b + 0.5*nu*exp(v[4:5])))
      Z.rand <- cbind(Z.rand, Z.rand)
      v.rand <- c(1,1)
      
    } else if (covar.ri == "Leroux"){
      neig.map <- spdep::poly2nb(map,row.names = map$MSOA11CD)
      Rn <- matrix(0, nrow = q.rand, ncol =q.rand)
      
      for (s in 1:q.rand) {
        # Diagonal elements (N_s)
        Rn[s, s] <- length(neig.map[[s]])
        
        # Off-diagonal elements (-1 for neighbors)
        for (u in neig.map[[s]]) {
          Rn[s, u] <- -1
        }
      }
      Rn <- Matrix::Matrix(Rn, sparse = TRUE)
      
      Gv <- function(v) exp(v[4])*(Matrix::Diagonal(n = q.rand, 
                                         x = 1 - exp(v[5])/(1+exp(v[5]))) + 
        Matrix::Matrix(exp(v[5])/(1+exp(v[5]))*Rn, sparse = T))
      Lv <- function(v) (Matrix::Diagonal(n = q.rand, 
                                                    x = 1 - exp(v[5])/(1+exp(v[5]))) + 
                                     Matrix::Matrix(exp(v[5])/(1+exp(v[5]))*Rn, sparse = T))
      logpv.rand <- function(v)  {
        value <- 0.5 * nu * v[4] - (0.5*nu + a)*log(b + 0.5*nu*exp(v[4])) + 
          0.5*sum(sapply(eigen(Lv(v),only.values = T)$values,log)) + 
          a.rho*v[5] - (a.rho + b.rho)*log(1 + exp(v[5]))
        return(as.numeric(value))}
      v.rand <- c(1,1)
      }
    # Global design mtarix
    X <- Matrix::Matrix(cbind(Z.linear,Wcb,Z.rand), sparse = TRUE)
    
    #Precision matrix for parameter xi
    Qv <- function(v) Matrix::bdiag(Matrix::Diagonal(n = p, x = zeta),
                            Pv(v),
                            Gv(v))
  } else {
    # Global design mtarix
    X <- Matrix::Matrix(cbind(Z.linear,Wcb), sparse = TRUE)
    #Precision matrix for parameter xi
    Qv <- function(v) Matrix::bdiag(Matrix::Diagonal(n = p, x = zeta),
                                    Pv(v))
  }
  
  # For Poisson GLM with log-link
  mu <- function(xi) exp(as.numeric(X %*% xi))
  W <- function(xi) Matrix::Diagonal(x = exp(as.numeric(X %*% xi)))
  s <- function(gam) exp(gam)
  
  log_pxi <- function(xi, Qv) {
    post <- sum(y * as.numeric(X %*% xi) -s(X %*% xi)) - .5 * Matrix::t(xi) %*% Qv %*% xi
    return(as.numeric(post))
  }
  
  Grad.logpxi <- function(xi, Qv){
    value <- Matrix::t(X)%*%(y-mu(xi)) - Qv%*%xi
    as.numeric(value)
  }
  
  Hess.logpxi <- function(xi,Qv){
    value <- -Matrix::t(X)%*%W(xi)%*%X - Qv 
    value
  }
  
  
  # Laplace approximation to conditional posterior of xi
  # using Newton-Raphson algorithm
  NR_xi <- function(xi0, Qv){
    
    epsilon <- 1e-03 # Stop criterion
    maxiter <- 100   # Maximum iterations
    iter <- 0        # Iteration counter
    
    for (k in 1:maxiter) {
      dxi <- as.numeric((-1) * solve(Hess.logpxi(xi0, Qv),
                                     Grad.logpxi(xi0, Qv)))
      xi.new <- xi0 + dxi
      step <- 1
      iter.halving <- 1
      logpxi.current <- log_pxi(xi0, Qv)
      while (log_pxi(xi.new, Qv) <= logpxi.current) {
        step <- step * .5
        xi.new <- xi0 + (step * dxi)
        iter.halving <- iter.halving + 1
        if (iter.halving > 30) {
          break
        }
      }
      dist <- sqrt(sum((xi.new - xi0) ^ 2))
      iter <- iter + 1
      xi0 <- xi.new
      if(dist < epsilon) break
    }
    
    xistar <- xi0
    return(xistar)
  }
  # Initial values for log-penalty and log-overdispersion parameter
  
  if(!missing(ID)){
    v_init <- c(rep(5,3), v.rand)
  } else{
    v_init <- c(rep(5,3))
  }
  
  Qv_init <- Qv(v_init)
  xi_init <- NR_xi(xi0 = rep(0,dim(X)[2]), Qv = Qv_init)
  
  # Log-posterior for the penalty- vector
  XWX <- Matrix::t(X)%*%W(xi_init)%*%X
  Xxi <- X%*%xi_init
  log_pv <- function(v){
    Qv <- Qv(v)
    e1 <- eigen(XWX + Qv,only.values = T)$values
    e2 <- eigen(Pv(v),only.values = T)$values
    
    a1 <- 0.5*sum(sapply(e1,log))
    a2 <- sum(y*(Xxi))
    a3 <- sum(s(Xxi))
    a4 <- 0.5 * sum((xi_init * Qv) %*% xi_init)
    a5 <- 0.5*sum(sapply(e2, log))
    a6 <- (0.5*nu + a)*(log(b + 0.5*nu*exp(v[1]))+log(b + 0.5*nu*exp(v[2]))
                        +log(b + 0.5*nu*exp(v[3])))
    a7 <- 0.5*nu*(v[1]+v[2]+v[3])
    
    value <- -a1+a2-a3-a4+a5-a6+a7+logpv.rand(v)
    return(as.numeric(value))
  }
  
  # Mode a posteriori estimate of v
  v_mode <- optim(par = v_init, 
                  fn = log_pv, 
                  method="Nelder-Mead", 
                  control = list(fnscale = -1, reltol = 1e-08))$par
  # Estimate of xi (regression parameter)
  Qv_mode <- Qv(v_mode)
  xi_mode <- NR_xi(xi0 = xi_init, Qv = Qv_mode)
  ind.cb <- seq(dim(Z.linear)[2]+1, by = 1, length.out = dim(crossbasis)[2])
  xi_mode_cb <- xi_mode[ind.cb]
  
  Sigma <- - solve(Hess.logpxi(xi = xi_mode,
                              Qv = Qv_mode))
  
  if(!missing(ID)){
    if(covar.ri == "Convolution"){
      str <- tail(xi_mode, 2*q.rand)[1:q.rand]
      unstr <- tail(xi_mode, q.rand)
      xispat <- (str + unstr) - mean(str + unstr)
    } else {
      xispat <- tail(xi_mode, q.rand) - mean(tail(xi_mode, q.rand))
    } 
  } else {
    xispat <- NULL
  }
  
  output <- list("xi_mode" = xi_mode,
                 "Sigma" = Sigma,
                 "ind.cb" = ind.cb,
                 "xi_mode_cb" = xi_mode_cb,
                 "crossbasis" = crossbasis,
                 "xispat" = xispat,
                 "v_mode" = v_mode)
  
  }
