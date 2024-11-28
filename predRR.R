predRR <- function(model,at_x,cen, L){
  crossbasis <- model$crossbasis
  model_laplace
  ind.cb <- model$ind.cb
  xi_mode_cb <- model$xi_mode_cb
  sigma_xi <- model$Sigma
  predlag <- 0:L # lag 0:L
  nrow <- length(at_x) # number of predictions
  at_x <- rep(at_x, length(predlag)) # Replicate values for every lag
  
  basisvar <- do.call("onebasis", c(list(x = at_x), attr(crossbasis,"argvar")))
  basislag <- do.call("onebasis", c(list(x = predlag), attr(crossbasis,"arglag")))
  
  # basis variables for center
  basiscen <- do.call("onebasis", c(list(x = cen), attr(crossbasis,"argvar")))
  # center data matrix
  Xpred_cen <- scale(basisvar,center = basiscen, scale = F)
  
  # Prediction matrix
  Xpred <- matrix(0, nrow=length(at_x), ncol = ncol(crossbasis))
  
  for (l in seq(length = length(predlag))){
    for (v in seq(length = ncol(basisvar))){
      for (k in 1:ncol(basislag)){
        Xpred[((l-1)*nrow+1):(nrow*l),(ncol(basislag)*(v-1)+k)] = Xpred_cen[((l-1)*nrow+1):(nrow*l),v]*basislag[l,k] # prediction for every variable at lag = l
      } 
    }
  }
  
  ##########################
  logpredX <- as.numeric(Xpred %*% xi_mode_cb)
  sd_logpredX <- as.numeric(sqrt(diag(Xpred %*% sigma_xi[ind.cb, ind.cb] %*% t(Xpred))))
  
  quantiles <- mapply(function(mean_val, sd_val) {
    qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
  }, logpredX, sd_logpredX)
  Qlower_logpredX <- quantiles[1,]
  Qupper_logpredX <- quantiles[2,]
  ###############
  
  # Overall RR
  Xpredall <- 0
  for (j in seq(L+1)) {
    ind_all <- seq(nrow) + nrow * (j - 1) 
    Xpredall <- Xpredall + Xpred[ind_all, , drop = FALSE] # add effect of lag period to cumulative effect
  }
  pred_all <- exp(as.numeric(Xpredall %*% xi_mode_cb))
  sd_all <- sqrt(diag(Xpredall %*% sigma_xi[ind.cb, ind.cb] %*% t(Xpredall)))
  quantiles_all <- mapply(function(mean_val, sd_val) {
    qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
  }, log(pred_all), sd_all)
  Qlower_all = exp(quantiles_all[1,])
  Qupper_all = exp(quantiles_all[2,])

  output <- list("Xpred" = Xpred,
                 "logpredX" = logpredX,
                 "sd_logpredX" = sd_logpredX,
                 "Qlower_logpredX" = Qlower_logpredX,
                 "Qupper_logpredX" = Qupper_logpredX,
                 "Xpredall" = Xpredall,
                 "Qlower_all" = Qlower_all,
                 "Qupper_all" = Qupper_all,
                 "pred_all" = pred_all,
                 "sd_all" = sd_all)
}
