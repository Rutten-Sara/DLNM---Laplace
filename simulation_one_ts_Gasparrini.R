rm(list = ls())

source('DLNM_Laplace.R')
source('predRR.R')

library(dlnm);library(mgcv)

# DEFINE THE EXPOSURE

# REAL TEMPERATURE SERIES, STANDARDIZED IN 0-10
x <- chicagoNMMAPS$temp
x <- (x-min(x))/diff(range(x))*10

# MATRIX Q OF EXPOSURE HISTORIES
Q <- tsModel::Lag(x,0:40)


# BASIC FUNCTIONS TO SIMULATE UNIDIMENSIONAL SHAPES
flin <- function(x) 0.1*x
fflex <- function(x) {
  coef <- c(0.2118881,0.1406585,-0.0982663,0.0153671,-0.0006265)
  as.numeric(outer(x,0:4,'^')%*%coef)
}
fdnorm <- function(x) (dnorm(x,1.5,2)+1.5*dnorm(x,7.5,1))

wconst <- function(lag) lag-lag+0.20
wdecay <- function(lag) exp(-lag/2)
wpeak1 <- function(lag) 12*dnorm(lag,8,5)
wpeak2 <- function(lag) 15*dnorm(lag,8,10)
wdnorm <- function(lag) 5*(dnorm(lag,4,6)+dnorm(lag,25,4))

# FUNCTIONS TO SIMULATE THE BI-DIMENSIONAL EXPOSURE-LAG-RESPONSE
fplane <- function(x,lag) 0.1 * (flin(x)-flin(2)) * wconst(lag)
ftemp <- function(x,lag) 0.1 * (fflex(x)-fflex(5)) * 
  ifelse(is.na(x),NA,ifelse(x>=5,wdecay(lag),wpeak1(lag)))
fcomplex <- function(x,lag) 0.1 * (fdnorm(x)-fdnorm(5)) * 
  ifelse(is.na(x),NA,ifelse(x>=5,wdnorm(lag),wpeak2(lag)))

# COMBINATIONS OF FUNCTIONS USED TO SIMULATE DATA
combsim <- c("fplane","ftemp","fcomplex")
names(combsim) <- c("Plane","Temperature","Complex")

# CENTERING POINT
cen_ind <- c(2,5,5)


# LIST WITH TRUE EFFECT SURFACES OF FOR EACH COMBINATIONS
trueeff <- lapply(combsim, function(fun) {
  temp <- outer(seq(0,10,0.25),0:40,fun)
  dimnames(temp) <- list(seq(0,10,0.25),paste("lag",0:40,sep=""))
  return(temp)
})
names(trueeff) <- names(combsim)

# FUNCTION TO COMPUTE THE CUMULATIVE EFFECT GIVEN AN EXPOSURE HISTORY
fcumeff <- function(hist,lag,fun) sum(do.call(fun,list(hist,lag)))

ind = 1
cumeff <- apply(Q,1,fcumeff,0:40,combsim[ind])

# Plotting the surfaces
library(plot3D)

pdf("scenarios.pdf",height=4,width=12)

layout(matrix(1:3,ncol=3,byrow=TRUE))
par(mar=c(1,1,3,1))

# scenario 1
df = data.frame(x = rep(seq(0,10,0.25)), y = rep(0:40, each = 41))

dens <- akima::interp(x = df$x, 
                      y = df$y, 
                      z = trueeff[[1]], 
                      duplicate = "mean", linear=FALSE,
                      xo=seq(min(df$x), max(df$x), length = 200),
                      yo=seq(min(df$y), max(df$y), length = 200))
persp3D(x=dens$x, y=dens$y, z=dens$z,ticktype="detailed",theta=230,
        ltheta=200,phi=30,lphi=30,xlab="exposure",ylab="lag",zlab="log-RR", zlim=c(-0.005,0.02),
        nticks = 4,cex.main = 2,
        shade = 0.75,r=sqrt(3),d=5,cex.axis=1.2, cex.lab=2,border=NA,
        col="steelblue", main = "Plane")

# scenario 2
dens <- akima::interp(x = df$x, 
                      y = df$y, 
                      z = trueeff[[2]], 
                      duplicate = "mean", linear= T,
                      xo=seq(min(df$x), max(df$x), length = 200),
                      yo=seq(min(df$y), max(df$y), length = 200))
persp3D(x=dens$x, y=dens$y, z=dens$z,ticktype="detailed",theta=230,
        ltheta=200,phi=30,lphi=30,xlab="exposure",ylab="lag",zlab="log-RR", zlim=c(-0.005,0.10),
        nticks = 4,cex.main = 2,
        shade = 0.75,r=sqrt(3),d=5,cex.axis=1.2, cex.lab=2,border=NA,
        col="steelblue", main = "Temp")


# scenario 3

dens <- akima::interp(x = df$x, 
                      y = df$y, 
                      z = trueeff[[3]], 
                      duplicate = "mean", linear=FALSE,
                      xo=seq(min(df$x), max(df$x), length = 200),
                      yo=seq(min(df$y), max(df$y), length = 200))
persp3D(x=dens$x, y=dens$y, z=dens$z, ticktype="detailed",theta=230,
        ltheta=200,phi=30,lphi=30,xlab="exposure",ylab="lag",zlab="log-RR", zlim=c(-0.001,0.025),
        nticks = 4,cex.main = 2,
        shade = 0.75,r=sqrt(3),d=5,cex.axis=1.2, cex.lab=2,border=NA,
        col="steelblue", main = "Complex")


dev.off()

# NUMBER OF ITERATIONS 
nsim <- 500
nsample <- 25
# BASELINE
base <- c(15,150,15)
# NOMINAL VALUE
qn <- qnorm(0.975)

L = 40
lagvar = 0:L
vx <- 9 # number of Bsplines for variable dimension
vl <- 10 # number of Bsplines for delay dimension
cen <- cen_ind[ind]
at_x <- seq(0,10,0.25)

# Store results
bias_Laplace <- matrix(0,ncol=dim(trueeff[[ind]])[2], nrow=dim(trueeff[[ind]])[1])
cov_Laplace <- matrix(0,ncol=dim(trueeff[[ind]])[2], nrow=dim(trueeff[[ind]])[1])
rmse_Laplace <- matrix(0,ncol=dim(trueeff[[ind]])[2], nrow=dim(trueeff[[ind]])[1])
cov_all_Laplace <- rep(0,dim(trueeff$Temperature)[1])
rmse_all_Laplace <- rep(0,dim(trueeff$Temperature)[1])
time_Laplace <- NULL
cov_mu_Laplace <- rmse_mu_Laplace <- numeric(sum(!is.na(cumeff)))

bias_gam <- matrix(0,ncol=dim(trueeff[[ind]])[2], nrow=dim(trueeff[[ind]])[1])
cov_gam <- matrix(0,ncol=dim(trueeff[[ind]])[2], nrow=dim(trueeff[[ind]])[1])
rmse_gam <- matrix(0,ncol=dim(trueeff[[ind]])[2], nrow=dim(trueeff[[ind]])[1])
cov_all_gam <- rep(0,dim(trueeff$Temperature)[1])
rmse_all_gam <- rep(0,dim(trueeff$Temperature)[1])
time_gam <- NULL
cov_mu_gam <- rmse_mu_gam <- numeric(sum(!is.na(cumeff)))


# Plots
pred_Laplace.meanx <- pred_gam.meanx <- rep(0,length(seq(0,10,0.25) ))
pred_Laplace.meanlag <- pred_gam.meanlag <- rep(0,L+1)
Laplace_x <- gam_x <- matrix(0,ncol=nsample, nrow=length(seq(0,10,0.25) ))
Laplace_lag <- gam_lag <- matrix(0,ncol=nsample, nrow=L+1)


for (i in 1:nsim){
  if (round(i/10)==i/10){print(i)}
  # print(i)
  
  set.seed(12805+i)
  
  suppressWarnings(y_all <- rpois(length(x),exp(log(base[ind])+cumeff)))
  mu <- exp(log(base[ind])+cumeff)
  
  # cor(na.omit(log(y_all+1)),na.omit(log(mu+1)))
  
  # Laplace
  crossbasis <- crossbasis(x, argvar=list(df=vx, fun="ps", intercept=F),
                           arglag=list(df=vl, fun="ps",intercept=T), lag=L)
  y_all[which(is.na(crossbasis[, 1]))] <- 0
  
  mtime <- proc.time()
  model_laplace <- DLNM_Laplace(y_all ~ 1,
                                crossbasis = crossbasis, 
                                vx = vx, vl = vl)
  time_Laplace[i] <- (proc.time()-mtime)[3]
  
  
  ####################################
  ################### Prediction #####
  pred_laplace <- predRR(model = model_laplace,
                         at_x = at_x,
                         cen = cen,
                         L = L)
  
  # STORE THE RESULTS
  bias_Laplace <- bias_Laplace + (matrix(pred_laplace$logpredX, ncol=L+1)-trueeff[[ind]])
  cov_Laplace <- cov_Laplace + (as.numeric(trueeff[[ind]]) >= pred_laplace$Qlower_logpredX &
                                  as.numeric(trueeff[[ind]]) <= pred_laplace$Qupper_logpredX)
  rmse_Laplace <- rmse_Laplace+(matrix(pred_laplace$logpredX, ncol = L+1)-trueeff[[ind]])^2
  
  cov_all_Laplace <- cov_all_Laplace + (as.numeric(apply(trueeff[[ind]],1,sum)) >= log(pred_laplace$Qlower_all) &
                                          as.numeric(apply(trueeff[[ind]],1,sum)) <= log(pred_laplace$Qupper_all))
  rmse_all_Laplace <- rmse_all_Laplace + (log(pred_laplace$pred_all) -as.numeric(apply(trueeff[[ind]],1,sum)))^2
  

  xvar <- 8
  xind <- which(seq(0,10,0.25)==xvar)
  pred_atxvar <- matrix(pred_laplace$logpredX, ncol = L+1)[xind,]
  pred_Laplace.meanlag <- pred_Laplace.meanlag + pred_atxvar
  pred_Laplace.meanx <- pred_Laplace.meanx + log(pred_laplace$pred_all)
  
  if(i<=nsample) {
    Laplace_x[,i] <- log(pred_laplace$pred_all)
    Laplace_lag[,i] <- pred_atxvar
  }
  
  # Predict outcome
  Xpred <- Matrix::Matrix(as.matrix(cbind(1,crossbasis[!is.na(crossbasis[,1]),])))
  
  mu_pred <-exp(as.numeric(Xpred%*%model_laplace$xi_mode))
  mu_true <- mu[!is.na(crossbasis[,1])]
  
  sd_mu <-  sqrt(pmax(0,Matrix::rowSums((Xpred%*%model_laplace$Sigma)*Xpred)))
  quantiles_mu <- mapply(function(mean_val, sd_val) {
    qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
  }, log(mu_pred), sd_mu)
  Qlower_mu = exp(quantiles_mu[1,])
  Qupper_mu = exp(quantiles_mu[2,])
  
  cov_mu_Laplace = cov_mu_Laplace + (mu_true >= Qlower_mu & mu_true <= Qupper_mu)
  rmse_mu_Laplace = rmse_mu_Laplace + (mu_true - mu_pred)^2
  
  
  
  
  # DEFINE THE PENALTY MATRICES
  # VARYING RIDGE PENALTY APPLIED TO COEFFICIENTS (Eq. 7a)
  Slag2 <-  diag((0:(vl-1))^2)
  
  cbPen <- cbPen(crossbasis,addSlag=list(Slag2)) 
  
  mtime <- proc.time()
  gam_try = tryCatch(
    {
    model_gam = model_gam<-bam(y_all ~ crossbasis,family=poisson(),
                               paraPen=list(crossbasis=cbPen))
    },
    error = function(err){
      model_gam<-bam(y_all ~ crossbasis,family=poisson(),
                     paraPen=list(crossbasis=cbPen), method="REML")
    }
  )


  time_gam[i] <- (proc.time()-mtime)[3]
  
  pred_gam <- crosspred(basis = crossbasis, model = model_gam, at = seq(0,10,0.25),
                        cen = cen, lag=L, bylag=1)
  # STORE THE RESULTS
  bias_gam <- bias_gam+ (pred_gam$matfit-trueeff[[ind]])
  cov_gam <- cov_gam + (trueeff[[ind]] >= pred_gam$matfit-qn*pred_gam$matse &
                          trueeff[[ind]] <= pred_gam$matfit+qn*pred_gam$matse)
  rmse_gam <- rmse_gam+(pred_gam$matfit-trueeff[[ind]])^2
  
  cov_all_gam <- cov_all_gam + (apply(trueeff[[ind]],1,sum) >= pred_gam$allfit-qn*pred_gam$allse &
                                  apply(trueeff[[ind]],1,sum) <= pred_gam$allfit+qn*pred_gam$allse)
  rmse_all_gam <- rmse_all_gam + (pred_gam$allfit-apply(trueeff[[ind]],1,sum))^2
  
  pred_delayvar <- pred_gam$allfit
  pred_gam.meanx <- pred_gam.meanx + pred_delayvar
  
  pred_atxvar <- pred_gam$matfit[xind,]
  pred_gam.meanlag <- pred_gam.meanlag + pred_atxvar
  
  if(i<=nsample) {
    gam_x[,i] <- pred_delayvar
    gam_lag[,i] <- pred_atxvar
  }
  
  # Predict outcome
  mu_pred <-exp(as.numeric(Xpred%*%coef(model_gam)))
  
  sd_mu <-  sqrt(pmax(0,Matrix::rowSums((Xpred%*%vcov(model_gam))*Xpred)))
  quantiles_mu <- mapply(function(mean_val, sd_val) {
    qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
  }, log(mu_pred), sd_mu)
  Qlower_mu = exp(quantiles_mu[1,])
  Qupper_mu = exp(quantiles_mu[2,])
  
  cov_mu_gam = cov_mu_gam + (mu_true >= Qlower_mu & mu_true <= Qupper_mu)
  rmse_mu_gam = rmse_mu_gam + (mu_true - mu_pred)^2
  
  
}

cor(log(y_all+1)[!is.na(crossbasis[,1])],log(predict(model_gam,type="response")+1))

results = data.frame(Metric = c("Bias", "Coverage", "Coverage all", "RMSE", "RMSE all",
                                "Coverage mu", "RMSE mu", "Time"),
                     Laplace = c(mean(bias_Laplace[seq(0,10,0.25)!=cen,]/nsim),
                                 mean((cov_Laplace[seq(0,10,0.25)!=cen,]/nsim)),
                                 mean((cov_all_Laplace[seq(0,10,0.25)!=cen]/nsim)),
                                 mean(sqrt(rmse_Laplace[seq(0,10,0.25)!=cen,]/nsim)),
                                 mean(sqrt(rmse_all_Laplace[seq(0,10,0.25)!=cen]/nsim)),
                                 mean(cov_mu_Laplace/nsim), mean(sqrt(rmse_mu_Laplace/nsim)),
                                 mean(time_Laplace)),
                     gam = c(mean(bias_gam[seq(0,10,0.25)!=cen,]/nsim),
                             mean((cov_gam[seq(0,10,0.25)!=cen,]/nsim)),
                             mean((cov_all_gam[seq(0,10,0.25)!=cen]/nsim)),
                             mean(sqrt(rmse_gam[seq(0,10,0.25)!=cen,]/nsim)),
                             mean(sqrt(rmse_all_gam[seq(0,10,0.25)!=cen]/nsim)),
                             mean(cov_mu_gam/nsim), mean(sqrt(rmse_mu_gam/nsim)),
                             mean(time_gam)))


print(results)


################################
#########################################
######################################

library(dplyr)
# Plot mean estimates
plot(0:L,Laplace_lag[,1],col=grey(0.8), type="l", ylim=c(-0.01,0.08),
     xlab="lag", ylab="log RR", main=paste("Laplace at x =", xvar))
for (m in 2:(nsample)){
  lines(0:L,Laplace_lag[,m],col=grey(0.8))
}
lines(0:L, trueeff[[ind]][xind,], col="red", lty=2)
lines(0:L, pred_Laplace.meanlag/nsim)

plot(0:L,gam_lag[,1],col=grey(0.8), type="l", ylim=c(-0.01,0.08),
     xlab="lag", ylab="log RR", main=paste("Gam at x =", xvar))
for (m in 2:(nsample)){
  lines(0:L,gam_lag[,m],col=grey(0.8))
}
lines(0:L, trueeff[[ind]][xind,], col="red", lty=2)
lines(0:L, pred_gam.meanlag/nsim)




# overall estimate
grid <- data.frame(x=rep(seq(0,10,0.25), each=41),lag=rep(0:40,41)) %>% mutate(result = do.call(combsim[ind], args=list(x,lag)))
trueeff <- matrix(grid$result, ncol=41, byrow=T)

plot(seq(0,10,0.25),Laplace_x[,1],col=grey(0.8), type="l", ylim=c(-0.3,0.8),
     xlab="var", ylab="log RR", main="Laplace overall risk")
for (m in 2:(nsample)){
  lines(seq(0,10,0.25),Laplace_x[,m],col=grey(0.8))
}
lines(seq(0,10,0.25), apply(trueeff,1,sum), col="red", lty=2)
lines(seq(0,10,0.25), pred_Laplace.meanx/nsim)

plot(seq(0,10,0.25),gam_x[,1],col=grey(0.8), type="l", ylim=c(-0.3,0.8),
     xlab="var", ylab="log RR", main="Gam overall risk")
for (m in 2:(nsample)){
  lines(seq(0,10,0.25),gam_x[,m],col=grey(0.8))
}
lines(seq(0,10,0.25), apply(trueeff,1,sum), col="red", lty=2)
lines(seq(0,10,0.25), pred_gam.meanx/nsim)
#########################################
######################################

