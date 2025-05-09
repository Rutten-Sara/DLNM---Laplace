rm(list = ls())
library(splines); library(dlnm); library(tsModel); library(biglm);  library(MASS); library(tidyverse); library(plot3D)
library(mgcv); library(nlraa); library(R.utils)
library(Matrix)
library(data.table)

################################################################################
# PREPARE TIME SERIES DATA
################################################################################
source('DLNM_Laplace_offset.R')
source('predRR.R')


# LOAD THE MORTALITY DATA
dataorig <- as.data.table(read.csv("data/lndmsoadeath.csv")) 

# RENAME AND CREATE DATE
names(dataorig) <- c("year", "month", "day", "MSOA11CD", "MSOA11NM", "d074",
                     "d75plus")
dataorig[, date:=as.Date(paste(year, month, day, sep="/"))]

# DEFINE SERIES OF UNIQUE MSOA AND DATES (CHECK LATTER IS COMPLETE IN 2 SUMMERS)
seqmsoa <- sort(unique(dataorig$MSOA11CD))
seqdate <- sort(unique(dataorig$date))
table(diff(seqdate))

# COMPLETE THE ORIGINAL SERIES (INCLUDING DATES WITH NO DEATH)
datafull <- expand.grid(MSOA11CD=seqmsoa, date=seqdate) |>
  data.table() |>
  merge(dataorig[,c(4, 6:8)], all.x=T) |>
  merge(unique(dataorig[,c("MSOA11CD", "MSOA11NM")]), by="MSOA11CD")

# PAD WITH 0 WHEN NO DEATH
datafull[is.na(datafull)] <- 0

# CREATE TOTAL DEATHS, RE-CREATE TIME VARS
datafull[, dtot:=d074+d75plus]
datafull[, `:=`(year=year(date), month=month(date), day=mday(date),
                doy=yday(date), dow=wday(date))]

# ORDER (IMPORTANT FOR KEEPING THE TIME SERIES SEQUENCE BY MSOA)
setkey(datafull, MSOA11CD, date)



library(sf);library(tmap); library(dplyr); library(SpatialEpi)

unzip("data/lndmsoashp.zip")
map <- st_read("lndmsoashp.shp")

total_per_area = datafull %>% group_by(MSOA11CD) %>% summarize(tot_death = sum(dtot)) %>%
  arrange(MSOA11CD)

map$tot_death = total_per_area$tot_death


# Define temperature 
library(terra)
library(exactextractr)
# LOAD MIN/MAX GRIDDED TEMPERATURE DATA (TWO NETCDF FILES) AND COMPUTE TMEAN
lndtmingrid <- rast("data/lndtmingrid.nc")
lndtmaxgrid <- rast("data/lndtmaxgrid.nc")
lndtmeangrid <- (lndtmingrid + lndtmaxgrid) /2


# COMPUTE THE AREA-WEIGHTED AVERAGE OF CELLS INTERSECTING EACH MSOA
lndtmeanmsoa <- exact_extract(lndtmeangrid, map, fun="mean")
dimnames(lndtmeanmsoa) <- list(seqmsoa, as.character(seqdate))

# MERGE WITH MAIN DATASETS (EXPLOIT ORDERED SEQUENCES)
datafull$tmean <- c(t(lndtmeanmsoa))


# DEFINE THE EXPOSURE

# REAL TEMPERATURE SERIES, STANDARDIZED IN 0-10
x <- datafull$tmean
x <- (x-min(x))/diff(range(x))*10

# MATRIX Q OF EXPOSURE HISTORIES
group <- factor(paste(datafull$MSOA11CD, datafull$year, sep="-"))
Q <- Lag(x,0:40, group=group)

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


ind = 3
cumeff <- apply(Q,1,fcumeff,0:40,combsim[ind])

# NUMBER OF ITERATIONS 
nsim <- 10
nsample <- 10

# BASELINE
base <- c(90,270,130)
base <- c(0.01,0.01,0.01)

# NOMINAL VALUE
qn <- qnorm(0.975)

L = 40
vx <- 9 # number of Bsplines for variable dimension
vl <- 10 # number of Bsplines for delay dimension
cen <- cen_ind[ind]
at_x <- seq(0,10,0.25)
datafull$ID = as.factor(datafull$MSOA11CD)


# Store results
bias_Laplace <- matrix(0,ncol=dim(trueeff$Temperature)[2], nrow=dim(trueeff$Temperature)[1])
cov_Laplace <- matrix(0,ncol=dim(trueeff$Temperature)[2], nrow=dim(trueeff$Temperature)[1])
rmse_Laplace <- matrix(0,ncol=dim(trueeff$Temperature)[2], nrow=dim(trueeff$Temperature)[1])
cov_all_Laplace <- numeric(dim(trueeff$Temperature)[1])
rmse_all_Laplace <- numeric(dim(trueeff$Temperature)[1])
time_Laplace <- rep(NA, nsim)

bias_gam <- matrix(0,ncol=dim(trueeff$Temperature)[2], nrow=dim(trueeff$Temperature)[1])
cov_gam <- matrix(0,ncol=dim(trueeff$Temperature)[2], nrow=dim(trueeff$Temperature)[1])
rmse_gam <- matrix(0,ncol=dim(trueeff[[ind]])[2], nrow=dim(trueeff[[ind]])[1])
cov_all_gam <- numeric(dim(trueeff$Temperature)[1])
rmse_gam <- matrix(0,ncol=dim(trueeff$Temperature)[2], nrow=dim(trueeff$Temperature)[1])
rmse_all_gam <- numeric(dim(trueeff$Temperature)[1])
time_gam <- rep(NA,nsim)


cov_re_Laplace <- rmse_re_Laplace <- numeric(length(unique(datafull$ID)))
cov_mu_Laplace <- rmse_mu_Laplace <- numeric(sum(!is.na(cumeff)))
cov_re_gam <- rmse_re_gam <- numeric(length(unique(datafull$ID)))
cov_mu_gam <- rmse_mu_gam <- numeric(sum(!is.na(cumeff)))


# Plots
pred_Laplace.meanx <- pred_gam.meanx <- numeric(length(seq(0,10,0.25) ))
pred_Laplace.meanlag <- pred_gam.meanlag <- numeric(L+1)
Laplace_x <- gam_x <- matrix(0,ncol=nsample, nrow=length(seq(0,10,0.25) ))
Laplace_lag <- gam_lag <- matrix(0,ncol=nsample, nrow=L+1)

spatial_effect <- numeric(nsim)

set.seed(1)
offset_true = ceiling(rgamma(length(unique(datafull$MSOA11CD)),1,0.05))

for (i in 1:nsim){
  if (round(i/10)==i/10){print(i)}
  print(i)
  
  set.seed(12805+i)
  
  unique_area = data.frame(area = unique(datafull$MSOA11CD), offset = offset_true)
  random_effect <- rnorm(length(unique_area$area),0,sqrt(0.5)) 
  unique_area$random_effect = random_effect
  
  random_area = inner_join(data.frame(area = datafull$MSOA11CD), unique_area,
                            by = "area")$random_effect
  
  offset = inner_join(data.frame(area = datafull$MSOA11CD), unique_area,
                      by = "area")$offset
  
  
  suppressWarnings(y_all <- rpois(length(x),offset*exp(log(base[ind])+cumeff+random_area))) 
  mu <- exp(log(base[ind])+cumeff+random_area)

  

  crossbasis <- crossbasis(x, argvar=list(df=vx, fun="ps", intercept=F),
                           arglag=list(df=vl, fun="ps",intercept=T), lag=L, group=group)
  
  # Fit bam model
  Slag2 <-  diag((0:(vl-1))^2)
  cbPen <- cbPen(crossbasis,addSlag=list(Slag2)) 
  
  ID_gam = as.factor(datafull$ID)
  
  mtime <- proc.time()
  
  gam_try = tryCatch(
    {
      model_gam<-withTimeout(
        bam(y_all ~ crossbasis+ s(ID_gam, bs="re")+offset(log(offset)),family=poisson(),
            paraPen=list(crossbasis=cbPen)),
        timeout = 30*60,
        onTimeout = "error"
      )
    },
    error = function(err){
      model_gam<-withTimeout(
        bam(y_all ~ crossbasis+ s(ID_gam, bs="re")+offset(log(offset)),family=poisson(),
            paraPen=list(crossbasis=cbPen), method = "REML"),
        timeout = 30*60,
        onTimeout = "silent"
      )
    }
  )
  
  if(is.null(model_gam)){
    next
  }
  
  time_gam[i] <- (proc.time()-mtime)[3]
  
  # Laplace
  y_all[which(is.na(crossbasis[, 1]))] <- 0
  
  mtime <- proc.time()
  model_laplace <- DLNM_Laplace(y_all ~ 1, 
                                crossbasis = crossbasis, 
                                ID = as.factor(datafull$ID),
                                covar.ri = "ind", offset = offset,
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
  Z.rand <- Matrix::sparse.model.matrix(~ as.factor(datafull$ID) + 0)[!is.na(crossbasis[,1]),] 
  Xpred <- Matrix::Matrix(as.matrix(cbind(1,crossbasis[!is.na(crossbasis[,1]),],Z.rand)))
  
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
  
  
  # Predict random effects
  ind_re <- (length(model_laplace$xi_mode)-length(unique(datafull$ID))+1) :(length(model_laplace$xi_mode))
  re <- as.numeric(Z.rand %*% model_laplace$xi_mode[ind_re])
  re_true <- random_area[!is.na(crossbasis[,1])]
  
  sd_re <- sqrt(pmax(0,Matrix::rowSums((Z.rand%*%model_laplace$Sigma[ind_re,ind_re])*Z.rand)))
  
  quantiles_re <- mapply(function(mean_val, sd_val) {
    qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
  }, re, sd_re)
  Qlower_re = quantiles_re[1,]
  Qupper_re = quantiles_re[2,]
  
  cov_re_Laplace <- cov_re_Laplace + (unique(re_true) >= unique(Qlower_re) & unique(re_true) <= unique(Qupper_re))
  rmse_re_Laplace <- rmse_re_Laplace + (unique(re_true) - unique(re))^2
  
  
  # Bam
  pred_gam <- crosspred(basis = crossbasis, model=model_gam, at = seq(0,10,0.25),
                        cen = cen, lag=L, bylag=1)
  # STORE THE RESULTS
  bias_gam <- bias_gam + (pred_gam$matfit-trueeff[[ind]])
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
  
  
  # Predict random effects
  ind_re <- grepl( "ID", names(coef(model_gam)))
  re <- as.numeric(Z.rand %*% coef(model_gam)[ind_re])
  sd_re <- sqrt(pmax(0,Matrix::rowSums((Z.rand%*%vcov(model_gam)[ind_re,ind_re])*Z.rand)))
  
  quantiles_re <- mapply(function(mean_val, sd_val) {
    qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
  }, re, sd_re)
  Qlower_re = quantiles_re[1,]
  Qupper_re = quantiles_re[2,]
  
  cov_re_gam <- cov_re_gam + (unique(re_true) >= unique(Qlower_re) & unique(re_true) <= unique(Qupper_re))
  rmse_re_gam <- rmse_re_gam + (unique(re_true) - unique(re))^2
  
  
}

#cor(na.omit(log(y_all+1)),log(predict(model_gam,type="response")+1))


results = data.frame(Metric = c("Bias", "Coverage", "Coverage all", "RMSE", "RMSE all",
                                "Coverage mu", "RMSE mu", "Coverage re", "RMSE re", "Time"),
                     Laplace = c(mean(bias_Laplace[seq(0,10,0.25)!=cen,]/nsim),
                                 mean((cov_Laplace[seq(0,10,0.25)!=cen,]/nsim)),
                                 mean((cov_all_Laplace[seq(0,10,0.25)!=cen]/nsim)),
                                 mean(sqrt(rmse_Laplace[seq(0,10,0.25)!=cen,]/nsim)),
                                 mean(sqrt(rmse_all_Laplace[seq(0,10,0.25)!=cen]/nsim)),
                                 mean(cov_mu_Laplace/nsim), mean(sqrt(rmse_mu_Laplace/nsim)),
                                 mean(cov_re_Laplace/nsim), mean(sqrt(rmse_re_Laplace/nsim)),
                                 mean(time_Laplace)),
                     gam = c(mean(bias_gam[seq(0,10,0.25)!=cen,]/nsim),
                             mean((cov_gam[seq(0,10,0.25)!=cen,]/nsim)),
                             mean((cov_all_gam[seq(0,10,0.25)!=cen]/nsim)),
                             mean(sqrt(rmse_gam[seq(0,10,0.25)!=cen,]/nsim)),
                             mean(sqrt(rmse_all_gam[seq(0,10,0.25)!=cen]/nsim)),
                             mean(cov_mu_gam/nsim), mean(sqrt(rmse_mu_gam/nsim)),
                             mean(cov_re_gam/nsim), mean(sqrt(rmse_re_gam/nsim)),
                             mean(time_gam)))
print(results)


#########################################
########################################
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
