rm(list = ls())
library(tidyverse) ; library(readxl)

source('DLNM_Laplace_offset.R')
source('predRR.R')

# LOAD THE ORIGINAL MORTALITY DATA 
library(data.table)
# LOAD THE MORTALITY DATA
dataorig <- data.table::as.data.table(read.csv("data/lndmsoadeath.csv"))

# RENAME AND CREATE DATE
names(dataorig) <- c("year", "month", "day", "MSOA11CD", "MSOA11NM", "d074",
                     "d75plus")
dataorig[, date:=as.Date(paste(year, month, day, sep="/"))]

# DEFINE SERIES OF UNIQUE MSOA AND DATES (CHECK LATTER IS COMPLETE IN 2 SUMMERS)
seqmsoa <- sort(unique(dataorig$MSOA11CD))
seqdate <- sort(unique(dataorig$date))
table(diff(seqdate))

################################################################################
# PREPARE THE CASE TIME SERIES DATASET (STRATIFIED BY MSOA)

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

data_plot = datafull %>% group_by(date, year) %>% summarize(tot_cases = sum(dtot)) %>%
  ungroup() %>%
  arrange(date)

par(mfrow=c(1,2))
plot(data_plot$date[data_plot$year==2006], data_plot$tot_cases[data_plot$year==2006], type = "h",
     xlab = "Date", ylab = "number of cases", main = "2006")
plot(data_plot$date[data_plot$year==2013], data_plot$tot_cases[data_plot$year==2013], type = "h",
     xlab = "Date", ylab = "number of cases", main = "2013")
par(mfrow=c(1,1))

#######################################################
# Define temperature 
library(terra)
library(exactextractr)
unzip("data/lndmsoashp.zip")
map <- sf::st_read("lndmsoashp.shp")
# LOAD MIN/MAX GRIDDED TEMPERATURE DATA (TWO NETCDF FILES) AND COMPUTE TMEAN
lndtmingrid <- rast("data/lndtmingrid.nc")
lndtmaxgrid <- rast("data/lndtmaxgrid.nc")
lndtmeangrid <- (lndtmingrid + lndtmaxgrid) /2

# COMPUTE THE AREA-WEIGHTED AVERAGE OF CELLS INTERSECTING EACH MSOA
lndtmeanmsoa <- exact_extract(lndtmeangrid, map, fun="mean")
dimnames(lndtmeanmsoa) <- list(seqmsoa, as.character(seqdate))

# MERGE WITH MAIN DATASETS (EXPLOIT ORDERED SEQUENCES)
datafull$tmean <- c(t(lndtmeanmsoa))


################################################################################
# LOAD POPULATION DENSITY
dens_London <- read_excel("data/land-area-population-density-lsoa11-msoa11.xlsx") %>%
  rename("MSOA11CD" = "MSOA11 Code", "pop" = "Mid-2014 population") %>%
  dplyr:::select("MSOA11CD", "pop")
datafull = merge(datafull, dens_London, by = "MSOA11CD")

################################################################################

# Plots

map_death = datafull %>% group_by(MSOA11CD, year)  %>%
  summarize(tot_cases = sum(dtot), tot_pop = mean(pop)) %>%
  mutate(tot_inc = tot_cases/tot_pop*1000)%>%
  ungroup()

map_temp = datafull %>% group_by(MSOA11CD, year) %>%
  summarize(avg_temp = mean(tmean)) %>%
  ungroup()


map$death_2006 = (map_death %>% filter(year == 2006) %>% arrange(MSOA11CD))$tot_cases
map$death_2013 = (map_death %>% filter(year == 2013) %>% arrange(MSOA11CD))$tot_cases
map$inc_2006 = (map_death %>% filter(year == 2006) %>%
                  arrange(MSOA11CD))$tot_inc
map$inc_2013 = (map_death %>% filter(year == 2013) %>%
                  arrange(MSOA11CD))$tot_inc


map$temp_2006 = (map_temp %>% filter(year == 2006) %>% arrange(MSOA11CD))$avg_temp
map$temp_2013 = (map_temp %>% filter(year == 2013) %>% arrange(MSOA11CD))$avg_temp

library(tmap)
# First map

map_2006 = tm_shape(map) +
  tm_polygons("death_2006", palette = "YlOrRd", title="2006", breaks = c(0,10,20,30,40,50)) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "bottom")
map_2013 = tm_shape(map) +
  tm_polygons("death_2013", palette = "YlOrRd", title="2013", breaks = c(0,10,20,30,40,50)) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "bottom")
tmap_arrange(map_2006, map_2013, ncol=2)


map_inc_2006 = tm_shape(map) +
  tm_polygons("inc_2006", palette = "YlOrRd", title="2006", breaks = c(0,1,2,3,4,5,6)) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "bottom")
map_inc_2013 = tm_shape(map) +
  tm_polygons("inc_2013", palette = "YlOrRd", title="2013", breaks = c(0,1,2,3,4,5,6)) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "bottom")
tmap_arrange(map_2013, map_inc_2013, ncol=2)
tmap_arrange(map_inc_2006, map_inc_2013, ncol=2)

map_temp_2006 = tm_shape(map) +
  tm_polygons("temp_2006", palette = "YlOrRd", title="2006", breaks = c(16,17,18,19,20), style="cont") +
  tm_layout(legend.outside = TRUE, legend.outside.position = "bottom")
map_temp_2013 = tm_shape(map) +
  tm_polygons("temp_2013", palette = "YlOrRd", title="2013", breaks = c(16,17,18,19,20), style="cont") +
  tm_layout(legend.outside = TRUE, legend.outside.position = "bottom")
tmap_arrange(map_temp_2006, map_temp_2013, ncol=2)


################################################################################
#Prepare crossbasis matrix 

library(dlnm)
L <- 7 # maximum lag
vx <- 9 # number of basis for exposure var
vl <- 10 # number of basis for lag var
group <- factor(paste(datafull$MSOA11CD, datafull$year, sep="-"))
crossbasis <- crossbasis(datafull$tmean, lag=L, 
                         argvar=list(fun="ps",df = vx, intercept=F),
                         arglag=list(fun="ps",df = vl, intercept=T), group=group)

#datafull$dtot_more = (datafull$dtot+1)*2
y_all = datafull$dtot
# DEFINE SPLINES OF DAY OF THE YEAR
spldoy <- onebasis(datafull$doy, "ns", df=3)
ID = as.factor(datafull$MSOA11CD) # random effect for every area


covar.rand <- c("ind",
           "ICAR",
           "Convolution",
           "Leroux")
tictoc::tic()
model_laplace <- DLNM_Laplace(y_all ~ spldoy:factor(year) + factor(dow),
                              crossbasis = crossbasis,
                              ID = ID, 
                              covar.ri = covar.rand[4], offset = datafull$pop, 
                              vx = vx, vl = vl,
                              data = datafull) # for ICAR and Convolution:
                                               # scale offset by factor 1/1000 to avoid numerical issues
tictoc::toc()

#save(model_laplace, file = "Models/model_Leroux_offset.RData")
library(mgcv)

# BAM

# DEFINE THE PENALTY MATRICES
Slag2 <-  diag((0:(vl-1))^2)

cbPen <- cbPen(crossbasis,addSlag=list(Slag2))

mtimegam <- proc.time()
model_gam = bam(y_all ~ crossbasis + s(ID, bs="re") + spldoy:factor(year) + factor(dow) + offset(log(pop)),
                data=datafull, family="poisson",
                paraPen=list(crossbasis=cbPen))
time_gam <- (proc.time()-mtimegam)[3]

#save(model_gam, file = "Models/model_gam_offset.RData")
#load("Models/model_gam.RData")

####################################
################### Prediction #####
at_x = seq(9,27, by=0.5)
cen = 14
pred_laplace<- predRR(model = model_laplace,
                       at_x = at_x,
                       cen = cen, L = L)
pred_gam <- crosspred(crossbasis, model_gam, 
                      at = at_x , cen=14)

# PLOT
library(ggplot2)
col <- c("darkgoldenrod3", "aquamarine3")
parold <- par(no.readonly=T)
par(mar=c(4,4,1,0.5), las=1, mgp=c(2.5,1,0))
# Plot gam
plot(pred_gam, "overall", ylim=c(0.8,1.8), ylab="RR", col=col[1], lwd=1.5,
     xlab=expression(paste("Temperature ("*degree,"C)")), 
     ci.arg=list(col=alpha(col[1], 0.2)))

# Plot Laplace
plot.arg <- list(type = "l",col=col[2],  lwd=1.5)
fci <- function(x, high, low, ci.arg, plot.arg, noeff = NULL){
  polygon.arg <- modifyList(list(col = grey(0.9), border = NA), 
                            ci.arg)
  polygon.arg <- modifyList(polygon.arg, list(x = c(x, 
                                                    rev(x)), y = c(high, rev(low))))
  do.call(polygon, polygon.arg)
}
fci(x=at_x, high = pred_laplace$Qupper_all,
    low = pred_laplace$Qlower_all, ci.arg=list(col=alpha(col[2], 0.2)), plot.arg = plot.arg)
plot.arg <- modifyList(plot.arg, c(list(x = at_x, 
                                        y = pred_laplace$pred_all)))
do.call("lines",plot.arg)

legend("top", c("Bam", "Laplace"), lty=1, lwd=1.5, col=col, bty="n",
       inset=0.05, y.intersp=2, cex=0.8)
par(parold)


# Lag specific plot

lag_pred_Laplace= matrix(exp(pred_laplace$logpredX), ncol = L+1)
lag_pred_Laplace_lower= matrix(exp(pred_laplace$Qlower_logpredX), ncol = L+1)
lag_pred_Laplace_higher= matrix(exp(pred_laplace$Qupper_logpredX), ncol = L+1)

col <- c("darkgoldenrod3", "aquamarine3")
parold <- par(no.readonly=T)
par(mar=c(4,4,1,0.5), las=1, mgp=c(2.5,1,0))
# Plot gam
plot(pred_gam, "slices", ylim=c(0.9,1.25), ylab="RR", col=col[1], lwd=1.5,
     xlab=" Lag (days)", var = 24,
     ci.arg=list(col=alpha(col[1], 0.2)))

# Plot Laplace
plot.arg <- list(type = "l",col=col[2],  lwd=1.5)
fci(x=0:L, high = lag_pred_Laplace_higher[31,],
    low = lag_pred_Laplace_lower[31,], ci.arg=list(col=alpha(col[2], 0.2)), plot.arg = plot.arg)
plot.arg <- modifyList(plot.arg, c(list(x = 0:L, 
                                        y = lag_pred_Laplace[31,])))
do.call("lines",plot.arg)

legend("top", c("Bam", "Laplace"), lty=1, lwd=1.5, col=col, bty="n",
       inset=0.05, y.intersp=2, cex=0.8)
par(parold)


col <- c("darkgoldenrod3", "aquamarine3")
parold <- par(no.readonly=T)
par(mar=c(4,4,1,0.5), las=1, mgp=c(2.5,1,0))
# Plot gam
plot(pred_gam, "slices", ylim=c(0.9,1.25), ylab="RR", col=col[1], lwd=1.5,
     xlab=" Lag (days)", var = "20",
     ci.arg=list(col=alpha(col[1], 0.2)))

# Plot Laplace
plot.arg <- list(type = "l",col=col[2],  lwd=1.5)
fci(x=0:L, high = lag_pred_Laplace_higher[23,],
    low = lag_pred_Laplace_lower[23,], ci.arg=list(col=alpha(col[2], 0.2)), plot.arg = plot.arg)
plot.arg <- modifyList(plot.arg, c(list(x = 0:L, 
                                        y = lag_pred_Laplace[23,])))
do.call("lines",plot.arg)

legend("top", c("Bam", "Laplace"), lty=1, lwd=1.5, col=col, bty="n",
       inset=0.05, y.intersp=2, cex=0.8)
par(parold)


################################################
library("gratia")
gamspat <- smooth_coefs(model_gam, "s(ID)")
xispat <- model_laplace$xispat

plot(gamspat,xispat)
abline(a = 0, b = 1, col = "red", lwd = 3)

map$laplacespat <- xispat
map$bamspat <- gamspat

library(tmap)
# First map
map1 <- tm_shape(map) +
  tm_polygons("laplacespat", palette = "YlOrRd", title="Random intercept (Laplace)",midpoint = 0,
              breaks = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2)) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "bottom")

# Second map
map2 <- tm_shape(map) +
  tm_polygons("bamspat", palette = "YlOrRd", title="Random intercept (bam())",midpoint=0,
              breaks = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2)) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "bottom")

# Arrange maps in a grid
tmap_arrange(map1, map2, ncol = 2)



# Exceedance probability
at_x_exceedance = seq(9,27, by=0.01)
pred_laplace_exceedance<- predRR(model = model_laplace,
                      at_x = at_x_exceedance,
                      cen = cen, L = L)
pred_gam_exceedance <- crosspred(crossbasis, model_gam, 
                      at = at_x_exceedance , cen=14)

pred_all = pred_laplace_exceedance$pred_all
sd_all = pred_laplace_exceedance$sd_all


exceedance_1 <- mapply(function(mean_val, sd_val) {
  pnorm(q = 0, mean = mean_val, sd = sd_val, lower.tail = F)
}, log(pred_all), sd_all)

plot(at_x_exceedance, exceedance_1, type = "l", col = "blue", xlab = expression(paste("Temperature ("*degree,"C)")),
     ylab = "P(RR>1)")

#map$ex = exceedance_1
#tm_shape(map) +
#  tm_polygons("ex", palette = "YlOrRd", title="Exceedance probability")+
# tm_layout(legend.outside = TRUE, legend.outside.position = "bottom")


#q_high_risk = quantile(pred_all,0.95)
#exceedance_quantile <- mapply(function(mean_val, sd_val) {
#  pnorm(q = q_high_risk , mean = mean_val, sd = sd_val, lower.tail = F)
#}, log(pred_all), sd_all)

#map$high_risk = exceedance_quantile


################ Exceedance probability attributable fraction
source("af_Laplace.R")

temp_per_area = datafull %>%
  arrange(MSOA11CD, date)

ex_prob_2006 = NULL  
af_2006 = NULL
for (i in 1:length(unique(temp_per_area$MSOA11CD))){
  temp_per_area_i = temp_per_area %>% filter(MSOA11CD == unique(temp_per_area$MSOA11CD)[i] & year==2006)
  #ID = temp_per_area_i$year
  af_i <- attrdl_Laplace(temp_per_area_i$tmean,crossbasis,temp_per_area_i$dtot,
                               xi_mode=model_laplace$xi_mode[model_laplace$ind.cb], 
                               sigma_xi=model_laplace$Sigma[model_laplace$ind.cb, model_laplace$ind.cb],
                               type="af",dir="back",tot=TRUE,cen = 14,
                               sim=T, group=NULL, ex.prob = 0)
  ex_prob_2006[i] <- af_i$probs
  af_2006[i] = af_i$res

}

map_exceedance_2006 <- map %>% arrange(MSOA11CD)
map_exceedance_2006$prob = ex_prob_2006
ex_2006 <- tm_shape(map_exceedance_2006) +
  tm_polygons("prob", palette = "Blues", title="P(af>0) in 2006",
              breaks = c(-Inf,0.9,0.95,1)) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "bottom")



ex_prob_2013 = NULL  
af_2013 = NULL
for (i in 1:length(unique(temp_per_area$MSOA11CD))){
  temp_per_area_i = temp_per_area %>% filter(MSOA11CD == unique(temp_per_area$MSOA11CD)[i] & year==2013)
  #ID = temp_per_area_i$year
  af_i <- attrdl_Laplace(temp_per_area_i$tmean,crossbasis,temp_per_area_i$dtot,
                         xi_mode=model_laplace$xi_mode[model_laplace$ind.cb], 
                         sigma_xi=model_laplace$Sigma[model_laplace$ind.cb, model_laplace$ind.cb],
                         type="af",dir="back",tot=TRUE,cen = 14,
                         sim=T, group=NULL, ex.prob = 0)
  ex_prob_2013[i] <- af_i$probs
  af_2013[i] = af_i$res
  
}

map_exceedance_2013<- map %>% arrange(MSOA11CD)
map_exceedance_2013$prob = ex_prob_2013
ex_2013 <- tm_shape(map_exceedance_2013) +
  tm_polygons("prob", palette = "Blues", title="P(af>0) in 2013",
              breaks = c(-Inf,0.9,0.95,1)) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "bottom")


tmap_arrange(ex_2006, ex_2013, ncol = 2)




################### Correlations between different models

load("Models/model_Leroux_offset.RData")
pred_laplace_Leroux<- predRR(model = model_laplace,
                      at_x = at_x,
                      cen = cen, L = L)


load("Models/model_BYM_offset.RData")
pred_laplace_BYM<- predRR(model = model_laplace,
                             at_x = at_x,
                             cen = cen, L = L)


load("Models/model_ICAR_offset.RData")
pred_laplace_ICAR<- predRR(model = model_laplace,
                             at_x = at_x,
                             cen = cen, L = L)


load("Models/model_ind_offset.RData")
pred_laplace_ind<- predRR(model = model_laplace,
                             at_x = at_x,
                             cen = cen, L = L)


data_cor = data.frame(Leroux = as.matrix(exp(pred_laplace_Leroux$logpredX)), Convolution = as.matrix(exp(pred_laplace_BYM$logpredX)),
                      ICAR = as.matrix(exp(pred_laplace_ICAR$logpredX)), Independent = as.matrix(exp(pred_laplace_ind$logpredX)))

library(GGally)
ggpairs(data_cor,columns = 1:4, 
        title = "Correlation between different Laplace models", 
       axisLabels = "show") 



