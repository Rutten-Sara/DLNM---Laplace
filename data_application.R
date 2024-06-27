rm(list = ls())
library(tidyverse)

setwd("G:/My Drive/Onderzoek/DLNM/DLNM Laplace/Final code")
source('DLNM_Laplace.R')
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

map_death = datafull %>% group_by(MSOA11CD, year) %>%
  summarize(tot_cases = sum(dtot)) %>%
  ungroup()

map_temp = datafull %>% group_by(MSOA11CD, year) %>%
  summarize(avg_temp = mean(tmean)) %>%
  ungroup()


map$death_2006 = (map_death %>% filter(year == 2006) %>% arrange(MSOA11CD))$tot_cases
map$death_2013 = (map_death %>% filter(year == 2013) %>% arrange(MSOA11CD))$tot_cases

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
tmap_arrange(map_2006, map_2013)

map_temp_2006 = tm_shape(map) +
  tm_polygons("temp_2006", palette = "YlOrRd", title="2006", breaks = c(16,17,18,19,20), style="cont") +
  tm_layout(legend.outside = TRUE, legend.outside.position = "bottom")
map_temp_2013 = tm_shape(map) +
  tm_polygons("temp_2013", palette = "YlOrRd", title="2013", breaks = c(16,17,18,19,20), style="cont") +
  tm_layout(legend.outside = TRUE, legend.outside.position = "bottom")
tmap_arrange(map_temp_2006, map_temp_2013)

################################################################################
#Prepare crossbasis matrix 

library(dlnm)
L <- 3 # maximum lag
vx <- 10 # number of basis for exposure var
vl <- 5 # number of basis for lag var
group <- factor(paste(datafull$MSOA11CD, datafull$year, sep="-"))
crossbasis <- crossbasis(datafull$tmean, lag=3, 
                         argvar=list(fun="ps",df = vx, intercept=F),
                         arglag=list(fun="ps",df = vl, intercept=T), group=group)
y_all = datafull$dtot
# DEFINE SPLINES OF DAY OF THE YEAR
spldoy <- onebasis(datafull$doy, "ns", df=3)
datafull$ID = as.factor(datafull$MSOA11CD) # random effect for every area

covar.rand <- c("ind",
           "ICAR",
           "Convolution",
           "Leroux")
tictoc::tic()
model_laplace <- DLNM_Laplace(y_all ~ spldoy:factor(year) + factor(dow),
                              crossbasis = crossbasis,
                              ID = ID, 
                              covar.ri = covar.rand[1], 
                              vx = vx, vl = vl,
                              data = datafull)
tictoc::toc()

#save(model_laplace, file = "Models/model_Leroux.RData")
library(mgcv)

# BAM

# DEFINE THE PENALTY MATRICES
Slag2 <-  diag((0:(vl-1))^2)

cbPen <- cbPen(crossbasis,addSlag=list(Slag2))

datafull$ID = as.factor(datafull$MSOA11CD) # random effect for every area

mtimegam <- proc.time()
model_gam = bam(y_all ~ crossbasis + s(ID, bs="re") + spldoy:factor(year) + factor(dow),
                data=datafull, family="poisson",
                paraPen=list(crossbasis=cbPen))
time_gam <- (proc.time()-mtimegam)[3]

save(model_gam, file = "model_gam.RData")
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
  tm_polygons("laplacespat", palette = "YlOrRd", title="Random intercept (Laplace)",
              breaks = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2)) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "bottom")

# Second map
map2 <- tm_shape(map) +
  tm_polygons("bamspat", palette = "YlOrRd", title="Random intercept (bam())",
              breaks = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2)) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "bottom")

# Arrange maps in a grid
tmap_arrange(map1, map2, ncol = 2)


#data_cor = data.frame(Leroux = as.matrix(exp(pred_laplace_Leroux$logpredX)), Convolution = as.matrix(exp(pred_laplace_BYM$logpredX)),
#                      ICAR = as.matrix(exp(pred_laplace_ICAR$logpredX)), Independent = as.matrix(exp(pred_laplace_ind$logpredX)))

#library(GGally)
#ggpairs(data_cor,columns = 1:4, 
#        title = "Correlation between different Laplace models", 
#       axisLabels = "show") 
