rm(list=ls())
library(dplyr)
library(ggplot2)
library(sp)
library(useful)
library(gstat)
library(spacetime)
library(raster)
library(doParallel)
library(rgdal)
library(rgeos)
library(readr)
library(GSIF)
library(clifro)
library(MLmetrics)
library(mgcv)
library(flexmix)
theme_set(theme_minimal())

#######################################
# READ DATA
#######################################
datos=read.table("output_14_10_2016_07_44.csv",
                 header = TRUE,dec=".",sep=";")
datos=datos[c("Time", "ICAO", "Lat", "Long", "Alt", "TAS",
              "GSP", "Course", "Heading", "Mach", "WindSpeed",
              "WindDir", "Temperature")]

datos$Time=strptime(as.character(datos$Time),"%Y-%m-%d_%H.%M.%S")

datos=datos[order(datos$Time),]
datos$dtm=round((-as.numeric(difftime(datos$Time[1],datos$Time)))/60,3)

int.ws=round(quantile(datos$WindSpeed,p=c(1/4,0.5,3/4)),0)
datos$speed=factor(findInterval(datos$WindSpeed,int.ws),labels =
                     c("f1","fmed","f3","max"))

int.dir=seq(0,360,by=90)
datos$dir=factor(findInterval(datos$WindDir,int.dir),labels =
                   as.character(paste("deg",int.dir[-1],sep="")))

# datos=(datos[c("Time","WindSpeed","speed", "dir","Lat", "Long", "Alt")])
hist(datos$WindSpeed)

# log transformation to achieve normality
hist(log(datos$WindSpeed))

shapiro.test(log(datos$WindSpeed))

# North and East Wind components:
datos$North_comp = log(datos$WindSpeed)*cos(datos$WindDir*pi/180)
datos$East_comp = log(datos$WindSpeed)*sin(datos$WindDir*pi/180)

#########################
# EXPLORATORY ANALYSIS
#########################

# Little data in first and second quarter
findInterval(datos$WindDir,int.dir)
table(findInterval(datos$WindDir,int.dir))

# Lets see the data on a map:
library(rworldmap)
library(ggplot2)
library("mapproj")
data(countriesLow)

world <- fortify(countriesLow)
## this converts any spatialobjectdataframe to a dataframe to be plotted in ggplot

map <- ggplot() +
  geom_polygon(data = world,
               aes(x=long, y=lat, group=group),
               color = "black", fill = "lightgreen") +
  geom_point(data = datos, 
             mapping = aes(x = Long, y = Lat))

xxlim <- c(-8,4)   ## selected range North Sea
yylim <- c(36,44)

map +
  coord_cartesian(xlim = xxlim, ylim = yylim)

map +
  coord_map("ortho", xlim = xxlim, ylim = yylim, orientation=c(55, 10, 0))

map +
  coord_map("stereographic", xlim = xxlim, ylim = yylim, orientation=c(55, 10, 0)) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

library(plotly)
# 3D-Scatterplot
plot_ly(as.data.frame(datos),x=~datos$Lat,y=~datos$Long,z=~datos$Alt, type="scatter3d", mode="markers")


# Quiver plot of the vectorial field
arrowtips = datos[c(3,4,5)]
arrowtips[1] = arrowtips[1]+cos(datos$WindDir*pi/180)
arrowtips[2] = arrowtips[2]+sin(datos$WindDir*pi/180)

datos$color=factor(findInterval(datos$WindSpeed,int.ws),labels =
                     c("green","blue","red","black"))
arrows3D(x0=as.matrix(datos[3]), y0=as.matrix(datos[4]), z0=as.matrix(datos[5]), x1=as.matrix(arrowtips[1]), y1 = as.matrix(arrowtips[2]), z1 = as.matrix(arrowtips[3]),  lwd = 2, d = 3, 
         main = "Wind speed and direction", bty ="g", ticktype = "detailed", col=datos$color, xlab = "Latitude", ylab= "Longitude", zlab="Altitude")


# Windrose diagram
library(clifro)
windrose(speed = log(datos$WindSpeed),
         direction = datos$WindDir,
         ggtheme='minimal', legend_title = "Log of Wind Speed")


#   Any trends?
par(mfrow=c(2,2))
plot(datos$North_comp, datos$Lat, xlab="North velocity (log kt)", ylab="Latitude (º)", main="North velocity vs Latitude")
plot(datos$North_comp, datos$Long, xlab="North velocity (log kt)", ylab="Longitude (º)", main="North velocity vs Longitude")
plot(datos$North_comp, datos$Alt, xlab="North velocity (log kt)", ylab="Alt (ft)", main="North velocity vs Alt")
ts.plot(datos$North_comp,  xlab="Time(s)", ylab="North velocity (log kt)", main="Evolution of North velocity")
par(mfrow=c(1,1))


par(mfrow=c(2,2))
plot(datos$East_comp, datos$Lat, xlab="East velocity (log kt)", ylab="Latitude (º)", main="East velocity vs Latitude")
plot(datos$East_comp, datos$Long, xlab="East velocity (log kt)", ylab="Longitude (º)", main="East velocity vs Longitude")
plot(datos$East_comp, datos$Alt, xlab="East velocity (log kt)", ylab="Alt (ft)", main="East velocity vs Alt")
ts.plot(datos$East_comp,  xlab="Time(s)", ylab="East velocity (log kt)", main="Evolution of East velocity")

#################################################################
# REGRESSION KRIGING 
#################################################################

# We are going to fit a trend to the data, and then perform ordinary kriging 
# of the residuals to predict the residuals at an unknown locations.
# Then the predicted residuals are added back to the trend



### Lets try using the first 4 minutes to predict the next minutes
# data with the first 4 minutes:
train_interval = datos$Time>=as.POSIXct('2016-10-14 07:44:45 CET')&datos$Time<=as.POSIXct('2016-10-14 07:47:59 CET')
sub1 <- datos[train_interval,]

# Create a SpatialPointsDataFrame. 
# transform the object sub into a spatial object, and then I changed its projection into UTM so that the variogram 
# will be calculated on meters and not degrees
coordinates(sub1)=~Lat+Long+Alt
projection(sub1)=CRS("+init=epsg:4326")

# ~ 7 minutes to be predicted:
time_interval = datos$Time>=as.POSIXct('2016-10-14 07:48:00 CET')&datos$Time<=as.POSIXct('2016-10-14 07:54:44 CET')


# We can try different models to fit a trend to the data


# Generalized Additive Model with cubic splines
gam_model <- gam(North_comp~s(Lat,bs="cs",m=3,k=25) + 
                   s(Long,bs="cs",m=2,k=20) + 
                   s(Alt,bs="cs",m=2,k=30), data = sub1,
                 family=gaussian(), select = T)

# Generalized Additive Mixed model with gaussian splines and cubic splines
# gam_model = gam(North_comp ~ s(Lat, Long, bs='gp',m=2, k=120) + s(dtm,bs="cs",m=2,k=30) ,family=gaussian, data=sub1)

# Generalized Additive mixed model with splines on a sphere
# gam_model = gam(North_comp ~ s(Lat, Long, bs='sos',m=2, k=120) + s(Alt, dtm, bs='sos', k=115),family=gaussian, data=sub1)

# gam_model  =  gam(North_comp ~ s(Lat, Long, bs='sos',m=2, k=120), data=sub1)


summary(gam_model)
par(mfrow=c(2,2))
gam.check(gam_model)
concurvity(gam_model)
AIC(gam_model)
par(mfrow=c(1,1))

# Get the residuals
north_residuals = gam_model$residuals
hist(north_residuals)

# KRIGING OF THE RESIDUALS

#Transform into Mercator Projection
Wind.UTM <- spTransform(sub1,CRS("+init=epsg:3395")) 


WindSP <- SpatialPoints(Wind.UTM@coords,CRS("+init=epsg:3395"))

North_DF <- data.frame(north_residuals)

WindTM <- as.POSIXct(Wind.UTM$Time,tz="CET")

NorthSTIDF <- STIDF(WindSP,WindTM,data=North_DF)

stplot(NorthSTIDF)

no_cores <- detectCores()

# spatio-temporal variogram of North Wind.
North_var <- variogramST(north_residuals~Alt,data=NorthSTIDF,assumeRegular=T,na.omit=T, cores=no_cores, tunit='secs')

plot(North_var ,map=F)
plot(North_var ,wireframe=T)
plot(North_var ,map=T)


# To save time we can load the file next times
write.csv(North_var,"North_var_sec_MS2.csv", row.names = TRUE)
# North_var <- read.csv("North_var_sec_MS2.csv")
# North_var = North_var[-1]
# class(North_var) = c("StVariogram", "data.frame")

##########################################################
# Trying to fit a spatio-temporal variogram model
# to spatio-temporal sample variogram:

####
#SEPARABLE MODEL
separable <- vgmST("separable", space = vgm(1,"Sph", 200000, 0),time = vgm(1,"Sph", 50, 1), sill=1) 

separable_Vgm <- fit.StVariogram(North_var, separable, fit.method=6, method="L-BFGS-B", tunit="secs")
attr(separable_Vgm,"MSE")

plot(North_var, separable_Vgm, wireframe=T, all=T)
plot(North_var, separable_Vgm, map=F, all=T)
####

###
#PRODUCT SUM MODEL
prodSumModel <- vgmST("productSum",space = vgm(1, "Gau", 200000, 0),time = vgm(1, "Gau", 50, 0),k = 30) 
prodSumModel_Vgm <- fit.StVariogram(North_var, prodSumModel,method = "L-BFGS-B", fit.method=6)
attr(prodSumModel_Vgm, "MSE")

plot(North_var, prodSumModel_Vgm, wireframe=T, all=T)
plot(North_var, prodSumModel_Vgm, map=F, all=T)
###

####
# METRIC MODEL
metric <- vgmST("metric", joint = vgm(1,"Gau", 200000, 0), stAni=0.1) 
metric_Vgm <- fit.StVariogram(North_var, metric, method="L-BFGS-B", fit.method = 6)
attr(metric_Vgm, "MSE")

plot(North_var, metric_Vgm, wireframe=T, all=T)
plot(North_var, metric_Vgm, map=F, all=T)
####

####
# SUM-METRIC MODEL
sumMetric <- vgmST("sumMetric", space = vgm(psill=1,"Gau", range=250000, nugget=0.1),time = vgm(psill=0.5,"Gau", range=50, nugget=0), joint = vgm(psill=0,"Sph", range=200000, nugget=0.01), stAni=10) 

sumMetric_Vgm <- fit.StVariogram(North_var, sumMetric, method="L-BFGS-B",tunit="secs",  fit.method = 6)
attr(sumMetric_Vgm, "MSE")

plot(North_var, sumMetric, wireframe=T, all=T)
plot(North_var, sumMetric_Vgm, wireframe=T, all=T)
plot(North_var, sumMetric_Vgm, wireframe=T, diff=TRUE)
plot(North_var, sumMetric_Vgm, map=F, all=T)

####
# SIMPLE SUM-METRIC MODEL
SimplesumMetric <- vgmST("simpleSumMetric",space = vgm(1,"Gau", 250000, 0),time = vgm(0.5,"Gau", 50, 0), joint = vgm(0,"Sph", 200000, 0.1), nugget=0.1, stAni=0.1) 
SimplesumMetric_Vgm <- fit.StVariogram(North_var, SimplesumMetric,method = "L-BFGS-B")
attr(SimplesumMetric_Vgm, "MSE")

plot(North_var, SimplesumMetric_Vgm, wireframe=T, all=T)
plot(North_var, SimplesumMetric_Vgm, map=F, all=T)



plot(North_var,list(separable_Vgm, prodSumModel_Vgm, metric_Vgm ,sumMetric_Vgm, SimplesumMetric_Vgm), wireframe=T, all=T)

# SPATIO-TEMPORAL PREDICTION GRID

time_interval = datos$Time>=as.POSIXct('2016-10-14 07:48:00 CET')&datos$Time<=as.POSIXct('2016-10-14 07:48:10 CET')

time.grid = unique(datos$Time[time_interval,])

grid = datos[time_interval, c(3,4,5)]

coordinates(grid)=~Lat+Long+Alt
projection(grid)=CRS("+init=epsg:3395")
bbox(grid)

attr(separable_Vgm , "temporal unit") <- units(abs(outer(index(NorthSTIDF@time[1]), index(NorthSTIDF@time[1]), "-")))
attr(sumMetric_Vgm, "temporal unit") <- units(abs(outer(index(NorthSTIDF@time[1]), index(NorthSTIDF@time[1]), "-")))


# grid.stf = STF(grid, time.grid)
# pred = krigeST(north_residuals~1, data=NorthSTIDF, modelList=sumMetric_Vgm, newdata=grid.stf, computeVar = T, nmax=15) 
# stplot(pred)


grid.sti = STI(grid, datos$Time[time_interval,])

# PREDICTING THE RESIDUALS WITH KRIGING

# pred = krigeST(north_residuals~1, data=NorthSTIDF, modelList=separable_Vgm, newdata=grid.sti, computeVar = T, stAni = 10)
# pred = krigeST(north_residuals~1, data=NorthSTIDF, modelList=sumMetric_Vgm, newdata=grid.sti, nmax=450, stAni=0.0001) 
pred = krigeST(north_residuals~1, data=NorthSTIDF, modelList=sumMetric_Vgm, newdata=grid.sti, stAni=0.0001) 
plot(as.numeric(unlist(pred@data[1])), datos$North_comp[time_interval])


# SPEED PREDICTIONS WITH THE TREND
trend_pred = predict(gam_model, datos[time_interval, c(3,4,5, 14)])

# ADDING PREDICTED RESIDUALS TO THE TREND
prediction = trend_pred + pred@data[1]

library(MLmetrics)
MSE(prediction$var1.pred,datos$North_comp[time_interval])
MSE(trend_pred,datos$North_comp[time_interval])

aux = which(prediction$var1.pred>-6 & prediction$var1.pred < 4)
plot(prediction$var1.pred[aux], datos$North_comp[time_interval][aux], xlab="Predicted values of north velocity", ylab="Observed values of north velocity",main="Predicted vs Observed values of north velocity")
abline(0, 1, col=2)
plot(trend_pred[aux], datos$North_comp[time_interval][aux], xlab="Predicted values of north velocity with the trend", ylab="Observed values of north velocity",main="Predicted vs Observed values of north velocity")
abline(0, 1, col=2)

MSE(prediction$var1.pred[aux], datos$North_comp[time_interval][aux])
MSE(trend_pred[aux], datos$North_comp[time_interval][aux])

plot(prediction$var1.pred, datos$North_comp[time_interval],  xlab="Predicted values of north velocity", ylab="Observed values of north velocity",main="Predicted vs Observed values of north velocity")
abline(0,1,col=2)
plot(trend_pred, datos$North_comp[time_interval])
abline(0,1,col=2)


#### East component ####

# We repeat the same analysis for the east component of the wind


# gam_model2 <- gam(East_comp~s(Lat,bs="cs",m=3,k=30) + 
#                     s(Long,bs="cs",m=3,k=15) + 
#                     s(Alt,bs="cs",m=3,k=26), data = sub1,
#                   family=gaussian(), select = T
# )

gam_model2 <- gam(East_comp~s(Lat,bs="cs",m=3,k=30) + 
                    s(Long,bs="cs",m=3,k=20) + 
                    s(Alt,bs="cs",m=3,k=26), data = sub1,
                  family=gaussian(), select = T
)

# gam_model2 = gam(East_comp ~ s(Long, Lat, bs='sos', k=90) +
#                   s(dtm, Alt, bs='sos', k=90) , data=sub1)



summary(gam_model2)
gam.check(gam_model2)
AIC(gam_model2)


East_residuals = gam_model2$residuals
hist(East_residuals)

shapiro.test(East_residuals)

East_DF <- data.frame(East_residuals)

EastSTIDF <- STIDF(WindSP,WindTM,data=East_DF)

stplot(EastSTIDF)

no_cores <- detectCores()

# spatio-temporal variogram of North Wind. We can't really say there is a trend with the altitude
East_var <- variogramST(East_residuals~1,data=EastSTIDF,assumeRegular=T,na.omit=T, cores=no_cores, tunit='secs')


plot(East_var ,map=F)
plot(East_var ,map=T)
plot(East_var ,wireframe=T)


# To save time we can load the file next times
write.csv(East_var,"East_var_sec_MS2.csv", row.names = TRUE)
# East_var <- read.csv("East_var_sec_MS2.csv")
# East_var = East_var[-1]
# class(East_var) = c("StVariogram", "data.frame")


####
#SEPARABLE
separable2 <- vgmST("separable", space = vgm(0.08,"Sph", 28000, 0.01),time = vgm(0.08,"Sph", 500, 0.01), sill=0.08) 

separable_Vgm2 <- fit.StVariogram(East_var, separable2, fit.method=6, method="L-BFGS-B", tunit="secs")
attr(separable_Vgm2,"MSE")

plot(East_var, separable_Vgm2, wireframe=T, all=T)
plot(East_var, separable_Vgm2, map=F, all=T)
####

###
#PRODUCT SUM
prodSumModel2 <- vgmST("productSum",space = vgm(0.2, "Gau", 28000, 0.1),time = vgm(0.2, "Gau", 50, 1),k = 10) 
prodSumModel_Vgm2 <- fit.StVariogram(East_var, prodSumModel2,method = "L-BFGS-B", fit.method=6)
attr(prodSumModel_Vgm2, "MSE")

plot(East_var, prodSumModel_Vgm2, wireframe=T, all=T)
plot(East_var, prodSumModel_Vgm2, map=F, all=T)
###

####
# METRIC
metric2 <- vgmST("metric", joint = vgm(0.2,"Mat", 28000, 0.1), stAni=0.1) 
metric_Vgm2 <- fit.StVariogram(East_var, metric2, method="L-BFGS-B", fit.method = 6)
attr(metric_Vgm2, "MSE")

plot(East_var, metric_Vgm2, wireframe=T, all=T)
plot(East_var, metric_Vgm2, map=F, all=T)
####

####
# SUM-METRIC
sumMetric2 <- vgmST("sumMetric", space = vgm(psill=0.2,"Sph", range=28000, nugget=0),time = vgm(psill=0.03,"Sph", range=82, nugget=0), joint = vgm(psill=0,"Sph", range=10, nugget=0), stAni=0.1) 

sumMetric_Vgm2 <- fit.StVariogram(East_var, sumMetric2, method="L-BFGS-B",tunit="secs",  fit.method = 6)
attr(sumMetric_Vgm2, "MSE")

plot(East_var, sumMetric2, wireframe=T, all=T)
plot(East_var, sumMetric_Vgm2, wireframe=T, all=T)
plot(East_var, sumMetric_Vgm2, wireframe=T, diff=TRUE)
plot(East_var, sumMetric_Vgm2, map=F, all=T)

####
# SIMPLE SUM-METRIC
SimplesumMetric2 <- vgmST("simpleSumMetric",space = vgm(0.2,"Sph", 28000, 0),time = vgm(0.2,"Sph", 10, 0), joint = vgm(0.2,"Sph", 10, 0), nugget=2, stAni=0.1) 
SimplesumMetric_Vgm2 <- fit.StVariogram(East_var, SimplesumMetric2,method = "L-BFGS-B")
attr(SimplesumMetric_Vgm2, "MSE")

plot(East_var, SimplesumMetric_Vgm2, wireframe=T, all=T)
plot(East_var, SimplesumMetric_Vgm2, map=F, all=T)



plot(East_var,list(separable_Vgm2, prodSumModel_Vgm2, metric_Vgm2 ,sumMetric_Vgm2, SimplesumMetric_Vgm2), wireframe=T, all=T)


sumMetric2 <- vgmST("sumMetric", space = vgm(psill=0.2,"Sph", range=100000, nugget=0),time = vgm(psill=0.1,"Sph", range=10, nugget=0), joint = vgm(psill=0.1,"Sph", range=10, nugget=0), stAni=0.1) 

sumMetric_Vgm2 <- fit.StVariogram(East_var, sumMetric2, method="L-BFGS-B",tunit="secs")
attr(sumMetric_Vgm2, "MSE")

plot(East_var, sumMetric2, wireframe=T, all=T)
plot(East_var, sumMetric_Vgm2, wireframe=T, all=T)
plot(East_var, sumMetric_Vgm2, wireframe=T, diff=TRUE)
plot(East_var, sumMetric_Vgm2, map=F, all=T)


# Universal Kriging
attr(sumMetric_Vgm2, "temporal unit") <- units(abs(outer(index(EastSTIDF@time[1]), index(EastSTIDF@time[1]), "-")))
attr(separable_Vgm2, "temporal unit") <- units(abs(outer(index(EastSTIDF@time[1]), index(EastSTIDF@time[1]), "-")))

pred2 = krigeST(East_residuals~1, data=EastSTIDF, modelList=sumMetric_Vgm2, newdata=grid.sti, computeVar = T, nmax=10) 
plot(as.numeric(unlist(pred2@data[1])), datos$East_comp[time_interval])
stplot(pred2)

trend_pred2 = predict(gam_model2, datos[time_interval, c(3,4,5, 14)])
prediction2 = trend_pred2 + pred2@data[1]


MSE(prediction2$var1.pred,datos$East_comp[time_interval])
MSE(trend_pred2,datos$East_comp[time_interval])


plot(prediction2$var1.pred, datos$East_comp[time_interval], xlab="Predicted values of east velocity", ylab="Observed values of east velocity",main="Predicted vs Observed values of east velocity")
abline(0,1, col=2)
plot(trend_pred2, datos$East_comp[time_interval])
abline(0,1, col=2)
stplot(pred2)


cbind(data.frame(prediction2), data.frame(datos$East_comp[time_interval]))
cbind(trend_pred2, data.frame(datos$East_comp[time_interval]))


# CARTESIAN TO POLAR COORDINATES

speed = c()
direction = c()
for(i in 1:length(datos$East_comp[time_interval])){
  speed[i] = exp(cart2pol(prediction$var1.pred[i], prediction2$var1.pred[i], degrees = TRUE)$r)
  direction[i] = cart2pol(prediction$var1.pred[i], prediction2$var1.pred[i], degrees = TRUE)$theta
}

results = cbind(data.frame(speed), data.frame(direction)) 

true_values = datos[time_interval,c(11,12,16,17)]


plot(results$speed, true_values$WindSpeed)
abline(0,1,col=2)
plot(results$direction, true_values$WindDir, xlab="Predicted angles", ylab="Observed angles", main="Predicted vs Observed angles of the wind direction")
abline(0,1, col=2)

require(gridExtra)

plot1  = windrose(speed = log(results$speed),
                  direction = results$direction,
                  ggtheme='minimal') + labs(title="Predicted wind speed") + theme(legend.position = "none")


plot2 = windrose(speed = log(true_values$WindSpeed),
                 direction = true_values$WindDir,
                 ggtheme='minimal') + labs(title="Observed wind speed") + theme(legend.position = "none")

grid.arrange(plot1, plot2, ncol=2)

aux = which(true_values$WindDir<350 & true_values$WindDir > 50 )
plot(results$direction[aux], true_values$WindDir[aux], xlab="Predicted angles", ylab="Observed angles", main="Predicted vs Observed angles of the wind direction")
abline(0,1, col=2)

options(scipen=0)

str(factor(findInterval(results$direction,int.dir)))

# Only directions in the third and fourth quadrant are predicted
results$dir=factor(findInterval(results$direction,int.dir),labels = as.character(paste("deg",int.dir[c(-1,-2)],sep="")))

str(factor(findInterval(results$speed,int.ws)))
results$Windpeed=factor(findInterval(results$speed,int.ws),labels = c("f1","fmed","f3","max"))

aux = which(prediction$var1.pred>-6 & prediction$var1.pred < 4)
par(mfrow=c(1,2))
plot(prediction$var1.pred[aux], datos$North_comp[time_interval][aux], xlab="Predicted values of north velocity", ylab="Observed values of north velocity",main="Predicted vs Observed values of north velocity")
abline(0, 1, col=2)
plot(prediction2$var1.pred[aux], datos$East_comp[time_interval][aux], xlab="Predicted values of east velocity", ylab="Observed values of east velocity",main="Predicted vs Observed values of east velocity")
abline(0, 1, col=2)
par(mfrow=c(1,1))

library(caret)
confusionMatrix(results$dir, true_values$dir)$table

confusionMatrix(results$Windpeed, true_values$speed)

true_values$WindSpeed-results$speed

aux = which(results$speed<100 & results$speed>-100 &datos$WindSpeed[time_interval]<100)

aux = which(datos$WindSpeed[time_interval]<200)
plot(results$speed[aux], datos$WindSpeed[time_interval][aux])

plot(results$direction, true_values$WindDir)
results$direction - true_values$WindDir

#########################################################
# REPEATING THE PROCESS WITH ANOTHER TREND: RANDOM FOREST
#########################################################



# NORTH

train_interval = datos$Time>=as.POSIXct('2016-10-14 07:44:45 CET')&datos$Time<=as.POSIXct('2016-10-14 07:47:59 CET')
sub1 <- datos[train_interval,]

time_interval = datos$Time>=as.POSIXct('2016-10-14 07:48:00 CET')&datos$Time<=as.POSIXct('2016-10-14 07:54:44 CET')

coordinates(sub1)=~Lat+Long+Alt
projection(sub1)=CRS("+init=epsg:4326")

library(randomForest)
rf.train <- randomForest(North_comp~Lat + Long + Alt + dtm, data = sub1,
                         ntree=length(train_interval),cutoff=c(0.7,0.3),mtry=3,importance=TRUE, do.trace=F)

#Find residuals by subtracting predicted from acutal values
err <- rf.train$predicted - sub1$North_comp

mean(err)

hist(err)


summary(err)

#Make data frame holding residuals and fitted values
df <- data.frame(Residuals=err, Fitted.Values=rf.train$predicted)

#Sort data by fitted values
df2 <- df[order(df$Fitted.Values),]

#Create plot
plot(Residuals~Fitted.Values, data=df2)

#Add origin line at (0,0) with grey color #8
abline(0,0, col=8)

#Add the same smoothing line from lm regression with color red #2
lines(lowess(df2$Fitted.Values, df2$Residuals), col=2)

aux = which(err> -0.1  & err< 0.1)
mean(err[aux])
hist(err[aux])

rf.pred <- predict(rf.train, newdata=datos[time_interval, ])

library(MLmetrics)
MSE(rf.pred, datos$North_comp[time_interval])

plot(rf.pred, datos$North_comp[time_interval])
abline(0,1, col=2)

Wind.UTM <- spTransform(sub1,CRS("+init=epsg:3395")) 


WindSP <- SpatialPoints(Wind.UTM@coords,CRS("+init=epsg:3395"))

North_DF <- data.frame(north_err=err)

WindTM <- as.POSIXct(Wind.UTM$Time,tz="CET")

NorthSTIDF <- STIDF(WindSP,WindTM,data=North_DF)

stplot(NorthSTIDF)

no_cores <- detectCores()

# spatio-temporal variogram of North Wind. I NEED A BETTER TREND MODEL
Forest_North_var <- variogramST(north_err~1,data=NorthSTIDF,assumeRegular=T,na.omit=T, cores=no_cores, tunit='secs')
plot(Forest_North_var ,map=F)
plot(Forest_North_var ,wireframe=T)
plot(Forest_North_var ,map=T)


write.csv(Forest_North_var,"Random_forest_var.csv", row.names = TRUE)
# Forest_North_var <- read.csv("Random_forest_var.csv")
# Forest_North_var = Forest_North_var[-1]
# class(Forest_North_var) = c("StVariogram", "data.frame")

# SUM-METRIC
sumMetric_Forest <- vgmST("sumMetric", space = vgm(psill=0.02,"Sph", range=200000, nugget=0),time = vgm(psill=0.05,"Sph", range=1000, nugget=0), joint = vgm(psill=0.05,"Sph", range=20000, nugget=0), stAni=3) 
sumMetric_Vgm_Forest <- fit.StVariogram(Forest_North_var, sumMetric_Forest, method="L-BFGS-B",tunit="secs",  fit.method = 6)
attr(sumMetric_Vgm_Forest, "MSE")

plot(Forest_North_var, sumMetric_Forest, wireframe=T, all=T)
plot(Forest_North_var, sumMetric_Vgm_Forest, wireframe=T, all=T)
plot(Forest_North_var, sumMetric_Vgm_Forest, wireframe=T, diff=TRUE)
plot(Forest_North_var, sumMetric_Vgm_Forest, map=F, all=T)

attr(sumMetric_Forest, "temporal unit") <- units(abs(outer(index(NorthSTIDF@time[1]), index(NorthSTIDF@time[1]), "-")))
attr(sumMetric_Vgm_Forest, "temporal unit") <- units(abs(outer(index(NorthSTIDF@time[1]), index(NorthSTIDF@time[1]), "-")))

### Checking accuracy:

time.grid = unique(datos$Time[time_interval,])

grid = datos[time_interval, c(3,4,5)]

coordinates(grid)=~Lat+Long+Alt
projection(grid)=CRS("+init=epsg:3395")
bbox(grid)

grid.sti = STI(grid, datos$Time[time_interval,])

pred_forest = krigeST(north_err~1, data=NorthSTIDF, modelList=sumMetric_Vgm_Forest, newdata=grid.sti, nmax=5) 
plot(as.numeric(unlist(pred_forest@data[1])), datos$North_comp[time_interval])


rf.pred <- predict(rf.train, newdata=datos[time_interval, ])
prediction = rf.pred + pred_forest@data[1]

library(MLmetrics)
MSE(prediction$var1.pred,datos$North_comp[time_interval])
MSE(rf.pred,datos$North_comp[time_interval])

plot(prediction$var1.pred, datos$North_comp[time_interval],  xlab="Predicted values of north speed", ylab="Observed values of north speed",main="Predicted vs Observed values of north speed")
abline(0,1, col=2)
plot(rf.pred, datos$North_comp[time_interval], xlab="Predicted values of north velocity", ylab="Observed values of north velocity",main="Predicted vs Observed values of north velocity")
abline(0,1, col=2)
cbind(data.frame(prediction$var1.pred), data.frame(datos$North_comp[time_interval]))




# EAST

rf.train2 <- randomForest(East_comp~Lat + Long + Alt, data = sub1,
                          ntree=length(train_interval),cutoff=c(0.7,0.3),mtry=3,importance=TRUE, do.trace=F)

#Find residuals by subtracting predicted from acutal values
err2 <- rf.train2$predicted - sub1$East_comp

mean(err2)

hist(err2)

#Make data frame holding residuals and fitted values
df3 <- data.frame(Residuals=err2, Fitted.Values=rf.train2$predicted)

#Sort data by fitted values
df4 <- df3[order(df3$Fitted.Values),]

#Create plot
plot(Residuals~Fitted.Values, data=df4)

#Add origin line at (0,0) with grey color #8
abline(0,0, col=8)

#Add the same smoothing line from lm regression with color red #2
lines(lowess(df4$Fitted.Values, df4$Residuals), col=2)


rf.pred2 <- predict(rf.train2, newdata=datos[time_interval, ])

MSE(rf.pred2, datos$East_comp[time_interval])

plot(rf.pred2, datos$East_comp[time_interval])
abline(0,1, col=2)
cbind(data.frame(rf.pred2), data.frame(datos$East_comp[time_interval]))



shapiro.test(East_residuals)

East_DF <- data.frame(East_err = err2)

EastSTIDF <- STIDF(WindSP,WindTM,data=East_DF)

stplot(EastSTIDF)

no_cores <- detectCores()


Forest_East_var <- variogramST(East_err~1,data=EastSTIDF,assumeRegular=T,na.omit=T, cores=no_cores, tunit='secs')


plot(Forest_East_var ,map=F)
plot(Forest_East_var ,map=T)
plot(Forest_East_var ,wireframe=T)


# To save time we can load the file next times
write.csv(Forest_East_var,"Forest_East_var_sec_MS2.csv", row.names = TRUE)
# Forest_East_var <- read.csv("Forest_East_var_sec_MS2.csv")
# Forest_East_var = Forest_East_var[-1]
# class(Forest_East_var) = c("StVariogram", "data.frame")

sumMetric2_Forest <- vgmST("sumMetric", space = vgm(psill=0.005,"Sph", range=1000, nugget=0.001),time = vgm(psill=0.005,"Sph", range=10, nugget=0), joint = vgm(psill=0.005,"Sph", range=10, nugget=0), stAni=0.1) 

sumMetric_Vgm2_Forest <- fit.StVariogram(Forest_East_var, sumMetric2_Forest, method="L-BFGS-B",tunit="secs")
attr(sumMetric_Vgm2_Forest, "MSE")

plot(Forest_East_var, sumMetric2_Forest, wireframe=T, all=T)
plot(Forest_East_var, sumMetric_Vgm2_Forest, wireframe=T, all=T)
plot(Forest_East_var, sumMetric_Vgm2_Forest, wireframe=T, diff=TRUE)
plot(Forest_East_var, sumMetric_Vgm2_Forest, map=F, all=T)

attr(sumMetric_Vgm2_Forest, "temporal unit") <- units(abs(outer(index(EastSTIDF@time[1]), index(EastSTIDF@time[1]), "-")))

pred2_forest = krigeST(err2~1, data=EastSTIDF, modelList=sumMetric_Vgm2_Forest, newdata=grid.sti) 
hist(as.numeric(unlist(pred2_forest@data[1])))
stplot(pred2_forest)

rf.pred2 <- predict(rf.train2, newdata=datos[time_interval, ])
prediction2 = rf.pred2  + pred2_forest@data[1]

library(MLmetrics)
MSE(prediction2$var1.pred,datos$East_comp[time_interval])
MSE(rf.pred2,datos$East_comp[time_interval])



plot(prediction2$var1.pred, datos$East_comp[time_interval], xlab="Predicted values of east velocity", ylab="Observed values of east velocity",main="Predicted vs Observed values of east velocity")
abline(0,1, col=2)
plot(rf.pred2, datos$East_comp[time_interval], xlab="Predicted values of east velocity", ylab="Observed values of east velocity",main="Predicted vs Observed values of east velocity")
abline(0,1, col=2)

speed = c()
direction = c()
for(i in 1:length(datos$East_comp[time_interval])){
  speed[i] = exp(cart2pol(prediction$var1.pred[i], prediction2$var1.pred[i], degrees = TRUE)$r)
  direction[i] = cart2pol(prediction$var1.pred[i], prediction2$var1.pred[i], degrees = TRUE)$theta
}

results = cbind(data.frame(speed), data.frame(direction)) 

true_values = datos[time_interval,c(11,12,16,17)]


aux=which(true_values$WindSpeed < 200)
plot(results$speed[aux], true_values$WindSpeed[aux])
abline(0,1,col=2)
plot(results$direction, true_values$WindDir, xlab="Predicted angles", ylab="Observed angles", main="Predicted vs Observed angles of the wind direction")
abline(0,1, col=2)


require(gridExtra)

plot1  = windrose(speed = log(results$speed),
         direction = results$direction,
         ggtheme='minimal') + labs(title="Predicted wind speed") + theme(legend.position = "none")


plot2 = windrose(speed = log(true_values$WindSpeed),
         direction = true_values$WindDir,
         ggtheme='minimal') + labs(title="Observed wind speed")+ theme(legend.position = "none")

grid.arrange(plot1, plot2, ncol=2)




plot1  = windrose(speed = log(results$speed[aux]),
                  direction = results$direction[aux],
                  ggtheme='minimal') + labs(title="Predicted wind speed") + theme(legend.position = "none")


plot2 = windrose(speed = log(true_values$WindSpeed[aux]),
                 direction = true_values$WindDir[aux],
                 ggtheme='minimal') + labs(title="Observed wind speed")+ theme(legend.position = "none")

grid.arrange(plot1, plot2, ncol=2)

str(factor(findInterval(results$direction,int.dir)))
str(factor(findInterval(true_values$WindDir,int.dir)))

# Only directions in the third and fourth quadrant are predicted
results$dir=factor(findInterval(results$direction,int.dir),labels = as.character(paste("deg",int.dir[c(-1)],sep="")))

true_values$dir=factor(findInterval(true_values$WindDir,int.dir),labels =
                   as.character(paste("deg",int.dir[-1],sep="")))

confusionMatrix(results$dir, true_values$dir)$table

plot(Forest_North_var, Forest_East_var, all=T)

