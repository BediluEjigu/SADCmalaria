############################################################
#  Modeling Malaria  Data  using  Nonstationary MBG        #
#      Using Data from five SDAC countries                 #
############################################################
rm(list=ls())
# Loading libraries and our function to fit the non-stationary model
library(geoR)
library(PrevMap)
library(ncf)
library(sf)
library(dplyr)
library(raster)
source("NSgeo_MLE.R") # Calling our non-stationary R-module

# Load  clean dataset 
 # In this version of the data child level malaria status 
 # aggregated at cluster level and combined with cluster level environmental factors
ao=read.csv("aof.csv")
bf=read.csv("bff.csv")
mw=read.csv("mwf.csv")
mz=read.csv("mzf.csv")
tz=read.csv("tzf.csv")

# Loading Shapefile of SDAC countries  
library(sf)
ao_admn0=st_read("Shapefiles/AGO/ago_admbnda_adm2_gadm_ine_ocha_20180904.shp")
ao_admn1=st_read("Shapefiles/AGO/geoBoundaries-AGO-ADM1_simplified.shp")
bf_admn0=st_read("Shapefiles/BFA/bfa_admbnda_adm1_igb_20200323.shp")
mw_admn0=st_read("Shapefiles/MWI/mwi_admbnda_adm1_nso_hotosm_20230405.shp")
mz_admn0=st_read("Shapefiles/MOZ/moz_admbnda_adm1_ine_20190607.shp")
tz_admn0=st_read("Shapefiles/TZA/tza_admbnda_adm1_20181019.shp")

par(mfrow=c(2,3),mar=c(4,4,2,2))
plot(ao_admn0$geometry,main="AGO")
plot(bf_admn0$geometry,main="BFA")
plot(mw_admn0$geometry,main="MWI")
plot(mz_admn0$geometry,main="MOZ")
plot(tz_admn0$geometry,main="TZA")

# Map the prevalence of malaria at each cluster 
library(ggplot2)
library(viridis)
library(leaflet)
par(mfrow=c(2,3),mar=c(4,4,2,2))
plot(ao_admn1$geometry, main="ANG");points(ao$lon,ao$lat,cex=ao$Prevalence,ylab ="Latitude",xlab="Longitude",col="red")
plot(bf_admn0$geometry,main="BFS");points(bf$lon,bf$lat,cex=bf$Prevalence,ylab ="Latitude",xlab="Longitude",col="red")
plot(mw_admn0$geometry, main="MLW");points(mw$lon,mw$lat,cex=mw$Prevalence,ylab ="Latitude",xlab="Longitude",col="red")
plot(mz_admn0$geometry,main="MZQ");points(mz$lon,mz$lat,cex=mz$Prevalence,ylab ="Latitude",xlab="Longitude",col="red")
plot(tz_admn0$geometry,main="TZA");points(tz$lon,tz$lat,cex=tz$Prevalence,ylab ="Latitude",xlab="Longitude",col="red")



# Case 1: leaflet plot showing the prevalence at each reporting facility
pal <- colorBin("viridis", bins = c(0, 0.2,0.4,0.6,0.8,1.0))
leaflet(ao)%>%addProviderTiles(providers$CartoDB.Positron) %>%
  addCircles(lng = ~lon, lat = ~lat, color = ~ pal(Prevalence)) %>%
  addLegend("bottomleft", pal = pal, values = ~Prevalence,
            title = "U5 malaria prevance" ) %>% addScaleBar(position = c("bottomleft"))
# Arrange prediction location
library(sf)
ao.utm <- st_transform(ao_admn0,32736)
st_crs(ao.utm)
ao.union <- st_union(ao.utm)
ao.grid.sq <-st_make_grid(ao.utm,cellsize =10000,what="centers")
ao.inout <- st_intersects(ao.grid.sq, ao.union, sparse = FALSE)
ao.grid <- ao.grid.sq[ao.inout]
# Save the grid as a matrix
aod=matrix(data=NA,nrow=934,ncol=2)
for(i in 1:934)
{
  aod[i,]=rbind(ao.grid[[i]])
}
colnames(aod) <- c("utm_x","utm_y")
write.csv(aod,file="mw.pred10km.csv")

# Load Prediction locations
ao.pd=read.csv("ao.pred10km.csv")
bf.pd=read.csv("bf.pred10km.csv")
mw.pd=read.csv("mw.pred10km.csv")
mz.pd=read.csv("mz.pred10km.csv")
tz.pd=read.csv("tz.pred10km.csv")

# Extract geospatial covariates at each production locations
library(terra)
ANC <- read.csv("ANC2015.csv")
# raster population density in 1km grid
SA_pop=terra::rast("Shapefiles/PopDensity/pd_2019_1km.tif")
# Survey locations
sdp=as.data.frame(cbind(ANC$lon, ANC$lat))
# Get  the values of the raster at each point
valueatpoints=extract(SA_pop,sdp)
cbind(sdp,valueatpoints)

library("raster")
dataset <- getData (name ="worldclim", var = "tmax", res = 10)
install.packages("MODIStsp")
library(MODIStsp)
MODIStsp_get_prodlayers("M*D13Q1")

# Transform lat-log to UTM:AO:32S, BF:30N,
# AO=32634S, BF=32630N, WW=32736,MZ=32736,TZ=32736
library(sf)
bf_sf <-st_as_sf(bf,coords = c("lon", "lat"));st_crs(bf_sf)<-4326
# Converting lat lon to UTM
projUTM <- "+proj=utm +zone=32 +north +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"
d <-st_transform(bf_sf, crs = projUTM)
coo <- as.data.frame(st_coordinates(d))
bf$utm_x=coo$X;bf$utm_y=coo$Y;

##################################################
#           Exploring association                # 
##################################################
par(mfrow=c(3,2),mar=c(4,4,2,2))
plot(Logit ~ temp, data = bf,main="(a)",xlab="Average temprature",pch=20,cex=0.5)
plot(Logit ~ alt, data = bf,main="(b)",xlab="Altitude",pch=20,cex=0.5)
plot(Logit ~ prec, data=bf,main="(c)",xlab="Precipitation",pch=20,cex=0.5)
plot(Logit ~ evi, data=bf,main="(d)",xlab="Enhanced vegation index",pch=20,cex=0.5)
plot(Logit ~ itncoverage, data=bf,main="(e)",xlab="ITN coverage",pch=20,cex=0.5)

# To See how evi have an impact of the variability of malaria prevalence
 # 1st split the data by quartiles
# Split the data into different groups/quartailes to see how covariates affect the variability
ao$eviq<-cut(ao$evi,quantile(ao$evi),include.lowest=TRUE,labels=FALSE)
mw$eviq<-cut(mw$evi,quantile(mw$evi),include.lowest=TRUE,labels=FALSE)

par(mfrow=c(1,2),mar=c(4,4,2,2))
boxplot(Logit ~ eviq, data=ao,main="Angola",xlab="Quartiles of enhanced vegation index",ylab="logit of malaria prevalence",pch=20,cex=0.5)
boxplot(Logit ~ eviq, data=mw,main="Malawi",xlab="Quartiles of enhanced vegation index",ylab="logit of malaria prevalence",pch=20,cex=0.5)
ao$EVI=ao$evi
mw$EVI=mw$evi
ao.evi=as.geodata(ao,coords=c(1,16), data.col=7)
evi.vgao=variog(ao.evi);plot(evi.vgao,xlab="Vegetation index difference", main="")
mw.evi=as.geodata(mw,coords=c(1,16), data.col=7)
evi.vgmw=variog(mw.evi);plot(evi.vgmw,xlab="Vegetation index difference", main="")

library(ggplot2)
# Temperature
plot.temp <-ggplot(bf, aes(x=temp,y=Logit))+ geom_point() +
  labs(x="Temperature",y="Empirical logit")+
  stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x,
              col="green",lty="dashed",se=FALSE)
# EVI
plot.evi <- ggplot(bf, aes(x=evi, y=Logit)) + geom_point() +
  labs(x="Enhanced vegatation index",y="Empirical logit")+
  stat_smooth(method ="gam", formula = y ~ s(x),se=FALSE)+
  stat_smooth(method ="lm", formula = y ~ x,
              col="red",lty="dashed",se=FALSE)

# ITN Coverage 
plot.itn <- ggplot(bf, aes(x=itncoverage, y =Logit))+ geom_point() +
  labs(x="ITN Coverage",y="Empirical logit") +
  stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x ,
              col="red",lty="dashed",se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x ,
              col="green",lty="dashed",se=FALSE)
# Precipitation
plot.prec <- ggplot(bf, aes(x=prec,y = Logit)) + geom_point() +
  labs(x="Precipitation",y="Empirical logit")+
  stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x ,
              col="green",lty="dashed",se=FALSE)
# Altitude
plot.alt <- ggplot(bf, aes(x=alt, y=Logit)) + geom_point() +
  labs(x="Altitude",y="Empirical logit")+
  stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE)+
  stat_smooth(method ="lm", formula = y~x,
              col="red",lty="dashed",se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x,
              col="green",lty="dashed",se=FALSE)

library(ggcorrplot)
var.names <- c("Logit","evi","temp","itncoverage","alt","prec")
plot.cor <- ggcorrplot(corr=cor(bf[,var.names]), type ="lower",
                       ggtheme=ggplot2::theme_minimal,
                       hc.order = FALSE, show.diag = FALSE,
                       outline.col = "white",
                       lab=TRUE, legend.title="Correlation",
                       tl.cex =11, tl.srt = 55)
plot.cor

library(gridExtra)
grid.arrange(plot.alt,
             plot.itn,
             plot.prec,
             plot.evi,
             plot.temp,
             plot.cor)

# Split the data into different groups/quartailes to see how covariates affect the variability
tz$altq<-cut(tz$alt,quantile(tz$alt),include.lowest=TRUE,labels=FALSE)
tz$tempq<-cut(tz$temp,quantile(tz$temp),include.lowest=TRUE,labels=FALSE)
tz$precq<-cut(tz$prec,quantile(tz$prec),include.lowest=TRUE,labels=FALSE)
tz$eviq<-cut(tz$evi,quantile(tz$evi),include.lowest=TRUE,labels=FALSE)
tz$itnq<-cut(tz$itncoverage,quantile(tz$itncoverage),include.lowest=TRUE,labels=FALSE)

par(mfrow=c(2,2),mar=c(4,4,2,2))
boxplot(Logit ~altq, data=tz, xlab ="Altitude ",
        ylab = "logit of malaria prevalence", main = "(a)")
boxplot(Logit ~tempq, data=tz, xlab ="Average temprature",
        ylab = "logit of malaria prevalence", main = "(b)")
boxplot(Logit ~eviq, data=tz, xlab = "Enhanced vegation index",
        ylab = "logit of malaria prevalence", main = "(c)")
boxplot(Logit ~itnq, data=tz, xlab = "ITN Coverage",
        ylab = "logit of malaria prevalence", main = "(d)")

# Testing for residual spatial correlation using the variogram
library(PrevMap)
ao.vaiog=spat.corr.diagnostic(pos ~ prec+temp + evi+ itncoverage + alt, units.m = ~ tested,
                              coords = ~utm_x+utm_y,likelihood = "Binomial",which.test="variogram", data=ao)
bf.vaiog=spat.corr.diagnostic(pos ~prec+temp + evi+ itncoverage + alt, units.m = ~ tested,
                              coords = ~utm_x+utm_y,likelihood = "Binomial",which.test="variogram", data=bf)
mw.vaiog=spat.corr.diagnostic(pos ~prec+temp + evi+ itncoverage +alt,units.m = ~ tested,
                              coords = ~utm_x+utm_y,likelihood = "Binomial",which.test="variogram", data=mw)
mz.vaiog=spat.corr.diagnostic(pos ~prec+temp + evi+ itncoverage +alt,units.m = ~ tested,
                              coords = ~utm_x+utm_y,likelihood = "Binomial", which.test="variogram", data=mz)
tz.vaiog=spat.corr.diagnostic(pos ~ prec+temp + evi+itncoverage + alt, units.m = ~ tested,
                              coords = ~utm_x+utm_y,likelihood = "Binomial", which.test="variogram",data=tz)
# Plot Variograms demonstarting the need to consider spatial modeling
par(mfrow=c(3,2),mar=c(4,4,2,2))
plot(ao.vaiog,xlab="Distance (km)",ylab="Variogram", main="(a)")
plot(bf.vaiog,xlab="Distance (km)",ylab="Variogram", main="(b)")
plot(mw.vaiog,xlab="Distance (km)",ylab="Variogram", main="(c)")
plot(mz.vaiog,xlab="Distance (km)",ylab="Variogram", main="(d)")
plot(tz.vaiog,xlab="Distance (km)",ylab="Variogram", main="(e)")

############################################
# Variogram and Correlogram Exploration    #
############################################
# a) variogam
library(geoR)
par(mfrow=c(2,2),mar=c(4,4,2,2))
# Variogam based on geographical distance and covariates
tz$EVI=tz$evi
tz$Alt=tz$alt
tz$ITN=tz$itncoverage
tz$Temp=tz$temp
names(tz)
geodata.vg=as.geodata(tz,coords=12:13, data.col=7)
dist.vg=variog(geodata.vg);plot(dist.vg, xlab="Euclidean distance(km)",main="(a)")

geodata.evi=as.geodata(tz,coords=c(1,14), data.col=7)
evi.vg=variog(geodata.evi);plot(evi.vg,xlab="Vegetation index difference", main="(b)")

geodata.alt=as.geodata(tz,coords=c(11,15), data.col=7)
alt.vg=variog(geodata.alt);plot(alt.vg,xlab="Altitude difference", main="(c)")

geodata.itn=as.geodata(tz,coords=c(2,16), data.col=7)
itn.vg=variog(geodata.itn);plot(itn.vg, xlab="Difference in ITN coverage",main="(d)")

##############################################
#      Geostatistical Modeling               #
##############################################

#Angola: temp+prec+alt+itncoverage+evi

# Standard MBG 
ID.coords0 <-create.ID.coords(ao,~utm_x + utm_y)
MBGao <- linear.model.MLE(formula=Logit~temp +evi_centered, coords = ~utm_x + utm_y, 
                          data =ao, start.cov.pars = c(2, 1),ID.coords=ID.coords0,fixed.rel.nugget=0.5, kappa=0.5)
summary(MBGao, log.cov.pars = FALSE)
# 95% CI for parameter estimates
beta.MBGao <- summary(MBGao)$coefficients[,1:2]
cbind(beta.MBGao, beta.MBGao[,1]-qnorm(0.975)*beta.MBGao[,2],
      beta.MBGao[,1]+qnorm(0.975)*beta.MBGao[,2])
# Model validation for linear geostatistical model using Monte Carlo methods based on the variogram
variog.diag.ao <-variog.diagnostic.lm(MBGao, n.sim =1000,uvec = NULL,plot.results = TRUE,
                                      range.fact = 1, which.test = "variogram",param.uncertainty = FALSE)
# Non-stationary MBG (NS-MBG)
ID.coords_evi <- create.ID.coords(ao,~utm_x + utm_x+evi)
NS_MBGao<- NSgeo.MLE(Logit ~temp+prec+alt+itncoverage+evi,coords = ~utm_x + utm_y,
                     rand.effect.domain=~evi, start.cov.pars=c(2,1),fixed.rel.nugget=0.5,
                     data=ao,ID.coords=ID.coords_evi,method="nlminb",
                     messages = TRUE,return_se=TRUE)
#Parameter estimates with 95%CI
beta.ao=NS_MBGao$par; se.beta.ao=NS_MBGao$se
ll.ao=beta.ao-1.96*se.beta.ao; ul.ao=beta.ao+1.96*se.beta.ao
CI.ao=cbind(beta.ao,se.beta.ao,ll.ao,ul.ao);round(CI.ao,3)

# Brukina Faso
ID.coords0 <- create.ID.coords(bf,~utm_x + utm_y)
MBGbf <- linear.model.MLE(formula=Logit~temp+prec+alt+itncoverage+evi_centered, coords = ~utm_x + utm_y, 
                          data =bf, start.cov.pars = c(2, 1),ID.coords=ID.coords0,fixed.rel.nugget=0.5, kappa=0.5, method="nlminb")
summary(MBGbf, log.cov.pars = FALSE)
# 95% CI for parameter estimates
beta.MBGbf <- summary(MBGbf)$coefficients[,1:2]
bfbeta=cbind(beta.MBGbf, beta.MBGbf[,1]-qnorm(0.975)*beta.MBGbf[,2],
             beta.MBGbf[,1]+qnorm(0.975)*beta.MBGbf[,2]);round(bfbeta,3)
# Model validation for linear geostatistical model using Monte Carlo methods based on the variogram
variog.diag.bf <-variog.diagnostic.lm (MBGbf, n.sim =1000,uvec = NULL,plot.results = TRUE,range.fact = 1, which.test = "variogram",param.uncertainty = FALSE)

# Non-stationary MBG (NS-MBG)
ID.coords_evi <- create.ID.coords(bf,~utm_x + utm_y+evi_centered)

NS_MBGbf<- NSgeo.MLE(Logit ~temp+prec+alt+itncoverage+ evi_centered,coords = ~utm_x + utm_y,
                     rand.effect.domain=~evi_centered, start.cov.pars=c(2,1),fixed.rel.nugget=0.5,
                     data=bf,ID.coords=ID.coords_evi,method="nlminb",
                     messages = TRUE,return_se=TRUE)
#Parameter estimates with 95%CI
beta.bf=NS_MBGbf$par; se.beta.bf=NS_MBGbf$se
ll.bf=beta.bf-1.96*se.beta.bf; ul.bf=beta.bf+1.96*se.beta.bf
CI.bf=cbind(beta.bf,se.beta.bf,ll.bf,ul.bf);round(CI.bf,3)

# Malawi
ID.coords0 <- create.ID.coords(mw,~utm_x + utm_y)
MBGmw <- linear.model.MLE(formula=Logit~temp+prec+alt+itncoverage+evi_centered, coords = ~utm_x + utm_y, 
                          data =mw, start.cov.pars = c(2, 1),ID.coords=ID.coords0,fixed.rel.nugget=0.5, kappa=0.5, method="nlminb")
summary(MBGmw, log.cov.pars = FALSE)
# 95% CI for parameter estimates
beta.MBGmw <- summary(MBGmw)$coefficients[,1:2]
mwbeta=cbind(beta.MBGmw, beta.MBGmw[,1]-qnorm(0.975)*beta.MBGmw[,2],
             beta.MBGmw[,1]+qnorm(0.975)*beta.MBGmw[,2]);round(mwbeta,3)
# Model validation for linear geostatistical model using Monte Carlo methods based on the variogram
variog.diag.mw <-variog.diagnostic.lm (MBGmw, n.sim =1000,uvec = NULL,plot.results = TRUE,range.fact = 1, which.test = "variogram",param.uncertainty = FALSE)

# Non-stationary MBG (NS-MBG)
ID.coords_evi <- create.ID.coords(mw,~utm_x + utm_y+evi_centered)

NS_MBGmw1<- NSgeo.MLE(Logit ~temp+prec+alt+itncoverage,coords = ~utm_x + utm_y,
                      rand.effect.domain=~evi_centered, start.cov.pars=c(2,1),fixed.rel.nugget=0.5,
                      data=mw,ID.coords=ID.coords_evi,method="nlminb",
                      messages = TRUE,return_se=TRUE)
#Case 2
NS_MBGmw<- NSgeo.MLE(Logit ~temp+prec+alt+itncoverage+evi_centered,coords = ~utm_x + utm_y,
                     rand.effect.domain=~evi_centered, start.cov.pars=c(2,1),fixed.rel.nugget=0.5,
                     data=mw,ID.coords=ID.coords_evi,method="nlminb",
                     messages = TRUE,return_se=TRUE)
# Parameter estimates with 95%CI
beta.mw=NS_MBGmw$par; se.beta.mw=NS_MBGmw$se
ll.mw=beta.mw-1.96*se.beta.mw; ul.mw=beta.mw+1.96*se.beta.mw
CI.mw=cbind(beta.mw,se.beta.mw,ll.mw,ul.mw);round(CI.mw,3)

# Mozambique
ID.coords0 <- create.ID.coords(mz,~utm_x + utm_y)
MBGmz <- linear.model.MLE(formula=Logit~temp+prec+alt+itncoverage+evi, coords = ~utm_x + utm_y, 
                          data =mz, start.cov.pars = c(2, 1),ID.coords=ID.coords0,fixed.rel.nugget=0.5, kappa=1.5, method="nlminb")
summary(MBGmz, log.cov.pars = FALSE)
# 95% CI for parameter estimates
beta.MBGmz <- summary(MBGmz)$coefficients[,1:2]
mz=cbind(beta.MBGmz, beta.MBGmz[,1]-qnorm(0.975)*beta.MBGmz[,2],
         beta.MBGmz[,1]+qnorm(0.975)*beta.MBGmz[,2]);round(mz,3)
# Model validation for linear geostatistical model using Monte Carlo methods based on the variogram
variog.diag.mz <-variog.diagnostic.lm (MBGmz, n.sim =1000,uvec = NULL,plot.results = TRUE,range.fact = 1, which.test = "variogram",param.uncertainty = FALSE)
# Non-stationary MBG (NS-MBG)
ID.coords_evi <- create.ID.coords(mz,~utm_x + utm_y+evi)
NS_MBGmz<- NSgeo.MLE(Logit ~temp+prec+alt+itncoverage+evi,coords = ~utm_x + utm_y,
                     rand.effect.domain=~evi, start.cov.pars=c(2,1),fixed.rel.nugget=0.5,
                     data=mz,ID.coords=ID.coords_evi,method="nlminb",
                     messages = TRUE,return_se=TRUE)
#Parameter estimates with 95%CI
beta.mz=NS_MBGmz$par; se.beta.mz=NS_MBGmz$se
ll.mz=beta.mz-1.96*se.beta.mz; ul.mz=beta.mz+1.96*se.beta.mz
CI.mz=cbind(beta.mz,se.beta.mz,ll.mz,ul.mz);round(CI.mz,3)

#Tanzania (precipitation, temp )
ID.coords0 <- create.ID.coords(tz,~utm_x + utm_y)
MBGtz <- linear.model.MLE(formula=Logit~temp+prec+alt+itncoverage+evi_centered, coords = ~utm_x + utm_y, 
                          data =tz, start.cov.pars = c(2, 1),ID.coords=ID.coords0,fixed.rel.nugget=0.5, kappa=1.5, method="nlminb")
summary(MBGtz, log.cov.pars = FALSE)
# 95% CI for parameter estimates
beta.MBGtz <- summary(MBGtz)$coefficients[,1:2]
tzbeta=cbind(beta.MBGtz, beta.MBGtz[,1]-qnorm(0.975)*beta.MBGtz[,2],
             beta.MBGtz[,1]+qnorm(0.975)*beta.MBGtz[,2]);round(tzbeta,3)
# Model validation for linear geostatistical model using Monte Carlo methods based on the variogram
variog.diag.tz <-variog.diagnostic.lm (MBGtz, n.sim =1000,uvec = NULL,plot.results = TRUE,range.fact = 1, which.test = "variogram",param.uncertainty = FALSE)
# Non-stationary MBG (NS-MBG)
ID.coords_evi <- create.ID.coords(tz,~utm_x + utm_y+evi)
NS_MBGtz<- NSgeo.MLE(Logit ~temp+itncoverage+evi,coords = ~utm_x + utm_y,
                     rand.effect.domain=~evi, start.cov.pars=c(2,1),fixed.rel.nugget=0.5,
                     data=tz,ID.coords=ID.coords_evi,method="nlminb",
                     messages = TRUE,return_se=TRUE)
#Parameter estimates with 95%CI
beta.tz=NS_MBGtz$par; se.beta.tz=NS_MBGtz$se
ll.tz=beta.tz-1.96*se.beta.tz; ul.tz=beta.tz+1.96*se.beta.tz
CI.tz=cbind(beta.tz,se.beta.tz,ll.tz,ul.tz);round(CI.tz,3)

# Spatial error diagnosis
par(mfrow=c(3,2),mar=c(4,4,2,2))
plot(variog.diag.ao,xlab="Distance (km)",ylab="Variogram", main="(a)")
plot(variog.diag.bf,xlab="Distance (km)",ylab="Variogram", main="(b)")
plot(variog.diag.mw,xlab="Distance (km)",ylab="Variogram", main="(c)")
plot(variog.diag.mz,xlab="Distance (km)",ylab="Variogram", main="(d)")
plot(variog.diag.tz,xlab="Distance (km)",ylab="Variogram", main="(e)")

# Model Comparison using AIC
cat(paste ("Stationary MBGao AIC: ",  round(-2*MBGao$log.lik + 2*length(MBGao$par),3), "\n"))
cat(paste ("Nonstationary Model NS_MBGao : ", round(2*length(NS_MBGao$par)+2*NS_MBGao$objective,3), "\n"))

cat(paste ("Stationary MBGbf AIC: ",  round(-2*MBGbf$log.lik + 2*length(MBGbf$par),3), "\n"))
cat(paste ("Nonstationary Model NS_MBGbf : ", round(2*length(NS_MBGbf$par)+2*NS_MBGbf$objective,3), "\n"))

cat(paste ("Stationary MBGmw AIC: ",  round(-2*MBGmw$log.lik + 2*length(MBGmw$par),3), "\n"))
cat(paste ("Nonstationary Model NS_MBGmw : ", round(2*length(NS_MBGmw$par)+2*NS_MBGmw$objective,3), "\n"))

cat(paste ("Stationary MBGmz AIC: ",  round(-2*MBGmz$log.lik + 2*length(MBGmz$par),3), "\n"))
cat(paste ("Nonstationary Model NS_MBGmz : ", round(2*length(NS_MBGmz$par)+2*NS_MBGmz$objective,3), "\n"))

cat(paste ("Stationary MBGtz AIC: ",  round(-2*MBGtz$log.lik + 2*length(MBGtz$par),3), "\n"))
cat(paste ("Nonstationary Model NS_MBGtz : ", round(2*length(NS_MBGtz$par)+2*NS_MBGtz$objective,3), "\n"))


#######################################################
#      Malaria Prediction  to unsampled location      #
#######################################################



