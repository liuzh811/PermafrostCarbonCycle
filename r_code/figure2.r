# 8/17/2022
# by Zhihua Liu, liuzh811@126.com

# fig 2, calculate the correlation between GPP and net carbon uptake

rm(list = ls())

library(rgdal)
library(raster)
library(rasterVis)
require(ncdf4)
library(ggplot2)
library(maps)

rasterOptions(tmpdir='R:/zhihua/modis_data/tmp2')

#NSIDC EASE-Grid 2.0 GloNDVIl: https://epsg.io/6933
proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
proj.easegrid2 = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 
proj.easegrid2_N = "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
proj.easegrid2_S = "+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "
 
proj_polar = "+proj=stere +lat_0=90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
 
#### read into spatial files
region.grd1d = raster("R:/zhihua/L4C_NH_analysis_Revision/spatialfires/region.grd1d.tif") 

pctbor.grd1d = raster("R:/zhihua/L4C_NH_analysis_Revision/spatialfires/pctbor.grd1d.tif")
pctbor.grd1d = pctbor.grd1d > 80
pctbor.grd1d[pctbor.grd1d == 0] = NA

frost.grd1d = raster("R:/zhihua/L4C_NH_analysis_Revision/spatialfires/frost.grd1d.tif")
frost.grd1d = frost.grd1d + 1

biome.grd1d = raster("R:/zhihua/L4C_NH_analysis_Revision/biome.grd1d.tif")

lc.grd1d = raster("R:/zhihua/L4C_NH_analysis_Revision/spatialfires/lc.grd1d.tif")

mask.1d = raster("R:/zhihua/L4C_NH_analysis_Revision/mask.1d.tif")
mask.1d = mask.1d*lc.grd1d

region.shp = readOGR( "R:/zhihua/L4C_NH_analysis", "region")
frost.shp = readOGR("R:/zhihua/L4C_NH_analysis_Revision/spatialfires", "frost.shp")

biome.grd1d = biome.grd1d*mask.1d
plot(biome.grd1d)
plot(region.shp, add = TRUE)
plot(frost.shp, border = "red", add = TRUE)

####################################################################################
tc.grd = raster("R:/zhihua/L4C_NH_analysis_Revision3/treecover_mean_1983-2016.tif")
perma.grd = raster("R:/zhihua/L4C_NH_analysis_Revision3/permaf_ext_mean_1997_2019_geo2.tif")
sv.grd = raster("R:/zhihua/L4C_NH_analysis_Revision3/shortvegetation_mean_1983-2016.tif")
mod_tc.grd = raster("R:/zhihua/L4C_NH_analysis_Revision3/modis_tc_cmg.tif")

tc.grd1d = aggregate(tc.grd, fact=20, fun=mean, na.rm = TRUE)
sv.grd1d = aggregate(sv.grd, fact=20, fun=mean, na.rm = TRUE)
perma.grd1d = aggregate(perma.grd, fact=20, fun=mean, na.rm = TRUE)
mod_tc.grd1d = aggregate(mod_tc.grd, fact=20, fun=mean, na.rm = TRUE)

#################################################################
# extract data and plot
#define a function to get all points for a raster
Ex.pts.all = function(x){
  proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
  pts = xyFromCell(x, seq(1, (nrow(x)*ncol(x))))
  pts.sp = SpatialPoints(coords = pts, proj4string = CRS(proj.geo))
  return(pts.sp)
}

#define a function to get all points for a raster
Ex.pts.all.nonNA = function(x){
  proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
  #pts = xyFromCell(x, seq(1, (nrow(x)*ncol(x))))
  pts = rasterToPoints(x) #get raster coordinate, from left to right, top to bottom
  pts = data.frame(pts)
  pts <- SpatialPoints(coords = cbind(pts$x,pts$y),proj4string = CRS(proj.geo))
  pts.sp = SpatialPoints(coords = pts, proj4string = CRS(proj.geo))
  return(pts.sp)
}

pts.sp = Ex.pts.all.nonNA(mask.1d)

Ex.pts.all.nonNA2 = function(x){
  proj.geo = projection(x)
  #pts = xyFromCell(x, seq(1, (nrow(x)*ncol(x))))
  pts = rasterToPoints(x) #get raster coordinate, from left to right, top to bottom
  pts = data.frame(pts)
  pts <- SpatialPoints(coords = cbind(pts$x,pts$y),proj4string = CRS(proj.geo))
  pts.sp = SpatialPoints(coords = pts, proj4string = CRS(proj.geo))
  return(pts.sp)
}


# define a function to change points to raster
Point2raster = function(points, raster){ #points is vector, raster the template
  nr = nrow(raster)
  nc = ncol(raster)
  ext = extent(raster)
  proj.geo = projection(raster)
  r = raster(matrix(points, nrow = nr, ncol = nc, byrow = TRUE),
             xmn=ext@xmin, xmx=ext@xmax,
             ymn=ext@ymin, ymx=ext@ymax, 
             crs=CRS(proj.geo))
  return(r)
}

#find correlations
COR.test <- function(x) {
  x = as.numeric(x)
  if (length(which(!is.na(x[1:(length(x)/2)]))) < 10 | length(which(!is.na(x[(1+length(x)/2):length(x)]))) < 10) 
  {
    c(NA,NA)
  } else 
  {
    test = cor.test(as.numeric(x[1:(length(x)/2)]), as.numeric(x[(1+length(x)/2):length(x)]),na.action = na.omit)
    c(test$estimate,test$p.value)
  }
}

# trend function
library(spatialEco)
# raster.kendall(x, tau = FALSE, intercept = FALSE, p.value = FALSE, z.value = FALSE, confidence = FALSE, autocorrelation = FALSE, ...)

#length of NAs
lenNA = function(x){return(length(which(is.na(x))))}

# Write a function to calcuate the trend
Trend1 = function(stkgrd){
require(spatialEco)
stkgrd_lenNA=calc(stkgrd,lenNA)
stkgrd[stkgrd_lenNA > 10] = -9999 #replace NAs with -9999 to use the raster.kendall function
stkgrd_trend= raster.kendall(stkgrd,p.value = TRUE, na.rm = TRUE)
stkgrd_trend[stkgrd_lenNA > 10] = NA
return(stkgrd_trend)
}

########GMISS NDVI

NDVIJan.grd = stack("R:/zhihua/Temperature_NDVI_sensitivity/NDVIJan.tif")
NDVIFeb.grd = stack("R:/zhihua/Temperature_NDVI_sensitivity/NDVIFeb.tif")
NDVIMar.grd = stack("R:/zhihua/Temperature_NDVI_sensitivity/NDVIMar.tif")
NDVIApr.grd = stack("R:/zhihua/Temperature_NDVI_sensitivity/NDVIApr.tif")
NDVIMay.grd = stack("R:/zhihua/Temperature_NDVI_sensitivity/NDVIMay.tif")
NDVIJun.grd = stack("R:/zhihua/Temperature_NDVI_sensitivity/NDVIJun.tif")
NDVIJul.grd = stack("R:/zhihua/Temperature_NDVI_sensitivity/NDVIJul.tif")
NDVIAug.grd = stack("R:/zhihua/Temperature_NDVI_sensitivity/NDVIAug.tif")
NDVISep.grd = stack("R:/zhihua/Temperature_NDVI_sensitivity/NDVISep.tif")
NDVIOct.grd = stack("R:/zhihua/Temperature_NDVI_sensitivity/NDVIOct.tif")
NDVINov.grd = stack("R:/zhihua/Temperature_NDVI_sensitivity/NDVINov.tif")
NDVIDec.grd = stack("R:/zhihua/Temperature_NDVI_sensitivity/NDVIDec.tif")
NDVIAnnual.grd = stack("R:/zhihua/Temperature_NDVI_sensitivity/NDVIAnnual.tif")

NDVIJan.grd = crop(NDVIJan.grd, c(-180,180,49,75))
NDVIFeb.grd = crop(NDVIFeb.grd, c(-180,180,49,75))
NDVIMar.grd = crop(NDVIMar.grd, c(-180,180,49,75))
NDVIApr.grd = crop(NDVIApr.grd, c(-180,180,49,75))
NDVIMay.grd = crop(NDVIMay.grd, c(-180,180,49,75))
NDVIJun.grd = crop(NDVIJun.grd, c(-180,180,49,75))
NDVIJul.grd = crop(NDVIJul.grd, c(-180,180,49,75))
NDVIAug.grd = crop(NDVIAug.grd, c(-180,180,49,75))
NDVISep.grd = crop(NDVISep.grd, c(-180,180,49,75))
NDVIOct.grd = crop(NDVIOct.grd, c(-180,180,49,75))
NDVINov.grd = crop(NDVINov.grd, c(-180,180,49,75))
NDVIDec.grd = crop(NDVIDec.grd, c(-180,180,49,75))
NDVIAnnual.grd = crop(NDVIAnnual.grd, c(-180,180,49,75))

NDVIJan.grd[NDVIJan.grd < 0.1] = NA
NDVIFeb.grd[NDVIFeb.grd < 0.1] = NA
NDVIMar.grd[NDVIMar.grd < 0.1] = NA
NDVIApr.grd[NDVIApr.grd < 0.1] = NA
NDVIMay.grd[NDVIMay.grd < 0.1] = NA
NDVIJun.grd[NDVIJun.grd < 0.1] = NA
NDVIJul.grd[NDVIJul.grd < 0.1] = NA
NDVIAug.grd[NDVIAug.grd < 0.1] = NA
NDVISep.grd[NDVISep.grd < 0.1] = NA
NDVIOct.grd[NDVIOct.grd < 0.1] = NA
NDVINov.grd[NDVINov.grd < 0.1] = NA
NDVIDec.grd[NDVIDec.grd < 0.1] = NA
NDVIAnnual.grd[NDVIAnnual.grd < 0.1] = NA
 
names(NDVIJan.grd) <- paste0("Y", 1982:2015)
names(NDVIFeb.grd) <- paste0("Y", 1982:2015)
names(NDVIMar.grd) <- paste0("Y", 1982:2015)
names(NDVIApr.grd) <- paste0("Y", 1982:2015)
names(NDVIMay.grd) <- paste0("Y", 1982:2015)
names(NDVIJun.grd) <- paste0("Y", 1982:2015)
names(NDVIJul.grd) <- paste0("Y", 1982:2015)
names(NDVIAug.grd) <- paste0("Y", 1982:2015)
names(NDVISep.grd) <- paste0("Y", 1982:2015)
names(NDVIOct.grd) <- paste0("Y", 1982:2015)
names(NDVINov.grd) <- paste0("Y", 1982:2015)
names(NDVIDec.grd) <- paste0("Y", 1982:2015)
names(NDVIAnnual.grd) <- paste0("Y", 1982:2015)
 

# only do early, late
SpringNDVI.grd = stack()
FallNDVI.grd = stack()
WinterNDVI.grd = stack()
AnnualNDVI.grd = stack()
GrowingNDVI.grd = stack()

for (i in 1:nlayers(NDVIMay.grd)){
t1 = calc(stack(NDVIMay.grd[[i]], NDVIJun.grd[[i]], NDVIJul.grd[[i]], NDVIAug.grd[[i]]),
           mean, na.rm = TRUE)

t2 = calc(stack(NDVISep.grd[[i]],NDVIOct.grd[[i]]), 
           mean, na.rm = TRUE)

t3 = calc(stack(NDVIJan.grd[[i]], NDVIFeb.grd[[i]], NDVIMar.grd[[i]], NDVIApr.grd[[i]], NDVINov.grd[[i]],NDVIDec.grd[[i]]),
           mean, na.rm = TRUE)
SpringNDVI.grd = addLayer(SpringNDVI.grd, t1)
FallNDVI.grd = addLayer(FallNDVI.grd, t2)
WinterNDVI.grd = addLayer(WinterNDVI.grd, t3)

t4 = calc(stack(NDVIJan.grd[[i]], NDVIFeb.grd[[i]], NDVIMar.grd[[i]], NDVIApr.grd[[i]], NDVIMay.grd[[i]], NDVIJun.grd[[i]], 
                NDVIJul.grd[[i]], NDVIAug.grd[[i]], NDVISep.grd[[i]], NDVIOct.grd[[i]], NDVINov.grd[[i]], NDVIDec.grd[[i]]),
           mean, na.rm = TRUE)

AnnualNDVI.grd = addLayer(AnnualNDVI.grd, t4)

t5 = calc(stack(NDVIMay.grd[[i]], NDVIJun.grd[[i]], 
                NDVIJul.grd[[i]], NDVIAug.grd[[i]], NDVISep.grd[[i]]),
           mean, na.rm = TRUE)

GrowingNDVI.grd = addLayer(GrowingNDVI.grd, t5)

}

SpringNDVI.grd = aggregate(SpringNDVI.grd, fact = 4, fun = "mean")
FallNDVI.grd = aggregate(FallNDVI.grd, fact = 4, fun = "mean")
WinterNDVI.grd = aggregate(WinterNDVI.grd, fact = 4, fun = "mean")
AnnualNDVI.grd = aggregate(AnnualNDVI.grd, fact = 4, fun = "mean")
NDVIAnnual.grd = aggregate(NDVIAnnual.grd, fact = 4, fun = "mean")

GrowingNDVI.grd = aggregate(GrowingNDVI.grd, fact = 4, fun = "mean")

############## NEE 
land.grd = stack("R:/zhihua/InversionCO2_GCB2020/NEE.monthly.mean1976.2017.grd")
land.grd = crop(land.grd, extent(-180,180, 50,75))

idx = which(as.numeric(as.character(substr(names(land.grd), 2,5))) >= 1980 & as.numeric(as.character(substr(names(land.grd), 2,5))) <= 2017) 
land = land.grd[[idx]]

land = land*mask.1d

M = c(paste0("M0", 1:9), paste0("M", 10:12) )
names(land) <- paste0("Y",rep(c(1980:2017), each = 12), rep(M, length(c(1980:2017))))

plot(land[[1]])

SpringC.grd = stack()
# SummerC.grd = stack()
FallC.grd = stack()
WinterC.grd = stack()
AnnualC.grd = stack()


for(yr in 1982:2015){
  
  idx = which((as.numeric(as.character(substr(names(land), 7,8))) >= 5 & as.numeric(as.character(substr(names(land), 7,8))) <= 8) &
             as.numeric(as.character(substr(names(land), 2,5))) == yr) 
  SpringC.grd = addLayer(SpringC.grd, calc(land[[idx]], mean, na.rm = TRUE)) 

  # idx = which((as.numeric(as.character(substr(names(land), 7,8))) >= 7 & as.numeric(as.character(substr(names(land), 7,8))) <= 8) &
             # as.numeric(as.character(substr(names(land), 2,5))) == yr) 
  # SummerC.grd = addLayer(SummerC.grd, calc(land[[idx]], mean, na.rm = TRUE)) 
  
   idx = which((as.numeric(as.character(substr(names(land), 7,8))) >= 9 & as.numeric(as.character(substr(names(land), 7,8))) <= 10) &
             as.numeric(as.character(substr(names(land), 2,5))) == yr) 
  FallC.grd = addLayer(FallC.grd, calc(land[[idx]], mean, na.rm = TRUE)) 
  
   idx = which(((as.numeric(as.character(substr(names(land), 7,8))) >= 1 & as.numeric(as.character(substr(names(land), 7,8))) <= 4) |  
                (as.numeric(as.character(substr(names(land), 7,8))) >= 11 & as.numeric(as.character(substr(names(land), 7,8))) <= 12))&
             as.numeric(as.character(substr(names(land), 2,5))) == yr) 
  WinterC.grd = addLayer(WinterC.grd, calc(land[[idx]], mean, na.rm = TRUE)) 
  
     idx = which(as.numeric(as.character(substr(names(land), 2,5))) == yr) 
  AnnualC.grd = addLayer(AnnualC.grd, calc(land[[idx]], mean, na.rm = TRUE)) 
 
 print(paste("Finish calculating Year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}

#######################################################
########calculate correlation from here
###################################################
pts.sp = Ex.pts.all(mask.1d)

AnnualC.df = extract(AnnualC.grd, pts.sp)
# AnnualNDVI.df = extract(AnnualNDVI.grd, pts.sp)
AnnualNDVI.df = extract(GrowingNDVI.grd, pts.sp)

data.df = cbind(AnnualC.df, AnnualNDVI.df)
data.cor = apply(data.df, 1, COR.test)

cor.grd = Point2raster(data.cor[1,], mask.1d)
pvalue.grd = Point2raster(data.cor[2,], mask.1d)
pvalue.grd = pvalue.grd <0.05
pvalue.grd[pvalue.grd ==0] = NA

pvalue.sp = Ex.pts.all.nonNA(pvalue.grd)



## plot
proj_polar = "+proj=stere +lat_0=90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
proj.sinu = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

country.shp = readOGR("R:/zhihua/AmazonFire/country_polygon", "country")
centroids <- getSpPPolygonsLabptSlots(country.shp)
country.shp@data$Lat = centroids[,2]
country.shp@data$Long = centroids[,1]
country.shp = country.shp[which(country.shp@data$Lat > 25), ]

country.shp.polar <- spTransform(country.shp, CRS(proj_polar))  #projection: geographic projection


r1.grd = cor.grd

Q1 =  quantile(as.matrix(r1.grd),probs = c(0.005,0.995), na.rm = TRUE)
min_ = Q1[1]
max_ = Q1[2]

r1.range = c(-max(abs(min_),abs(max_)),max(abs(min_),abs(max_))) 
r1.grd[r1.grd > r1.range[2]] = r1.range[2]
r1.grd[r1.grd < r1.range[1]] = r1.range[1]

# r1.grd[r1.grd > 15] = 15
# r1.grd[r1.grd < -15] = -15

r1.grd = projectRaster(r1.grd, crs = CRS(proj_polar), method = "ngb")
pvalue.grd1 = projectRaster(pvalue.grd, crs = CRS(proj_polar), method = "ngb")
# pvalue.sp.polar = spTransform(pvalue.sp, CRS(proj_polar))
pvalue.sp.polar = Ex.pts.all.nonNA2(pvalue.grd1)

plot(r1.grd)
plot(country.shp.polar, add = TRUE)

png("R:/zhihua/L4C_NH_analysis_Revision3/Figures/ndvi_carbon_correlation.png",height = 1500, width = 1500, res = 300, units = "px", bg = "transparent")

p.strip <- list(cex=1, lines=2, fontface='bold')

myTheme <- BTCTheme()
myTheme$panel.background$col = 'gray' 

levelplot(r1.grd, par.settings = myTheme,
          zlim = r1.range,
		  margin = FALSE,
		  # main = "C Sink Trend (gC m-2 d-1, 1980-2017)",
		  # main = expression("Net C Uptake Trend"~ "("~gC ~ m^{-2} ~ yr ^{-1}~ yr ^{-1}~ "1980-2017" ~")"),		  
          maxpixels = nrow(r1.grd)*ncol(r1.grd),
          col.regions = rev(colorRampPalette(c("blue", "white", "red"))(255)),
		  # col.regions = colorRampPalette(c("blue", "yellow", "red"))(255),
          # scales=list(x=list(cex=1),y=list(cex=1)),
          # xlab=list(label = "Longtitude", cex=1),ylab=list(label = "Latitude", cex=1),
		  xlab=NULL, ylab=NULL,
		  scales=list(draw=FALSE),
		  layout=c(1, 1),
		  colorkey=list(space="bottom", 
		                at=seq(r1.range[1], r1.range[2], length.out=100), 
						height=1, 
						labels = list(cex = 1.5)),
          par.strip.text=p.strip) + 
		  # latticeExtra::layer(sp.polygons(region.shp, col = "black", lwd = 1.5)) 
		  # latticeExtra::layer(sp.polygons(region.shp0, col = "black", lwd = 2))
		  latticeExtra::layer(sp.points(pvalue.sp.polar, pch=20, col=1, cex = 0.15))+
		  latticeExtra::layer(sp.polygons(country.shp.polar, col = "black", lwd = 0.75))


dev.off()

#```
 theme1 = theme(legend.position = c(0.85,0.75)) +
  theme(legend.text = element_text(size = 18))+
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
  theme(strip.text.x = element_text(size=18), strip.text.y = element_text(size=18)) 
#```{R}


#################################################################
## scatter plot
AnnualC_trend = Trend1(AnnualC.grd)
r5 = AnnualC_trend[[1]]

AnnualNDVI_trend = Trend1(AnnualNDVI.grd)
r50 = AnnualNDVI_trend[[1]]

# GrowingNDVI_trend = Trend1(GrowingNDVI.grd)
# r50 = GrowingNDVI_trend[[1]]

AnnualT.grd = stack("R:/zhihua/PermafrostCCI_v2/AnnualT_1980_2017.tif")
tmean.grd = calc(AnnualT.grd, mean, na.rm = TRUE)

pts.sp = Ex.pts.all.nonNA(mask.1d)

df1 = data.frame(
					  tc = extract(tc.grd1d, pts.sp),
					  pe = extract(perma.grd1d, pts.sp),
					  cor = extract(cor.grd, pts.sp),
					  tmean = extract(tmean.grd, pts.sp),
					  nee_trend = extract(r5, pts.sp),
					  gpp_trend = extract(r50, pts.sp)
					  )


plot(nee_trend~gpp_trend, data = df1)

library(ggplot2)
library(reshape2)
library(ggpmisc)


my.formula <- y ~ x

# plot density plot for cor, colored by tc (< 50 and >= 50)
df1$tc2 = "<50"
df1$tc2[which(df1$tc >= 50)] = ">=50"
df1$tc2 = factor(df1$tc2)

ggplot(df1, aes(x = cor, colour = tc2)) +
  geom_density(lwd = 1.2, linetype = 1)+
   geom_vline(aes(xintercept=0),
            color="blue", linetype="dashed", size=1) + 
  labs(x = "r",
       y = "Density") +		 
  theme_bw()+ theme1 +
  theme(legend.title=element_blank()) + # remove legend title
  guides(size = "none") 
  
ggplot(df1, aes(gpp_trend, nee_trend*365)) +
  geom_point(aes(color = tc)) +
  scale_color_viridis_c()  + 
  labs(x = "NDVI Trend",
       y = expression("Net C Uptake Trend")) +
  # coord_cartesian(xlim = c(0,800), ylim = c(0,800)) +
  # stat_smooth(se = TRUE,
              # colour = "black",
              # size = 0.5,
             # method = "lm") +
  # stat_poly_eq(formula = my.formula, 
               # aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")), 
               # parse = TRUE) + 			 
  theme_bw()+ theme1 +
  theme(legend.title=element_blank()) + # remove legend title
  guides(size = "none")  # remove size
  
#### plot along the 5 gradient

trend.df = df1

 # pe          cor     tmean    nee_trend     gpp_trend

steps = 5 
p1 = seq(0, 100, by = steps)
tc.df = data.frame()
pe.df = data.frame()

for (i in 1:(length(p1)-1)){
idx = which(trend.df$tc >= p1[i] & trend.df$tc < p1[i+1])
if (length(idx) >= 1) {
r0 = as.vector(apply(trend.df[idx,c("cor", "tmean", "nee_trend", "gpp_trend", "pe")], 2, median, na.rm = TRUE))
} else {
r0 = c(NA,NA,NA,NA,NA)
}
r0_1 = c(p1[i]+steps/2, length(idx),r0)
tc.df = rbind(tc.df, r0_1)

}

names(tc.df) <- c("tc", "obs", "cor", "tmean", "nee_trend", "gpp_trend", "pe")
tc.df$nee_trend = tc.df$nee_trend*365

### plot
library(reshape2)
library(ggpmisc)

my.formula <- y ~ x


ggplot(tc.df, 
# ggplot(subset(tc.df, gpp_trend >0.0005), 
aes(gpp_trend, nee_trend)) +
  geom_point(aes(color = tc), size = 5) +
  scale_color_viridis_c()  + 
  labs(x = "NDVI Trend (unit per year)",
       y = expression("Net C Uptake Trend (" ~ gC ~ m^{-2} ~ yr^{-2}~ ")")) +
  # coord_cartesian(xlim = c(0,800), ylim = c(0,800)) +
  stat_smooth(se = TRUE,
              colour = "black",
              size = 0.5,
             method = "lm") +
  stat_poly_eq(formula = my.formula, 
               #aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")), 
               aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 		   
               parse = TRUE) + 			 
  theme_bw()+ theme1 +
  theme(legend.title=element_blank()) + # remove legend title
  guides(size = "none")  # remove size
  
###################################################################################################
## looking at the LUE GPP trend 
############################################################### 

gpp.grd.rs = stack("R:/zhihua/BillSmithGPP/rs_gpp_1982_2015.grd")
gpp.grd.rs = crop(gpp.grd.rs, extent(-180,180, 50,75))
gpp.grd.rs[gpp.grd.rs < 0] = 0

gpp.grd.rs = aggregate(gpp.grd.rs, fact = 2, mean, na.rm = TRUE)
gpp.grd.rs[gpp.grd.rs < 0] = 0

# crop to northern land

#####################################################################################
## plot GPP trend
land = gpp.grd.rs*mask.1d

M = c(paste0("M0", 1:9), paste0("M", 10:12) )
names(land) <- paste0("Y",rep(c(1982:2015), each = 12), rep(M, length(c(1982:2015))))

## calculate T change for each season: Spring [5,6], Summer [7, 8], Fall [9,10], Winter [11-4], Annual [1-12]

SpringGPP.grd = stack()
# SummerC.grd = stack()
FallGPP.grd = stack()
WinterGPP.grd = stack()
AnnualGPP.grd = stack()

for(yr in 1982:2015){
  
  idx = which((as.numeric(as.character(substr(names(land), 7,8))) >= 5 & as.numeric(as.character(substr(names(land), 7,8))) <= 8) &
             as.numeric(as.character(substr(names(land), 2,5))) == yr) 
  SpringGPP.grd = addLayer(SpringGPP.grd, calc(land[[idx]], mean, na.rm = TRUE)) 

   idx = which((as.numeric(as.character(substr(names(land), 7,8))) >= 9 & as.numeric(as.character(substr(names(land), 7,8))) <= 10) &
             as.numeric(as.character(substr(names(land), 2,5))) == yr) 
  FallGPP.grd = addLayer(FallGPP.grd, calc(land[[idx]], mean, na.rm = TRUE)) 
  
   idx = which(((as.numeric(as.character(substr(names(land), 7,8))) >= 1 & as.numeric(as.character(substr(names(land), 7,8))) <= 4) |  
                (as.numeric(as.character(substr(names(land), 7,8))) >= 11 & as.numeric(as.character(substr(names(land), 7,8))) <= 12))&
             as.numeric(as.character(substr(names(land), 2,5))) == yr) 
  WinterGPP.grd = addLayer(WinterGPP.grd, calc(land[[idx]], mean, na.rm = TRUE)) 
  
     idx = which(as.numeric(as.character(substr(names(land), 2,5))) == yr) 
  AnnualGPP.grd = addLayer(AnnualGPP.grd, calc(land[[idx]], mean, na.rm = TRUE)) 
 

 print(paste("Finish calculating Year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}


#######################################################
########calculate correlation from here
###################################################
pts.sp = Ex.pts.all(mask.1d)

AnnualC.df = extract(AnnualC.grd, pts.sp)
AnnualNDVI.df = extract(AnnualGPP.grd, pts.sp)

data.df = cbind(AnnualC.df, AnnualNDVI.df)
data.cor = apply(data.df, 1, COR.test)

cor.grd = Point2raster(data.cor[1,], mask.1d)
pvalue.grd = Point2raster(data.cor[2,], mask.1d)
pvalue.grd = pvalue.grd <0.05
pvalue.grd[pvalue.grd ==0] = NA

pvalue.sp = Ex.pts.all.nonNA(pvalue.grd)

## plot
proj_polar = "+proj=stere +lat_0=90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
proj.sinu = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

country.shp = readOGR("R:/zhihua/AmazonFire/country_polygon", "country")
centroids <- getSpPPolygonsLabptSlots(country.shp)
country.shp@data$Lat = centroids[,2]
country.shp@data$Long = centroids[,1]
country.shp = country.shp[which(country.shp@data$Lat > 25), ]

country.shp.polar <- spTransform(country.shp, CRS(proj_polar))  #projection: geographic projection


r1.grd = cor.grd

Q1 =  quantile(as.matrix(r1.grd),probs = c(0.005,0.995), na.rm = TRUE)
min_ = Q1[1]
max_ = Q1[2]

r1.range = c(-max(abs(min_),abs(max_)),max(abs(min_),abs(max_))) 
r1.grd[r1.grd > r1.range[2]] = r1.range[2]
r1.grd[r1.grd < r1.range[1]] = r1.range[1]

# r1.grd[r1.grd > 15] = 15
# r1.grd[r1.grd < -15] = -15

r1.grd = projectRaster(r1.grd, crs = CRS(proj_polar), method = "ngb")
pvalue.grd1 = projectRaster(pvalue.grd, crs = CRS(proj_polar), method = "ngb")
# pvalue.sp.polar = spTransform(pvalue.sp, CRS(proj_polar))
pvalue.sp.polar = Ex.pts.all.nonNA2(pvalue.grd1)

plot(r1.grd)
plot(country.shp.polar, add = TRUE)

png("R:/zhihua/L4C_NH_analysis_Revision3/Figures/ndvi_carbon_correlation_gpp.png",height = 1500, width = 1500, res = 300, units = "px", bg = "transparent")

p.strip <- list(cex=1, lines=2, fontface='bold')

myTheme <- BTCTheme()
myTheme$panel.background$col = 'gray' 

levelplot(r1.grd, par.settings = myTheme,
          zlim = r1.range,
		  margin = FALSE,	  
          maxpixels = nrow(r1.grd)*ncol(r1.grd),
          col.regions = rev(colorRampPalette(c("blue", "white", "red"))(255)),
		  xlab=NULL, ylab=NULL,
		  scales=list(draw=FALSE),
		  layout=c(1, 1),
		  colorkey=list(space="bottom", 
		                at=seq(r1.range[1], r1.range[2], length.out=100), 
						height=1, 
						labels = list(cex = 1.5)),
          par.strip.text=p.strip) + 
		  latticeExtra::layer(sp.points(pvalue.sp.polar, pch=20, col=1, cex = 0.25))+
		  latticeExtra::layer(sp.polygons(country.shp.polar, col = "black", lwd = 0.75))


dev.off()


#################################################################
## scatterplot
AnnualC_trend = Trend1(AnnualC.grd)
r5 = AnnualC_trend[[1]]

AnnualNDVI_trend = Trend1(AnnualGPP.grd)
r50 = AnnualNDVI_trend[[1]]

AnnualT.grd = stack("R:/zhihua/PermafrostCCI_v2/AnnualT_1980_2017.tif")
tmean.grd = calc(AnnualT.grd, mean, na.rm = TRUE)

pts.sp = Ex.pts.all.nonNA(mask.1d)

df1 = data.frame(
					  tc = extract(tc.grd1d, pts.sp),
					  pe = extract(perma.grd1d, pts.sp),
					  cor = extract(cor.grd, pts.sp),
					  tmean = extract(tmean.grd, pts.sp),
					  nee_trend = extract(r5, pts.sp),
					  gpp_trend = extract(r50, pts.sp)
					  )


plot(nee_trend~gpp_trend, data = df1)

library(ggplot2)
library(reshape2)
library(ggpmisc)

my.formula <- y ~ x

# plot density plot for cor, colored by tc (< 50 and >= 50)
df1$tc2 = "<50"
df1$tc2[which(df1$tc >= 50)] = ">=50"
df1$tc2 = factor(df1$tc2)

ggplot(df1, aes(x = cor, colour = tc2)) +
  geom_density(lwd = 1.2, linetype = 1)+
   geom_vline(aes(xintercept=0),
            color="blue", linetype="dashed", size=1) + 
  labs(x = "r",
       y = "Density") + 
  theme_bw()+ theme1 +
  theme(legend.title=element_blank()) + # remove legend title
  guides(size = "none") 
  
ggplot(df1, aes(gpp_trend, nee_trend*365)) +
  geom_point(aes(color = tc)) +
  scale_color_viridis_c()  + 
  labs(x = "GPP Trend",
       y = expression("Net C Uptake Trend")) +
  # coord_cartesian(xlim = c(0,800), ylim = c(0,800)) +
  # stat_smooth(se = TRUE,
              # colour = "black",
              # size = 0.5,
             # method = "lm") +
  # stat_poly_eq(formula = my.formula, 
               # aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")), 
               # parse = TRUE) + 			 
  theme_bw()+ theme1 +
  theme(legend.title=element_blank()) + # remove legend title
  guides(size = "none")  # remove size
  


#### plot along the 5 gradient

trend.df = df1

steps = 5 
p1 = seq(0, 100, by = steps)
tc.df = data.frame()
pe.df = data.frame()

for (i in 1:(length(p1)-1)){
idx = which(trend.df$tc >= p1[i] & trend.df$tc < p1[i+1])
if (length(idx) >= 1) {
r0 = as.vector(apply(trend.df[idx,c("cor", "tmean", "nee_trend", "gpp_trend", "pe")], 2, median, na.rm = TRUE))
} else {
r0 = c(NA,NA,NA,NA,NA)
}
r0_1 = c(p1[i]+steps/2, length(idx),r0)
tc.df = rbind(tc.df, r0_1)

}

names(tc.df) <- c("tc", "obs", "cor", "tmean", "nee_trend", "gpp_trend", "pe")
tc.df$nee_trend = tc.df$nee_trend*365

### plot
library(reshape2)
library(ggpmisc)

my.formula <- y ~ x


ggplot(tc.df, 
aes(gpp_trend*365, nee_trend)) +
  geom_point(aes(color = tc), size = 5) +
  scale_color_viridis_c()  + 
  labs(x = expression("GPP Trend (" ~ gC ~ m^{-2} ~ yr^{-2}~ ")"),
       y = expression("Net C Uptake Trend (" ~ gC ~ m^{-2} ~ yr^{-2}~ ")")) +
  # coord_cartesian(xlim = c(0,800), ylim = c(0,800)) +
  stat_smooth(se = TRUE,
              colour = "black",
              size = 0.5,
             method = "lm") +
  stat_poly_eq(formula = my.formula, 
               # aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")), 
               aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 		   
               parse = TRUE) + 			 
  theme_bw()+ theme1 +
  theme(legend.title=element_blank()) + # remove legend title
  guides(size = "none")  # remove size
  
 
 
