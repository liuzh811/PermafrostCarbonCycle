# 4/25/2022
# calculate trend for different regions

# (2) try different P classification. 
# a) P, DisconP, and NonP; 
# (b) a + % tree cover from Song et al (2018 Nature); 
# (3) as a sup analysis, used ESA CCI %P (0-33%, 33% - 66%, 66-100%) + %TreeCover (0-33%, 33% - 66%, 66-100%)  

# based on ..Revision2/ncc_lc.r

rm(list = ls())

# load("F:/L4C_NH_analysis/Explanation.TrendyAndRS_2.RData")

## plot the net carbon change rate and percentage of permafrost
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

NDVIGrowing.grd = stack("R:/zhihua/L4C_NH_analysis_Revision/NDVIGrowing.tif")
NDVIannual.grd = stack("R:/zhihua/L4C_NH_analysis_Revision/NDVIAnnual.tif")

NDVIGrowing.grd = calc(NDVIGrowing.grd, mean, na.rm = TRUE)
NDVIannual.grd = calc(NDVIannual.grd, mean, na.rm = TRUE)

# aggregate to 1 degree
NDVIGrowing.grd = aggregate(NDVIGrowing.grd, fact = 4, fun = "mean")
NDVIAnnual.grd = aggregate(NDVIannual.grd, fact = 4, fun = "mean")

mask.1d[NDVIAnnual.grd < 0.1] = NA

## read into land cover and % of permafrost extent
tc.grd = raster("R:/zhihua/L4C_NH_analysis_Revision3/treecover_mean_1983-2016.tif")
perma.grd = raster("R:/zhihua/L4C_NH_analysis_Revision3/permaf_ext_mean_1997_2019_geo2.tif")
sv.grd = raster("R:/zhihua/L4C_NH_analysis_Revision3/shortvegetation_mean_1983-2016.tif")
mod_tc.grd = raster("R:/zhihua/L4C_NH_analysis_Revision3/modis_tc_cmg.tif")

tc.grd1d = aggregate(tc.grd, fact=20, fun=mean, na.rm = TRUE)
sv.grd1d = aggregate(sv.grd, fact=20, fun=mean, na.rm = TRUE)
perma.grd1d = aggregate(perma.grd, fact=20, fun=mean, na.rm = TRUE)
mod_tc.grd1d = aggregate(mod_tc.grd, fact=20, fun=mean, na.rm = TRUE)

# reclassify permafrost into NoP, DisconP, and ConP
frost.shp = readOGR("O:/Zhihua/E_drive_backup/NTSGdata/Permafrost", "perfafrost")
frost.grd = readGDAL("O:/Zhihua/E_drive_backup/NTSGdata/Permafrost/permafrost5k")
frost.grd = raster(frost.grd)
projection(frost.grd) <- projection(frost.shp)

# 5: continuous, 4: spiratic; 3: discontinuous 2: isolated; 1: no Permafrost
frost.grd2 = frost.grd
frost.grd2[is.na(frost.grd)] = NA
frost.grd2[frost.grd == 1 | frost.grd == 2] = 2 #no permafrost
frost.grd2[frost.grd <= 4 & frost.grd >= 3] = 3 #DisconP
frost.grd2[frost.grd == 5] = 4 #ConP

frost.grd2 = projectRaster(frost.grd2, crs = CRS(proj.geo), res = 0.25, method = "ngb")
frost.grd2 = crop(frost.grd2, extent(-180,180, 49,76))
frost.grd2 = aggregate(frost.grd2, fact = 4, fun = modal, na.rm = TRUE)

frost.grd2 = frost.grd2*mask.1d
extent(frost.grd2) <- extent(-180,180, 50,75)

########### read into ACIs, and calculate trend

ct.nee = stack("R:/zhihua/InversionCO2_GCB2020/NEE.ct2019b.2000.2018.grd")
cte.nee = stack("R:/zhihua/InversionCO2_GCB2020/NEE.cte.2000.2019.grd")
jena.nee = stack("R:/zhihua/InversionCO2_GCB2020/NEE.jena.1976.2017.grd") # not updated
jena.nee2 = stack("R:/zhihua/InversionCO2_GCB2020/NEE.jena2.1985.2017.grd") # not updated
jams.nee = stack("R:/zhihua/InversionCO2_GCB2020/NEE.jams.1996.2017.grd")# not updated
# cams.nee = stack("F:/InversionCO2_GCB2020/NEE.cams.1979.2019.grd")
nism.nee = stack("R:/zhihua/InversionCO2_GCB2020/NEE.nism.1990.2019.grd")
cams.nee = stack("R:/zhihua/InversionCO2_GCB2020/NEE.cams.1979.2019_1d.grd")

ensemble.nee = stack("R:/zhihua/InversionCO2_GCB2020/NEE.monthly.mean1976.2017.grd")

## calculate T change for each season: Spring [5,6], Summer [7, 8], Fall [9,10], Winter [11-4], Annual [1-12]
#land.grd = stack("R:/zhihua/InversionCO2/NEE.monthly.mean1976.2017.grd")
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

for(yr in 1980:2017){
  
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

############### PART1: for different regions:: NoP, DisconP, vs P      ###############################

d1.df01 = data.frame(cbind(t(data.frame(zonal(SpringC.grd, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
                   # t(data.frame(zonal(SummerC.grd, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
				   t(data.frame(zonal(FallC.grd, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
				   t(data.frame(zonal(WinterC.grd, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
				   t(data.frame(zonal(AnnualC.grd, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
                          c(1980:2017)))


names(d1.df01) = c("SpringNoP", "SpringDisconP", "SpringP", 
				   "FallNoP", "FallDisconP", "FallP",
				   "WinterNoP", "WinterDisconP", "WinterP",
				   "AnnualNoP", "AnnualDisconP", "AnnualP",
				   "Year")
				   
ct.df = data.frame(t(zonal(ct.nee*mask.1d, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(ct.df) <- c("AnnualNoP", "AnnualDisconP", "AnnualP")
ct.df$Year = as.numeric(as.character(substr(names(ct.nee), 2,5)))
ct.df = aggregate(cbind(AnnualNoP, AnnualDisconP,AnnualP) ~ Year, FUN=mean, na.rm = TRUE,data = ct.df)

cte.df = data.frame(t(zonal(cte.nee*mask.1d, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(cte.df) <- c("AnnualNoP", "AnnualDisconP", "AnnualP")
cte.df$Year = as.numeric(as.character(substr(names(cte.nee), 2,5)))
cte.df = aggregate(cbind(AnnualNoP, AnnualDisconP,AnnualP) ~ Year, FUN=mean, na.rm = TRUE,data = cte.df)

jena2.df = data.frame(t(zonal(jena.nee2*mask.1d, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(jena2.df) <- c("AnnualNoP", "AnnualDisconP", "AnnualP")
jena2.df$Year = as.numeric(as.character(substr(names(jena.nee2), 2,5)))
jena2.df = aggregate(cbind(AnnualNoP, AnnualDisconP,AnnualP) ~ Year, FUN=mean, na.rm = TRUE,data = jena2.df)

jena.df = data.frame(t(zonal(jena.nee*mask.1d, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(jena.df) <- c("AnnualNoP", "AnnualDisconP", "AnnualP")
jena.df$Year = as.numeric(as.character(substr(names(jena.nee), 2,5)))
jena.df = aggregate(cbind(AnnualNoP, AnnualDisconP,AnnualP) ~ Year, FUN=mean, na.rm = TRUE,data = jena.df)

jams.df = data.frame(t(zonal(jams.nee*mask.1d, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(jams.df) <- c("AnnualNoP", "AnnualDisconP", "AnnualP")
jams.df$Year = as.numeric(as.character(substr(names(jams.nee), 2,5)))
jams.df = aggregate(cbind(AnnualNoP, AnnualDisconP,AnnualP) ~ Year, FUN=mean, na.rm = TRUE,data = jams.df)

cams.df = data.frame(t(zonal(cams.nee*mask.1d, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(cams.df) <- c("AnnualNoP", "AnnualDisconP", "AnnualP")
cams.df$Year = as.numeric(as.character(substr(names(cams.nee), 2,5)))
cams.df = aggregate(cbind(AnnualNoP, AnnualDisconP,AnnualP) ~ Year, FUN=mean, na.rm = TRUE,data = cams.df)

ensemble.df = data.frame(t(zonal(ensemble.nee*mask.1d, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(ensemble.df) <- c("AnnualNoP", "AnnualDisconP", "AnnualP")
ensemble.df$Year = as.numeric(as.character(substr(names(ensemble.nee), 2,5)))
ensemble.df = aggregate(cbind(AnnualNoP, AnnualDisconP,AnnualP) ~ Year, FUN=mean, na.rm = TRUE,data = ensemble.df)


# plot permafrost and non-permrost
dummy.df = data.frame(dummy = NA, Year = 1980:2017)
dat.df0001 = merge(ct.df, dummy.df, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "ctNoP", "ctDisconP", "ctP")

dat.df0001 = merge(cte.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "cteNoP", "cteDisconP", "cteP")

dat.df0001 = merge(jena2.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "jena2NoP", "jena2DisconP", "jena2P")

dat.df0001 = merge(jena.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "jenaNoP", "jenaDisconP", "jenaP")

dat.df0001 = merge(jams.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "jamsNoP", "jamsDisconP", "jamsP")

dat.df0001 = merge(cams.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "camsNoP", "camsDisconP", "camsP")

dat.df0001 = merge(ensemble.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "EnsemNoP", "EnsemDisconP", "EnsemP")

sde = function(x){
x1 = sd(x, na.rm = TRUE)
n = length(which(!is.na(x)))
return(x1/(sqrt(n)))
}

dat.df0001$EnsemNoPsd = apply(dat.df0001[,c("ctNoP", "cteNoP", "jena2NoP","jenaNoP", "jamsNoP", "camsNoP")], 1, sd, na.rm = TRUE)
dat.df0001$EnsemDisconPsd = apply(dat.df0001[,c("ctDisconP", "cteDisconP", "jena2DisconP","jenaDisconP", "jamsDisconP", "camsDisconP")], 1, sd, na.rm = TRUE)
dat.df0001$EnsemPsd = apply(dat.df0001[,c("ctP", "cteP", "jena2P","jenaP", "jamsP", "camsP")], 1, sd, na.rm = TRUE)

## plot 
png("R:/zhihua/L4C_NH_analysis_Revision3/Figures/CsinkTrend_3_regions_Permafrost.png",height = 2100, width = 2100, res = 300, units = "px", bg = "transparent")

par(mar=c(4,5,5,0)+.1)

# P
plot(EnsemP~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017, main = "Blue = Continous Permafrost  Region\n Green = Discontinous Permafrost  Region\n Red = Non-Permafrost Region",
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "blue",lwd = 3, 
								  xlab = "Year", ylab = expression("Net C Uptake (" ~ gC ~ m^{-2} ~ yr ^{-1}~ ")"),
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)

# add uncertaines
polygon(c(dat.df0001$Year, rev(dat.df0001$Year)), 
		c(dat.df0001$EnsemP+dat.df0001$EnsemPsd, rev(dat.df0001$EnsemP-dat.df0001$EnsemPsd)), 
        col=rgb(0, 0, 1,0.25),
		border = NA)
	

# DisconP
par(new=TRUE)

plot(EnsemDisconP~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017,
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "green",lwd = 3, 
								  xlab = "Year", ylab = expression("Net C Uptake (" ~ gC ~ m^{-2} ~ yr ^{-1}~ ")"),
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)

# add uncertaines
polygon(c(dat.df0001$Year, rev(dat.df0001$Year)), 
		c(dat.df0001$EnsemDisconP+dat.df0001$EnsemDisconPsd, rev(dat.df0001$EnsemDisconP-dat.df0001$EnsemDisconPsd)), 
        col=rgb(0, 1, 0,0.25),
		border = NA)
	

# NonP
par(new=TRUE)
plot(EnsemNoP~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017,
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "red",lwd = 3, 
								  ylab = "", 
								  xlab = "", 
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)


# add uncertaines
polygon(c(dat.df0001$Year, rev(dat.df0001$Year)), 
		c(dat.df0001$EnsemNoP+dat.df0001$EnsemNoPsd, rev(dat.df0001$EnsemNoP-dat.df0001$EnsemNoPsd)), 
        col=rgb(1, 0, 0,0.25),
		border = NA)
	


## add statiscis 


lmf = function(x,y){

lm1 = lm(y~x)

newx = 1980:2017
conf_interval <- predict(lm1, newdata=data.frame(Year=newx), interval="confidence",level = 0.95)

modsum = summary(lm1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
sen = modsum$coefficients[2,1]
rp = vector('expression',1)
sd1 = modsum$coefficients[2,2]

# rp[1] = substitute(expression(""~ delta == MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

# rp[1] = substitute(expression("Trend:"~ MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

rp[1] = substitute(expression(""~ MYVALUE1 ~ "±" ~ MYVALUE2 ~ gC ~ m^{-2} ~ yr ^{-2}~""), 
		list(MYVALUE1 = format(sen, digits = 2), MYVALUE2 = format(sd1, digits = 2)))[2]

return(list(newx, conf_interval,lm1,rp))
}
				   
x = dat.df0001$Year
y = dat.df0001$EnsemP
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="blue", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="blue", lty=2)
abline(lm1[[3]], lwd = 3,col = "blue")
legend('bottomleft', legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "blue")
# text(x = 2000, y = 0.2*365, labels = expression("Permafrost: 1.00 ± 0.10"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "blue")


x = dat.df0001$Year
y = dat.df0001$EnsemDisconP
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="green", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="green", lty=2)
abline(lm1[[3]], lwd = 3,col = "green")
legend('topleft', legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "green")
# text(x = 2000, y = 0.22*365, labels = expression("Non-Permafrost: 0.22 ± 0.13"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "red")

				   
x = dat.df0001$Year
y = dat.df0001$EnsemNoP
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="red", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="red", lty=2)
abline(lm1[[3]], lwd = 3,col = "red")
legend('topright', legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "red")
# text(x = 2000, y = 0.22*365, labels = expression("Non-Permafrost: 0.22 ± 0.13"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "red")

# inset a subset to seasonsity
d1.df010 = d1.df01

par(new = TRUE)
par(fig = c(0.55, 1, 0.1, 0.5))
v1 = c()
v2 = c()
for(i in 1:12){

lm1 = lm(d1.df010[,i]~d1.df010[,13])
modsum = summary(lm1)
sen = modsum$coefficients[2,1]
v1 = c(v1, sen)

sd1 = modsum$coefficients[2,2]
v2 = c(v2, sd1)

}

v1.df = rbind(v1[c(1,4,7)],v1[c(2,5,8)],v1[c(3,6,9)])

rownames(v1.df) = c("NonP", "DisconP", "P")
colnames(v1.df) = c("EGS", "LGS", "Win")

# plot with sd
v2.df = rbind(v2[c(1,4,7)],v2[c(2,5,8)],v2[c(3,6,9)])

rownames(v2.df) = c("NonP", "DisconP", "P")
colnames(v2.df) = c("EGS", "LGS", "Win")

# us ggplot 2
# http://www.sthda.com/english/wiki/ggplot2-error-NDVIrs-quick-start-guide-r-software-and-data-visualization

v3.df = data.frame(mn = c(v1.df[,1], v1.df[,2], v1.df[,3]), 
                   sd1 = c(v2.df[,1], v2.df[,2], v2.df[,3]),
				   Season = c(rep("EGS",3),rep("LGS",3),rep("Win",3)),
				   Region = rep(c("NonP", "DisconP", "P"),3))
				   

v3.df$Season <- factor(v3.df$Season, levels = c("EGS", "LGS", "Win"))
v3.df$Region <- factor(v3.df$Region, levels = c("NonP", "DisconP", "P"))

require(ggplot2)
p <- ggplot(v3.df, aes(x=Season, y=mn, fill=Region)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mn-sd1, ymax=mn+sd1), width=.2,
                 position=position_dodge(.9)) + 
				     theme_bw() + 
  theme(legend.text = element_text(size = 18))+
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
 xlab("") + 
  ylab("") + 
 # ylab(expression(""~ gC ~ m^{-2} ~ yr ^{-2}~"")) + 
    theme(strip.text.x = element_text(size=18), strip.text.y = element_text(size=18)) + 
  theme(legend.position = "none") + 
  scale_fill_manual(values=c("red","green","blue"))

# print(p)


p = p +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

# overlay ggplot2 figure
require(grid)
par(new=TRUE)

print(p, vp=viewport(.78, .25, .45, .3)) #viewport, first two is the x,y coor; last 2 is the inset size

dev.off()

############### PART2: for different tree cover :: No tree (<30%), tree-short mix (30-50), tree (> 50)      ###############################

tc.grd1d2 = tc.grd1d*mask.1d
tc.grd1d3 = tc.grd1d2
tc.grd1d2[tc.grd1d3 < 30] = 3
tc.grd1d2[tc.grd1d3 >= 30 & tc.grd1d3 < 50] = 2
tc.grd1d2[tc.grd1d3 >= 50] = 1

tc.grd1d2 = tc.grd1d2*mask.1d

frost.grd2 = tc.grd1d2

d1.df01 = data.frame(cbind(t(data.frame(zonal(SpringC.grd, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
                   # t(data.frame(zonal(SummerC.grd, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
				   t(data.frame(zonal(FallC.grd, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
				   t(data.frame(zonal(WinterC.grd, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
				   t(data.frame(zonal(AnnualC.grd, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
                          c(1980:2017)))


names(d1.df01) = c("SpringNoP", "SpringDisconP", "SpringP", 
				   "FallNoP", "FallDisconP", "FallP",
				   "WinterNoP", "WinterDisconP", "WinterP",
				   "AnnualNoP", "AnnualDisconP", "AnnualP",
				   "Year")
				   
ct.df = data.frame(t(zonal(ct.nee*mask.1d, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(ct.df) <- c("AnnualNoP", "AnnualDisconP", "AnnualP")
ct.df$Year = as.numeric(as.character(substr(names(ct.nee), 2,5)))
ct.df = aggregate(cbind(AnnualNoP, AnnualDisconP,AnnualP) ~ Year, FUN=mean, na.rm = TRUE,data = ct.df)

cte.df = data.frame(t(zonal(cte.nee*mask.1d, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(cte.df) <- c("AnnualNoP", "AnnualDisconP", "AnnualP")
cte.df$Year = as.numeric(as.character(substr(names(cte.nee), 2,5)))
cte.df = aggregate(cbind(AnnualNoP, AnnualDisconP,AnnualP) ~ Year, FUN=mean, na.rm = TRUE,data = cte.df)

jena2.df = data.frame(t(zonal(jena.nee2*mask.1d, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(jena2.df) <- c("AnnualNoP", "AnnualDisconP", "AnnualP")
jena2.df$Year = as.numeric(as.character(substr(names(jena.nee2), 2,5)))
jena2.df = aggregate(cbind(AnnualNoP, AnnualDisconP,AnnualP) ~ Year, FUN=mean, na.rm = TRUE,data = jena2.df)

jena.df = data.frame(t(zonal(jena.nee*mask.1d, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(jena.df) <- c("AnnualNoP", "AnnualDisconP", "AnnualP")
jena.df$Year = as.numeric(as.character(substr(names(jena.nee), 2,5)))
jena.df = aggregate(cbind(AnnualNoP, AnnualDisconP,AnnualP) ~ Year, FUN=mean, na.rm = TRUE,data = jena.df)

jams.df = data.frame(t(zonal(jams.nee*mask.1d, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(jams.df) <- c("AnnualNoP", "AnnualDisconP", "AnnualP")
jams.df$Year = as.numeric(as.character(substr(names(jams.nee), 2,5)))
jams.df = aggregate(cbind(AnnualNoP, AnnualDisconP,AnnualP) ~ Year, FUN=mean, na.rm = TRUE,data = jams.df)

cams.df = data.frame(t(zonal(cams.nee*mask.1d, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(cams.df) <- c("AnnualNoP", "AnnualDisconP", "AnnualP")
cams.df$Year = as.numeric(as.character(substr(names(cams.nee), 2,5)))
cams.df = aggregate(cbind(AnnualNoP, AnnualDisconP,AnnualP) ~ Year, FUN=mean, na.rm = TRUE,data = cams.df)

ensemble.df = data.frame(t(zonal(ensemble.nee*mask.1d, frost.grd2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(ensemble.df) <- c("AnnualNoP", "AnnualDisconP", "AnnualP")
ensemble.df$Year = as.numeric(as.character(substr(names(ensemble.nee), 2,5)))
ensemble.df = aggregate(cbind(AnnualNoP, AnnualDisconP,AnnualP) ~ Year, FUN=mean, na.rm = TRUE,data = ensemble.df)


# plot permafrost and non-permrost
dummy.df = data.frame(dummy = NA, Year = 1980:2017)
dat.df0001 = merge(ct.df, dummy.df, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "ctNoP", "ctDisconP", "ctP")

dat.df0001 = merge(cte.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "cteNoP", "cteDisconP", "cteP")

dat.df0001 = merge(jena2.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "jena2NoP", "jena2DisconP", "jena2P")

dat.df0001 = merge(jena.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "jenaNoP", "jenaDisconP", "jenaP")

dat.df0001 = merge(jams.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "jamsNoP", "jamsDisconP", "jamsP")

dat.df0001 = merge(cams.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "camsNoP", "camsDisconP", "camsP")

dat.df0001 = merge(ensemble.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "EnsemNoP", "EnsemDisconP", "EnsemP")

sde = function(x){
x1 = sd(x, na.rm = TRUE)
n = length(which(!is.na(x)))
return(x1/(sqrt(n)))
}

dat.df0001$EnsemNoPsd = apply(dat.df0001[,c("ctNoP", "cteNoP", "jena2NoP","jenaNoP", "jamsNoP", "camsNoP")], 1, sd, na.rm = TRUE)
dat.df0001$EnsemDisconPsd = apply(dat.df0001[,c("ctDisconP", "cteDisconP", "jena2DisconP","jenaDisconP", "jamsDisconP", "camsDisconP")], 1, sd, na.rm = TRUE)
dat.df0001$EnsemPsd = apply(dat.df0001[,c("ctP", "cteP", "jena2P","jenaP", "jamsP", "camsP")], 1, sd, na.rm = TRUE)

## plot 
png("R:/zhihua/L4C_NH_analysis_Revision3/Figures/CsinkTrend_3_regions_TC.png",height = 2100, width = 2100, res = 300, units = "px", bg = "transparent")

par(mar=c(4,5,5,0)+.1)

# P
plot(EnsemP~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017, main = "Red: TC > 50; Green: TC = [30, 50); Blue: TC < 30",
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "blue",lwd = 3, 
								  xlab = "Year", ylab = expression("Net C Uptake (" ~ gC ~ m^{-2} ~ yr ^{-1}~ ")"),
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)

# add uncertaines
polygon(c(dat.df0001$Year, rev(dat.df0001$Year)), 
		c(dat.df0001$EnsemP+dat.df0001$EnsemPsd, rev(dat.df0001$EnsemP-dat.df0001$EnsemPsd)), 
        col=rgb(0, 0, 1,0.25),
		border = NA)
	

# DisconP
par(new=TRUE)

plot(EnsemDisconP~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017,
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "green",lwd = 3, 
								  xlab = "Year", ylab = expression("Net C Uptake (" ~ gC ~ m^{-2} ~ yr ^{-1}~ ")"),
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)

# add uncertaines
polygon(c(dat.df0001$Year, rev(dat.df0001$Year)), 
		c(dat.df0001$EnsemDisconP+dat.df0001$EnsemDisconPsd, rev(dat.df0001$EnsemDisconP-dat.df0001$EnsemDisconPsd)), 
        col=rgb(0, 1, 0,0.25),
		border = NA)
	

# NonP
par(new=TRUE)
plot(EnsemNoP~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017,
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "red",lwd = 3, 
								  ylab = "", 
								  xlab = "", 
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)


# add uncertaines
polygon(c(dat.df0001$Year, rev(dat.df0001$Year)), 
		c(dat.df0001$EnsemNoP+dat.df0001$EnsemNoPsd, rev(dat.df0001$EnsemNoP-dat.df0001$EnsemNoPsd)), 
        col=rgb(1, 0, 0,0.25),
		border = NA)
	


## add statiscis 


lmf = function(x,y){

lm1 = lm(y~x)

newx = 1980:2017
conf_interval <- predict(lm1, newdata=data.frame(Year=newx), interval="confidence",level = 0.95)

modsum = summary(lm1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
sen = modsum$coefficients[2,1]
rp = vector('expression',1)
sd1 = modsum$coefficients[2,2]

# rp[1] = substitute(expression(""~ delta == MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

# rp[1] = substitute(expression("Trend:"~ MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

rp[1] = substitute(expression(""~ MYVALUE1 ~ "±" ~ MYVALUE2 ~ gC ~ m^{-2} ~ yr ^{-2}~""), 
		list(MYVALUE1 = format(sen, digits = 2), MYVALUE2 = format(sd1, digits = 2)))[2]

return(list(newx, conf_interval,lm1,rp))
}
				   
x = dat.df0001$Year
y = dat.df0001$EnsemP
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="blue", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="blue", lty=2)
abline(lm1[[3]], lwd = 3,col = "blue")
legend('bottomleft', legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "blue")
# text(x = 2000, y = 0.2*365, labels = expression("Permafrost: 1.00 ± 0.10"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "blue")


x = dat.df0001$Year
y = dat.df0001$EnsemDisconP
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="green", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="green", lty=2)
abline(lm1[[3]], lwd = 3,col = "green")
legend('topleft', legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "green")
# text(x = 2000, y = 0.22*365, labels = expression("Non-Permafrost: 0.22 ± 0.13"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "red")

				   
x = dat.df0001$Year
y = dat.df0001$EnsemNoP
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="red", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="red", lty=2)
abline(lm1[[3]], lwd = 3,col = "red")
legend('topright', legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "red")
# text(x = 2000, y = 0.22*365, labels = expression("Non-Permafrost: 0.22 ± 0.13"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "red")

# inset a subset to seasonsity
d1.df010 = d1.df01

par(new = TRUE)
par(fig = c(0.55, 1, 0.1, 0.5))
v1 = c()
v2 = c()
for(i in 1:12){

lm1 = lm(d1.df010[,i]~d1.df010[,13])
modsum = summary(lm1)
sen = modsum$coefficients[2,1]
v1 = c(v1, sen)

sd1 = modsum$coefficients[2,2]
v2 = c(v2, sd1)

}

v1.df = rbind(v1[c(1,4,7)],v1[c(2,5,8)],v1[c(3,6,9)])

rownames(v1.df) = c("NonP", "DisconP", "P")
colnames(v1.df) = c("EGS", "LGS", "Win")

# plot with sd
v2.df = rbind(v2[c(1,4,7)],v2[c(2,5,8)],v2[c(3,6,9)])

rownames(v2.df) = c("NonP", "DisconP", "P")
colnames(v2.df) = c("EGS", "LGS", "Win")

# us ggplot 2
# http://www.sthda.com/english/wiki/ggplot2-error-NDVIrs-quick-start-guide-r-software-and-data-visualization

v3.df = data.frame(mn = c(v1.df[,1], v1.df[,2], v1.df[,3]), 
                   sd1 = c(v2.df[,1], v2.df[,2], v2.df[,3]),
				   Season = c(rep("EGS",3),rep("LGS",3),rep("Win",3)),
				   Region = rep(c("NonP", "DisconP", "P"),3))
				   

v3.df$Season <- factor(v3.df$Season, levels = c("EGS", "LGS", "Win"))
v3.df$Region <- factor(v3.df$Region, levels = c("NonP", "DisconP", "P"))

require(ggplot2)
p <- ggplot(v3.df, aes(x=Season, y=mn, fill=Region)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mn-sd1, ymax=mn+sd1), width=.2,
                 position=position_dodge(.9)) + 
				     theme_bw() + 
  theme(legend.text = element_text(size = 18))+
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
 xlab("") + 
  ylab("") + 
 # ylab(expression(""~ gC ~ m^{-2} ~ yr ^{-2}~"")) + 
    theme(strip.text.x = element_text(size=18), strip.text.y = element_text(size=18)) + 
  theme(legend.position = "none") + 
  scale_fill_manual(values=c("red","green","blue"))

# print(p)


p = p +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

# overlay ggplot2 figure
require(grid)
par(new=TRUE)

print(p, vp=viewport(.78, .25, .45, .3)) #viewport, first two is the x,y coor; last 2 is the inset size

dev.off()

###########################################################################################################
############### PART3: TC+%P
# for different tree cover :: No tree (<30%), tree-short mix (30-50), tree (> 50)      ###############################

# plot tree density in different P catogary

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

tc.grd1d2 = tc.grd1d*mask.1d

frost.grd2 = frost.grd
frost.grd2[is.na(frost.grd)] = NA
frost.grd2[frost.grd == 1 | frost.grd == 2] = 2 #no permafrost
frost.grd2[frost.grd <= 4 & frost.grd >= 3] = 3 #DisconP
frost.grd2[frost.grd == 5] = 4 #ConP

frost.grd2 = projectRaster(frost.grd2, crs = CRS(proj.geo), res = 0.25, method = "ngb")
frost.grd2 = crop(frost.grd2, extent(-180,180, 49,76))
frost.grd2 = aggregate(frost.grd2, fact = 4, fun = modal, na.rm = TRUE)

frost.grd2 = frost.grd2*mask.1d
extent(frost.grd2) <- extent(-180,180, 50,75)
frost.grd3 = frost.grd2 - 1


dat.df = data.frame(# tc = extract(mod_tc.grd1d, pts.sp),			  
					  tc = extract(tc.grd1d2, pts.sp),
					  pe = extract(frost.grd3, pts.sp)			  
					  )

#ggplot
library(ggplot2)
# Basic box plot
dat.df$pe <- as.factor(dat.df$pe)
levels(dat.df$pe) <- c("NoP", "DisconP", "P")

p <- ggplot(dat.df, aes(x=pe, y=tc, color=pe)) + 
  geom_boxplot() + 
				     theme_bw() + 
  theme(legend.text = element_text(size = 18))+
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
 xlab("Permafrost") + 
  ylab("% Tree Cover") + 
 # ylab(expression(""~ gC ~ m^{-2} ~ yr ^{-2}~"")) + 
    theme(strip.text.x = element_text(size=18), strip.text.y = element_text(size=18)) + 
  theme(legend.position = "none")
p

ggsave("R:/zhihua/L4C_NH_analysis_Revision3/Figures/TC_p.png", width = 5, height = 5, units = "in")


# classify TC
tc.grd1d2 = tc.grd1d*mask.1d
# tc.grd1d2 = mod_tc.grd1d*mask.1d
tc.grd1d3 = tc.grd1d2
tc.grd1d2[tc.grd1d3 < 30] = 3
tc.grd1d2[tc.grd1d3 >= 30 & tc.grd1d3 < 50] = 2
tc.grd1d2[tc.grd1d3 >= 50] = 1

# classify %p
frost.grd2 = frost.grd
frost.grd2[is.na(frost.grd)] = NA
frost.grd2[frost.grd == 1 | frost.grd == 2] = 2 #no permafrost
frost.grd2[frost.grd <= 4 & frost.grd >= 3] = 3 #DisconP
frost.grd2[frost.grd == 5] = 4 #ConP

frost.grd2 = projectRaster(frost.grd2, crs = CRS(proj.geo), res = 0.25, method = "ngb")
frost.grd2 = crop(frost.grd2, extent(-180,180, 49,76))
frost.grd2 = aggregate(frost.grd2, fact = 4, fun = modal, na.rm = TRUE)

frost.grd2 = frost.grd2*mask.1d
extent(frost.grd2) <- extent(-180,180, 50,75)
frost.grd3 = frost.grd2 - 1 #1 = NoP, 2 = DisconP, 3 = P

# combine together
frost.grd00 = frost.grd3*10+tc.grd1d2 ## 
frost.grd00 = frost.grd00*mask.1d
table(as.vector(as.matrix(frost.grd00)))

#11-NoP and TC > 50%
#12-NoP and TC [30-35]
#13-NoP and TC < 30 

#21-DisconP and TC > 50%
#22-DisconP and TC [30-35]
#23-DisconP and TC < 30 

#31-P and TC > 50%
#32-P and TC [30-35]
#33-P and TC < 30 

> table(as.vector(as.matrix(frost.grd00)))

  11   12   13   21   22   23   31   32   33 
 743  362  217  226  354  385   27  302 1519
## compare same tree cover in differ1 P regions
# 


frost.grd1d2 = mask.1d
frost.grd1d2[] = NA 
frost.grd1d2[frost.grd00 == 31] = 1 #Tree In P
frost.grd1d2[frost.grd00 == 32] = 2 #TundraTreeMix In P
frost.grd1d2[frost.grd00 == 33] = 3 #Tundra In P


d1.df01 = data.frame(cbind(t(data.frame(zonal(SpringC.grd, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
                   # t(data.frame(zonal(SummerC.grd, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
				   t(data.frame(zonal(FallC.grd, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
				   t(data.frame(zonal(WinterC.grd, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
				   t(data.frame(zonal(AnnualC.grd, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
                          c(1980:2017)))


ct.df = data.frame(t(zonal(ct.nee*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(ct.df) <- c("AnnualBorealNoP", "AnnualBorealP", "AnnualTundraP")
ct.df$Year = as.numeric(as.character(substr(names(ct.nee), 2,5)))
ct.df = aggregate(cbind(AnnualBorealNoP, AnnualBorealP,AnnualTundraP) ~ Year, FUN=mean, na.rm = TRUE,data = ct.df)

cte.df = data.frame(t(zonal(cte.nee*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(cte.df) <- c("AnnualBorealNoP", "AnnualBorealP", "AnnualTundraP")
cte.df$Year = as.numeric(as.character(substr(names(cte.nee), 2,5)))
cte.df = aggregate(cbind(AnnualBorealNoP, AnnualBorealP,AnnualTundraP) ~ Year, FUN=mean, na.rm = TRUE,data = cte.df)

jena2.df = data.frame(t(zonal(jena.nee2*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(jena2.df) <- c("AnnualBorealNoP", "AnnualBorealP", "AnnualTundraP")
jena2.df$Year = as.numeric(as.character(substr(names(jena.nee2), 2,5)))
jena2.df = aggregate(cbind(AnnualBorealNoP, AnnualBorealP,AnnualTundraP) ~ Year, FUN=mean, na.rm = TRUE,data = jena2.df)

jena.df = data.frame(t(zonal(jena.nee*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(jena.df) <- c("AnnualBorealNoP", "AnnualBorealP", "AnnualTundraP")
jena.df$Year = as.numeric(as.character(substr(names(jena.nee), 2,5)))
jena.df = aggregate(cbind(AnnualBorealNoP, AnnualBorealP,AnnualTundraP) ~ Year, FUN=mean, na.rm = TRUE,data = jena.df)

jams.df = data.frame(t(zonal(jams.nee*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(jams.df) <- c("AnnualBorealNoP", "AnnualBorealP", "AnnualTundraP")
jams.df$Year = as.numeric(as.character(substr(names(jams.nee), 2,5)))
jams.df = aggregate(cbind(AnnualBorealNoP, AnnualBorealP,AnnualTundraP) ~ Year, FUN=mean, na.rm = TRUE,data = jams.df)

cams.df = data.frame(t(zonal(cams.nee*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(cams.df) <- c("AnnualBorealNoP", "AnnualBorealP", "AnnualTundraP")
cams.df$Year = as.numeric(as.character(substr(names(cams.nee), 2,5)))
cams.df = aggregate(cbind(AnnualBorealNoP, AnnualBorealP,AnnualTundraP) ~ Year, FUN=mean, na.rm = TRUE,data = cams.df)

ensemble.df = data.frame(t(zonal(ensemble.nee*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(ensemble.df) <- c("AnnualBorealNoP", "AnnualBorealP", "AnnualTundraP")
ensemble.df$Year = as.numeric(as.character(substr(names(ensemble.nee), 2,5)))
ensemble.df = aggregate(cbind(AnnualBorealNoP, AnnualBorealP,AnnualTundraP) ~ Year, FUN=mean, na.rm = TRUE,data = ensemble.df)


# plot permafrost and non-permrost
dummy.df = data.frame(dummy = NA, Year = 1980:2017)
dat.df0001 = merge(ct.df, dummy.df, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "ctBorealNoP", "ctBorealP", "ctTundraP")

dat.df0001 = merge(cte.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "cteBorealNoP", "cteBorealP", "cteTundraP")

dat.df0001 = merge(jena2.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "jena2BorealNoP", "jena2BorealP", "jena2TundraP")

dat.df0001 = merge(jena.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "jenaBorealNoP", "jenaBorealP", "jenaTundraP")

dat.df0001 = merge(jams.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "jamsBorealNoP", "jamsBorealP", "jamsTundraP")

dat.df0001 = merge(cams.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "camsBorealNoP", "camsBorealP", "camsTundraP")

dat.df0001 = merge(ensemble.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:4)] <- c("Year", "EnsemBorealNoP", "EnsemBorealP", "EnsemTundraP")

sde = function(x){
x1 = sd(x, na.rm = TRUE)
n = length(which(!is.na(x)))
return(x1/(sqrt(n)))
}

dat.df0001$EnsemBorealNoPsd = apply(dat.df0001[,c("ctBorealNoP", "cteBorealNoP", "jena2BorealNoP","jenaBorealNoP", "jamsBorealNoP", "camsBorealNoP")], 1, sd, na.rm = TRUE)
dat.df0001$EnsemBorealPsd = apply(dat.df0001[,c("ctBorealP", "cteBorealP", "jena2BorealP","jenaBorealP", "jamsBorealP", "camsBorealP")], 1, sd, na.rm = TRUE)
dat.df0001$EnsemTundraPsd = apply(dat.df0001[,c("ctTundraP", "cteTundraP", "jena2TundraP","jenaTundraP", "jamsTundraP", "camsTundraP")], 1, sd, na.rm = TRUE)

## plot 
png("R:/zhihua/L4C_NH_analysis_Revision3/Figures/CsinkTrend_3_regions_Permafrost+TC.png",height = 2100, width = 2100, res = 300, units = "px", bg = "transparent")

par(mar=c(4,5,5,0)+.1)

# Tundra in P
plot(EnsemTundraP~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017, main = "Blue = Tundra In Permafrost \n Green = TundraTreeMix In Permafrost\n Red = Tree In NonPermafrost",
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "blue",lwd = 3, 
								  xlab = "Year", ylab = expression("Net C Uptake (" ~ gC ~ m^{-2} ~ yr ^{-1}~ ")"),
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)

# add uncertaines
polygon(c(dat.df0001$Year, rev(dat.df0001$Year)), 
		c(dat.df0001$EnsemTundraP+dat.df0001$EnsemTundraPsd, rev(dat.df0001$EnsemTundraP-dat.df0001$EnsemTundraPsd)), 
        col=rgb(0, 0, 1,0.25),
		border = NA)
	

# Boreal in P
par(new=TRUE)

plot(EnsemBorealP~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017,
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "green",lwd = 3, 
								  xlab = "Year", ylab = expression("Net C Uptake (" ~ gC ~ m^{-2} ~ yr ^{-1}~ ")"),
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)

# add uncertaines
polygon(c(dat.df0001$Year, rev(dat.df0001$Year)), 
		c(dat.df0001$EnsemBorealP+dat.df0001$EnsemBorealPsd, rev(dat.df0001$EnsemBorealP-dat.df0001$EnsemBorealPsd)), 
        col=rgb(0, 1, 0,0.25),
		border = NA)
	



# Boreal in NonP
par(new=TRUE)
plot(EnsemBorealNoP~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017,
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "red",lwd = 3, 
								  ylab = "", 
								  xlab = "", 
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)


# add uncertaines
polygon(c(dat.df0001$Year, rev(dat.df0001$Year)), 
		c(dat.df0001$EnsemBorealNoP+dat.df0001$EnsemBorealNoPsd, rev(dat.df0001$EnsemBorealNoP-dat.df0001$EnsemBorealNoPsd)), 
        col=rgb(1, 0, 0,0.25),
		border = NA)
	


## add statiscis 


lmf = function(x,y){

lm1 = lm(y~x)

newx = 1980:2017
conf_interval <- predict(lm1, newdata=data.frame(Year=newx), interval="confidence",level = 0.95)

modsum = summary(lm1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
sen = modsum$coefficients[2,1]
rp = vector('expression',1)
sd1 = modsum$coefficients[2,2]

# rp[1] = substitute(expression(""~ delta == MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

# rp[1] = substitute(expression("Trend:"~ MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

rp[1] = substitute(expression(""~ MYVALUE1 ~ "±" ~ MYVALUE2 ~ gC ~ m^{-2} ~ yr ^{-2}~""), 
		list(MYVALUE1 = format(sen, digits = 2), MYVALUE2 = format(sd1, digits = 2)))[2]

return(list(newx, conf_interval,lm1,rp))
}
				   
x = dat.df0001$Year
y = dat.df0001$EnsemTundraP
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="blue", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="blue", lty=2)
abline(lm1[[3]], lwd = 3,col = "blue")
legend('bottomleft', legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "blue")
# text(x = 2000, y = 0.2*365, labels = expression("Permafrost: 1.00 ± 0.10"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "blue")


x = dat.df0001$Year
y = dat.df0001$EnsemBorealP
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="green", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="green", lty=2)
abline(lm1[[3]], lwd = 3,col = "green")
legend('topleft', legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "green")
# text(x = 2000, y = 0.22*365, labels = expression("Non-Permafrost: 0.22 ± 0.13"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "red")

				   
x = dat.df0001$Year
y = dat.df0001$EnsemBorealNoP
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="red", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="red", lty=2)
abline(lm1[[3]], lwd = 3,col = "red")
legend('topright', legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "red")
# text(x = 2000, y = 0.22*365, labels = expression("Non-Permafrost: 0.22 ± 0.13"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "red")

# inset a subset to seasonsity
d1.df010 = d1.df01

par(new = TRUE)
par(fig = c(0.55, 1, 0.1, 0.5))
v1 = c()
v2 = c()
for(i in 1:12){

lm1 = lm(d1.df010[,i]~d1.df010[,13])
modsum = summary(lm1)
sen = modsum$coefficients[2,1]
v1 = c(v1, sen)

sd1 = modsum$coefficients[2,2]
v2 = c(v2, sd1)

}

v1.df = rbind(v1[c(1,4,7)],v1[c(2,5,8)],v1[c(3,6,9)])

rownames(v1.df) = c("BorealNonP", "BorealP", "TundraP")
colnames(v1.df) = c("EGS", "LGS", "Win")

# plot with sd
v2.df = rbind(v2[c(1,4,7)],v2[c(2,5,8)],v2[c(3,6,9)])

rownames(v2.df) = c("BorealNonP", "BorealP", "TundraP")
colnames(v2.df) = c("EGS", "LGS", "Win")

# us ggplot 2
# http://www.sthda.com/english/wiki/ggplot2-error-NDVIrs-quick-start-guide-r-software-and-data-visualization

v3.df = data.frame(mn = c(v1.df[,1], v1.df[,2], v1.df[,3]), 
                   sd1 = c(v2.df[,1], v2.df[,2], v2.df[,3]),
				   Season = c(rep("EGS",3),rep("LGS",3),rep("Win",3)),
				   Region = rep(c("BorealNonP", "BorealP", "TundraP"),3))
				   

v3.df$Season <- factor(v3.df$Season, levels = c("EGS", "LGS", "Win"))

require(ggplot2)
p <- ggplot(v3.df, aes(x=Season, y=mn, fill=Region)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mn-sd1, ymax=mn+sd1), width=.2,
                 position=position_dodge(.9)) + 
				     theme_bw() + 
  theme(legend.text = element_text(size = 18))+
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
 xlab("") + 
  ylab("") + 
 # ylab(expression(""~ gC ~ m^{-2} ~ yr ^{-2}~"")) + 
    theme(strip.text.x = element_text(size=18), strip.text.y = element_text(size=18)) + 
  theme(legend.position = "none") + 
  scale_fill_manual(values=c("red","green","blue"))

# print(p)


p = p +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

# overlay ggplot2 figure
require(grid)
par(new=TRUE)

print(p, vp=viewport(.78, .25, .45, .3)) #viewport, first two is the x,y coor; last 2 is the inset size

dev.off()

########################################
## different vegetation+P combination

## red for NoP
#11-NoP and TC > 50%
#12-NoP and TC [30-35]
#13-NoP and TC < 30 

# brown for DisconP
#21-DisconP and TC > 50%
#22-DisconP and TC [30-35]
#23-DisconP and TC < 30 

# blue for P
#31-P and TC > 50%
#32-P and TC [30-35]
#33-P and TC < 30 


frost.grd1d2 = mask.1d
frost.grd1d2[] = NA 
frost.grd1d2[frost.grd00 == 11] = 1 #Tree In NoP
frost.grd1d2[frost.grd00 == 12] = 2 #TundraTreeMix In NoP
frost.grd1d2[frost.grd00 == 13] = 3 #Tundra In NoP

frost.grd1d2[frost.grd00 == 21] = 4 #Tree In DisconPP
frost.grd1d2[frost.grd00 == 22] = 5 #TundraTreeMix In DisconPP
frost.grd1d2[frost.grd00 == 23] = 6 #Tundra In DisconPP

frost.grd1d2[frost.grd00 == 31] = 7 #Tree In P
frost.grd1d2[frost.grd00 == 32] = 8 #TundraTreeMix In P
frost.grd1d2[frost.grd00 == 33] = 9 #Tundra In P

d1.df01 = data.frame(cbind(t(data.frame(zonal(SpringC.grd, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
                   # t(data.frame(zonal(SummerC.grd, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
				   t(data.frame(zonal(FallC.grd, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
				   t(data.frame(zonal(WinterC.grd, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
				   t(data.frame(zonal(AnnualC.grd, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
                          c(1980:2017)))

names(d1.df01) = c("SpringNoP_TreeH", "SpringNoP_TreeM","SpringNoP_TreeL","SpringDisconP_TreeH", "SpringDisconP_TreeM","SpringDisconP_TreeL","SpringP_TreeH", "SpringP_TreeM", "SpringP_TreeL", 
				   "FallNoP_TreeH", "FallNoP_TreeM","FallNoP_TreeL","FallDisconP_TreeH", "FallDisconP_TreeM","FallDisconP_TreeL","FallP_TreeH", "FallP_TreeM", "FallP_TreeL", 
				   "WinterNoP_TreeH", "WinterNoP_TreeM","WinterNoP_TreeL","WinterDisconP_TreeH", "WinterDisconP_TreeM","WinterDisconP_TreeL","WinterP_TreeH", "WinterP_TreeM", "WinterP_TreeL", 
				   "AnnualNoP_TreeH", "AnnualNoP_TreeM","AnnualNoP_TreeL","AnnualDisconP_TreeH", "AnnualDisconP_TreeM","AnnualDisconP_TreeL","AnnualP_TreeH", "AnnualP_TreeM", "AnnualP_TreeL", 
				   "Year")
		

ct.df = data.frame(t(zonal(ct.nee*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(ct.df) <- c("AnnualNoP_TreeH", "AnnualNoP_TreeM", "AnnualNoP_TreeL","AnnualDisconP_TreeH", "AnnualDisconP_TreeM", "AnnualDisconP_TreeL","AnnualP_TreeH", "AnnualP_TreeM", "AnnualP_TreeL")
ct.df$Year = as.numeric(as.character(substr(names(ct.nee), 2,5)))
ct.df = aggregate(cbind(AnnualNoP_TreeH, AnnualNoP_TreeM, AnnualNoP_TreeL,AnnualDisconP_TreeH, AnnualDisconP_TreeM, AnnualDisconP_TreeL,AnnualP_TreeH, AnnualP_TreeM, AnnualP_TreeL) ~ Year, FUN=mean, na.rm = TRUE,data = ct.df)

cte.df = data.frame(t(zonal(cte.nee*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(cte.df) <- c("AnnualNoP_TreeH", "AnnualNoP_TreeM", "AnnualNoP_TreeL","AnnualDisconP_TreeH", "AnnualDisconP_TreeM", "AnnualDisconP_TreeL","AnnualP_TreeH", "AnnualP_TreeM", "AnnualP_TreeL")
cte.df$Year = as.numeric(as.character(substr(names(cte.nee), 2,5)))
cte.df = aggregate(cbind(AnnualNoP_TreeH, AnnualNoP_TreeM, AnnualNoP_TreeL,AnnualDisconP_TreeH, AnnualDisconP_TreeM, AnnualDisconP_TreeL,AnnualP_TreeH, AnnualP_TreeM, AnnualP_TreeL) ~ Year, FUN=mean, na.rm = TRUE,data = cte.df)

jena2.df = data.frame(t(zonal(jena.nee2*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(jena2.df) <- c("AnnualNoP_TreeH", "AnnualNoP_TreeM", "AnnualNoP_TreeL","AnnualDisconP_TreeH", "AnnualDisconP_TreeM", "AnnualDisconP_TreeL","AnnualP_TreeH", "AnnualP_TreeM", "AnnualP_TreeL")
jena2.df$Year = as.numeric(as.character(substr(names(jena.nee2), 2,5)))
jena2.df = aggregate(cbind(AnnualNoP_TreeH, AnnualNoP_TreeM, AnnualNoP_TreeL,AnnualDisconP_TreeH, AnnualDisconP_TreeM, AnnualDisconP_TreeL,AnnualP_TreeH, AnnualP_TreeM, AnnualP_TreeL) ~ Year, FUN=mean, na.rm = TRUE,data = jena2.df)

jena.df = data.frame(t(zonal(jena.nee*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(jena.df) <- c("AnnualNoP_TreeH", "AnnualNoP_TreeM", "AnnualNoP_TreeL","AnnualDisconP_TreeH", "AnnualDisconP_TreeM", "AnnualDisconP_TreeL","AnnualP_TreeH", "AnnualP_TreeM", "AnnualP_TreeL")
jena.df$Year = as.numeric(as.character(substr(names(jena.nee), 2,5)))
jena.df = aggregate(cbind(AnnualNoP_TreeH, AnnualNoP_TreeM, AnnualNoP_TreeL,AnnualDisconP_TreeH, AnnualDisconP_TreeM, AnnualDisconP_TreeL,AnnualP_TreeH, AnnualP_TreeM, AnnualP_TreeL) ~ Year, FUN=mean, na.rm = TRUE,data = jena.df)

jams.df = data.frame(t(zonal(jams.nee*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(jams.df) <- c("AnnualNoP_TreeH", "AnnualNoP_TreeM", "AnnualNoP_TreeL","AnnualDisconP_TreeH", "AnnualDisconP_TreeM", "AnnualDisconP_TreeL","AnnualP_TreeH", "AnnualP_TreeM", "AnnualP_TreeL")
jams.df$Year = as.numeric(as.character(substr(names(jams.nee), 2,5)))
jams.df = aggregate(cbind(AnnualNoP_TreeH, AnnualNoP_TreeM, AnnualNoP_TreeL,AnnualDisconP_TreeH, AnnualDisconP_TreeM, AnnualDisconP_TreeL,AnnualP_TreeH, AnnualP_TreeM, AnnualP_TreeL) ~ Year, FUN=mean, na.rm = TRUE,data = jams.df)

cams.df = data.frame(t(zonal(cams.nee*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(cams.df) <- c("AnnualNoP_TreeH", "AnnualNoP_TreeM", "AnnualNoP_TreeL","AnnualDisconP_TreeH", "AnnualDisconP_TreeM", "AnnualDisconP_TreeL","AnnualP_TreeH", "AnnualP_TreeM", "AnnualP_TreeL")
cams.df$Year = as.numeric(as.character(substr(names(cams.nee), 2,5)))
cams.df = aggregate(cbind(AnnualNoP_TreeH, AnnualNoP_TreeM, AnnualNoP_TreeL,AnnualDisconP_TreeH, AnnualDisconP_TreeM, AnnualDisconP_TreeL,AnnualP_TreeH, AnnualP_TreeM, AnnualP_TreeL) ~ Year, FUN=mean, na.rm = TRUE,data = cams.df)

ensemble.df = data.frame(t(zonal(ensemble.nee*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(ensemble.df) <- c("AnnualNoP_TreeH", "AnnualNoP_TreeM", "AnnualNoP_TreeL","AnnualDisconP_TreeH", "AnnualDisconP_TreeM", "AnnualDisconP_TreeL","AnnualP_TreeH", "AnnualP_TreeM", "AnnualP_TreeL")
ensemble.df$Year = as.numeric(as.character(substr(names(ensemble.nee), 2,5)))
ensemble.df = aggregate(cbind(AnnualNoP_TreeH, AnnualNoP_TreeM, AnnualNoP_TreeL,AnnualDisconP_TreeH, AnnualDisconP_TreeM, AnnualDisconP_TreeL,AnnualP_TreeH, AnnualP_TreeM, AnnualP_TreeL) ~ Year, FUN=mean, na.rm = TRUE,data = ensemble.df)


# plot permafrost and non-permrost
dummy.df = data.frame(dummy = NA, Year = 1980:2017)
dat.df0001 = merge(ct.df, dummy.df, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:10)] <- c("Year", "ctAnnualNoP_TreeH", "ctAnnualNoP_TreeM", "ctAnnualNoP_TreeL","ctAnnualDisconP_TreeH", "ctAnnualDisconP_TreeM", "ctAnnualDisconP_TreeL","ctAnnualP_TreeH", "ctAnnualP_TreeM", "ctAnnualP_TreeL")

dat.df0001 = merge(cte.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:10)] <- c("Year", "cteAnnualNoP_TreeH", "cteAnnualNoP_TreeM", "cteAnnualNoP_TreeL","cteAnnualDisconP_TreeH", "cteAnnualDisconP_TreeM", "cteAnnualDisconP_TreeL","cteAnnualP_TreeH", "cteAnnualP_TreeM", "cteAnnualP_TreeL")

dat.df0001 = merge(jena2.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:10)] <- c("Year", "jena2AnnualNoP_TreeH", "jena2AnnualNoP_TreeM", "jena2AnnualNoP_TreeL","jena2AnnualDisconP_TreeH", "jena2AnnualDisconP_TreeM", "jena2AnnualDisconP_TreeL","jena2AnnualP_TreeH", "jena2AnnualP_TreeM", "jena2AnnualP_TreeL")

dat.df0001 = merge(jena.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:10)] <- c("Year", "jenaAnnualNoP_TreeH", "jenaAnnualNoP_TreeM", "jenaAnnualNoP_TreeL","jenaAnnualDisconP_TreeH", "jenaAnnualDisconP_TreeM", "jenaAnnualDisconP_TreeL","jenaAnnualP_TreeH", "jenaAnnualP_TreeM", "jenaAnnualP_TreeL")

dat.df0001 = merge(jams.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:10)] <- c("Year", "jamsAnnualNoP_TreeH", "jamsAnnualNoP_TreeM", "jamsAnnualNoP_TreeL","jamsAnnualDisconP_TreeH", "jamsAnnualDisconP_TreeM", "jamsAnnualDisconP_TreeL","jamsAnnualP_TreeH", "jamsAnnualP_TreeM", "jamsAnnualP_TreeL")

dat.df0001 = merge(cams.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:10)] <- c("Year", "camsAnnualNoP_TreeH", "camsAnnualNoP_TreeM", "camsAnnualNoP_TreeL","camsAnnualDisconP_TreeH", "camsAnnualDisconP_TreeM", "camsAnnualDisconP_TreeL","camsAnnualP_TreeH", "camsAnnualP_TreeM", "camsAnnualP_TreeL")

dat.df0001 = merge(ensemble.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:10)] <- c("Year", "EnsemAnnualNoP_TreeH", "EnsemAnnualNoP_TreeM", "EnsemAnnualNoP_TreeL","EnsemAnnualDisconP_TreeH", "EnsemAnnualDisconP_TreeM", "EnsemAnnualDisconP_TreeL","EnsemAnnualP_TreeH", "EnsemAnnualP_TreeM", "EnsemAnnualP_TreeL")

sde = function(x){
x1 = sd(x, na.rm = TRUE)
n = length(which(!is.na(x)))
return(x1/(sqrt(n)))
}

dat.df0001$EnsemAnnualNoP_TreeHsd = apply(dat.df0001[,c("ctAnnualNoP_TreeH", "cteAnnualNoP_TreeH", "jena2AnnualNoP_TreeH","jenaAnnualNoP_TreeH", "jamsAnnualNoP_TreeH", "camsAnnualNoP_TreeH")], 1, sd, na.rm = TRUE)
dat.df0001$EnsemAnnualNoP_TreeMsd = apply(dat.df0001[,c("ctAnnualNoP_TreeM", "cteAnnualNoP_TreeM", "jena2AnnualNoP_TreeM","jenaAnnualNoP_TreeM", "jamsAnnualNoP_TreeM", "camsAnnualNoP_TreeM")], 1, sd, na.rm = TRUE)
dat.df0001$EnsemAnnualNoP_TreeLsd = apply(dat.df0001[,c("ctAnnualNoP_TreeL", "cteAnnualNoP_TreeL", "jena2AnnualNoP_TreeL","jenaAnnualNoP_TreeL", "jamsAnnualNoP_TreeL", "camsAnnualNoP_TreeL")], 1, sd, na.rm = TRUE)
dat.df0001$EnsemAnnualDisconP_TreeHsd = apply(dat.df0001[,c("ctAnnualDisconP_TreeH", "cteAnnualDisconP_TreeH", "jena2AnnualDisconP_TreeH","jenaAnnualDisconP_TreeH", "jamsAnnualDisconP_TreeH", "camsAnnualDisconP_TreeH")], 1, sd, na.rm = TRUE)
dat.df0001$EnsemAnnualDisconP_TreeMsd = apply(dat.df0001[,c("ctAnnualDisconP_TreeM", "cteAnnualDisconP_TreeM", "jena2AnnualDisconP_TreeM","jenaAnnualDisconP_TreeM", "jamsAnnualDisconP_TreeM", "camsAnnualDisconP_TreeM")], 1, sd, na.rm = TRUE)
dat.df0001$EnsemAnnualDisconP_TreeLsd = apply(dat.df0001[,c("ctAnnualDisconP_TreeL", "cteAnnualDisconP_TreeL", "jena2AnnualDisconP_TreeL","jenaAnnualDisconP_TreeL", "jamsAnnualDisconP_TreeL", "camsAnnualDisconP_TreeL")], 1, sd, na.rm = TRUE)
dat.df0001$EnsemAnnualP_TreeHsd = apply(dat.df0001[,c("ctAnnualP_TreeH", "cteAnnualP_TreeH", "jena2AnnualP_TreeH","jenaAnnualP_TreeH", "jamsAnnualP_TreeH", "camsAnnualP_TreeH")], 1, sd, na.rm = TRUE)
dat.df0001$EnsemAnnualP_TreeMsd = apply(dat.df0001[,c("ctAnnualP_TreeM", "cteAnnualP_TreeM", "jena2AnnualP_TreeM","jenaAnnualP_TreeM", "jamsAnnualP_TreeM", "camsAnnualP_TreeM")], 1, sd, na.rm = TRUE)
dat.df0001$EnsemAnnualP_TreeLsd = apply(dat.df0001[,c("ctAnnualP_TreeL", "cteAnnualP_TreeL", "jena2AnnualP_TreeL","jenaAnnualP_TreeL", "jamsAnnualP_TreeL", "camsAnnualP_TreeL")], 1, sd, na.rm = TRUE)


## plot 
png("R:/zhihua/L4C_NH_analysis_Revision3/Figures/CsinkTrend_3_regions_Permafrost+TC3.png",height = 3600, width = 3600, res = 300, units = "px", bg = "transparent")

par(mar=c(4,5,5,0)+.1)

# plot Permafrost region [deep blue = tree #08306b, mid blue = mix #2171b5, light blue = shrub #9ecae1]
plot(EnsemAnnualP_TreeH~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017, main = "",
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "#08306b",lwd = 3, 
								  xlab = "Year", ylab = expression("Net C Uptake (" ~ gC ~ m^{-2} ~ yr ^{-1}~ ")"),
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)


lmf = function(x,y){

lm1 = lm(y~x)

newx = 1980:2017
conf_interval <- predict(lm1, newdata=data.frame(Year=newx), interval="confidence",level = 0.95)

modsum = summary(lm1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
sen = modsum$coefficients[2,1]
rp = vector('expression',1)
sd1 = modsum$coefficients[2,2]

# rp[1] = substitute(expression(""~ delta == MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

# rp[1] = substitute(expression("Trend:"~ MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

rp[1] = substitute(expression("TC = [50,100) in P: "~ MYVALUE1 ~ "±" ~ MYVALUE2 ~ gC ~ m^{-2} ~ yr ^{-2}~""), 
		list(MYVALUE1 = format(sen, digits = 2), MYVALUE2 = format(sd1, digits = 2)))[2]

return(list(newx, conf_interval,lm1,rp))
}

x = dat.df0001$Year
y = dat.df0001$EnsemAnnualP_TreeH
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="#08306b", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="#08306b", lty=2)
abline(lm1[[3]], lwd = 3,col = "#08306b")
legend(x = 1980, y = 10+80, legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "#08306b")
# text(x = 2000, y = 0.2*365, labels = expression("Permafrost: 1.00 ± 0.10"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "blue")

par(new=TRUE)

plot(EnsemAnnualP_TreeM~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017,
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "#2171b5",lwd = 3, 
								  xlab = "Year", ylab = expression("Net C Uptake (" ~ gC ~ m^{-2} ~ yr ^{-1}~ ")"),
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)


lmf = function(x,y){

lm1 = lm(y~x)

newx = 1980:2017
conf_interval <- predict(lm1, newdata=data.frame(Year=newx), interval="confidence",level = 0.95)

modsum = summary(lm1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
sen = modsum$coefficients[2,1]
rp = vector('expression',1)
sd1 = modsum$coefficients[2,2]

# rp[1] = substitute(expression(""~ delta == MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

# rp[1] = substitute(expression("Trend:"~ MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

rp[1] = substitute(expression("TC = [30,50) in P: "~ MYVALUE1 ~ "±" ~ MYVALUE2 ~ gC ~ m^{-2} ~ yr ^{-2}~""), 
		list(MYVALUE1 = format(sen, digits = 2), MYVALUE2 = format(sd1, digits = 2)))[2]

return(list(newx, conf_interval,lm1,rp))
}

x = dat.df0001$Year
y = dat.df0001$EnsemAnnualP_TreeM
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="#2171b5", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="#2171b5", lty=2)
abline(lm1[[3]], lwd = 3,col = "#2171b5")
legend(x = 1980, y = 10+75, legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "#2171b5")
# text(x = 2000, y = 0.2*365, labels = expression("Permafrost: 1.00 ± 0.10"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "blue")

par(new=TRUE)

plot(EnsemAnnualP_TreeL~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017,
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "#9ecae1",lwd = 3, 
								  xlab = "Year", ylab = expression("Net C Uptake (" ~ gC ~ m^{-2} ~ yr ^{-1}~ ")"),
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)


lmf = function(x,y){

lm1 = lm(y~x)

newx = 1980:2017
conf_interval <- predict(lm1, newdata=data.frame(Year=newx), interval="confidence",level = 0.95)

modsum = summary(lm1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
sen = modsum$coefficients[2,1]
rp = vector('expression',1)
sd1 = modsum$coefficients[2,2]

# rp[1] = substitute(expression(""~ delta == MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

# rp[1] = substitute(expression("Trend:"~ MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

rp[1] = substitute(expression("TC = [0,30] in P: "~ MYVALUE1 ~ "±" ~ MYVALUE2 ~ gC ~ m^{-2} ~ yr ^{-2}~""), 
		list(MYVALUE1 = format(sen, digits = 2), MYVALUE2 = format(sd1, digits = 2)))[2]

return(list(newx, conf_interval,lm1,rp))
}

x = dat.df0001$Year
y = dat.df0001$EnsemAnnualP_TreeL
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="#9ecae1", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="#9ecae1", lty=2)
abline(lm1[[3]], lwd = 3,col = "#9ecae1")
legend(x = 1980, y = 10+70, legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "#9ecae1")
# text(x = 2000, y = 0.2*365, labels = expression("Permafrost: 1.00 ± 0.10"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "blue")
								  
# plot discontinuous Permafrost region [deep green = tree #00441b, mid blue = mix #41ab5d light  = shrub #c7e9c0]
par(new=TRUE)

plot(EnsemAnnualDisconP_TreeH~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017, main = "",
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "#00441b",lwd = 3, 
								  xlab = "Year", ylab = expression("Net C Uptake (" ~ gC ~ m^{-2} ~ yr ^{-1}~ ")"),
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)


lmf = function(x,y){

lm1 = lm(y~x)

newx = 1980:2017
conf_interval <- predict(lm1, newdata=data.frame(Year=newx), interval="confidence",level = 0.95)

modsum = summary(lm1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
sen = modsum$coefficients[2,1]
rp = vector('expression',1)
sd1 = modsum$coefficients[2,2]

# rp[1] = substitute(expression(""~ delta == MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

# rp[1] = substitute(expression("Trend:"~ MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

rp[1] = substitute(expression("TC = [50,100) in DisconP: "~ MYVALUE1 ~ "±" ~ MYVALUE2 ~ gC ~ m^{-2} ~ yr ^{-2}~""), 
		list(MYVALUE1 = format(sen, digits = 2), MYVALUE2 = format(sd1, digits = 2)))[2]

return(list(newx, conf_interval,lm1,rp))
}

x = dat.df0001$Year
y = dat.df0001$EnsemAnnualDisconP_TreeH
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="#00441b", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="#00441b", lty=2)
abline(lm1[[3]], lwd = 3,col = "#00441b")
legend(x = 1980, y = 10+65, legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "#00441b")
# text(x = 2000, y = 0.2*365, labels = expression("Permafrost: 1.00 ± 0.10"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "blue")

par(new=TRUE)

plot(EnsemAnnualDisconP_TreeM~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017,
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "#41ab5d",lwd = 3, 
								  xlab = "Year", ylab = expression("Net C Uptake (" ~ gC ~ m^{-2} ~ yr ^{-1}~ ")"),
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)

lmf = function(x,y){

lm1 = lm(y~x)

newx = 1980:2017
conf_interval <- predict(lm1, newdata=data.frame(Year=newx), interval="confidence",level = 0.95)

modsum = summary(lm1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
sen = modsum$coefficients[2,1]
rp = vector('expression',1)
sd1 = modsum$coefficients[2,2]

# rp[1] = substitute(expression(""~ delta == MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

# rp[1] = substitute(expression("Trend:"~ MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

rp[1] = substitute(expression("TC = [30,50) in DisconP: "~ MYVALUE1 ~ "±" ~ MYVALUE2 ~ gC ~ m^{-2} ~ yr ^{-2}~""), 
		list(MYVALUE1 = format(sen, digits = 2), MYVALUE2 = format(sd1, digits = 2)))[2]

return(list(newx, conf_interval,lm1,rp))
}
								  
x = dat.df0001$Year
y = dat.df0001$EnsemAnnualDisconP_TreeM
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="#41ab5d", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="#41ab5d", lty=2)
abline(lm1[[3]], lwd = 3,col = "#41ab5d")
legend(x = 1980, y = 10+60, legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "#41ab5d")
# text(x = 2000, y = 0.2*365, labels = expression("Permafrost: 1.00 ± 0.10"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "blue")

par(new=TRUE)

plot(EnsemAnnualDisconP_TreeL~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017,
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "#c7e9c0",lwd = 3, 
								  xlab = "Year", ylab = expression("Net C Uptake (" ~ gC ~ m^{-2} ~ yr ^{-1}~ ")"),
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)

lmf = function(x,y){

lm1 = lm(y~x)

newx = 1980:2017
conf_interval <- predict(lm1, newdata=data.frame(Year=newx), interval="confidence",level = 0.95)

modsum = summary(lm1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
sen = modsum$coefficients[2,1]
rp = vector('expression',1)
sd1 = modsum$coefficients[2,2]

# rp[1] = substitute(expression(""~ delta == MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

# rp[1] = substitute(expression("Trend:"~ MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

rp[1] = substitute(expression("TC = [0,30) in DisconP: "~ MYVALUE1 ~ "±" ~ MYVALUE2 ~ gC ~ m^{-2} ~ yr ^{-2}~""), 
		list(MYVALUE1 = format(sen, digits = 2), MYVALUE2 = format(sd1, digits = 2)))[2]

return(list(newx, conf_interval,lm1,rp))
}
								  
x = dat.df0001$Year
y = dat.df0001$EnsemAnnualDisconP_TreeL
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="#c7e9c0", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="#c7e9c0", lty=2)
abline(lm1[[3]], lwd = 3,col = "#c7e9c0")
legend(x = 1980, y = 10+55, legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "#c7e9c0")
# text(x = 2000, y = 0.2*365, labels = expression("Permafrost: 1.00 ± 0.10"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "blue")
	
# plot no Permafrost region [deep red = tree #67000d, mid  = mix #ef3b2c light  = shrub #fcbba1]
par(new=TRUE)

plot(EnsemAnnualNoP_TreeH~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017, main = "",
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "#67000d",lwd = 3, 
								  xlab = "Year", ylab = expression("Net C Uptake (" ~ gC ~ m^{-2} ~ yr ^{-1}~ ")"),
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)

lmf = function(x,y){

lm1 = lm(y~x)

newx = 1980:2017
conf_interval <- predict(lm1, newdata=data.frame(Year=newx), interval="confidence",level = 0.95)

modsum = summary(lm1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
sen = modsum$coefficients[2,1]
rp = vector('expression',1)
sd1 = modsum$coefficients[2,2]

# rp[1] = substitute(expression(""~ delta == MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

# rp[1] = substitute(expression("Trend:"~ MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

rp[1] = substitute(expression("TC = (50,100) in NoP: "~ MYVALUE1 ~ "±" ~ MYVALUE2 ~ gC ~ m^{-2} ~ yr ^{-2}~""), 
		list(MYVALUE1 = format(sen, digits = 2), MYVALUE2 = format(sd1, digits = 2)))[2]

return(list(newx, conf_interval,lm1,rp))
}

x = dat.df0001$Year
y = dat.df0001$EnsemAnnualNoP_TreeH
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="#67000d", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="#67000d", lty=2)
abline(lm1[[3]], lwd = 3,col = "#67000d")
legend(x = 1980, y = 10+50, legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "#67000d")
# text(x = 2000, y = 0.2*365, labels = expression("Permafrost: 1.00 ± 0.10"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "blue")
	
par(new=TRUE)

plot(EnsemAnnualNoP_TreeM~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017,
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "#ef3b2c",lwd = 3, 
								  xlab = "Year", ylab = expression("Net C Uptake (" ~ gC ~ m^{-2} ~ yr ^{-1}~ ")"),
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)

lmf = function(x,y){

lm1 = lm(y~x)

newx = 1980:2017
conf_interval <- predict(lm1, newdata=data.frame(Year=newx), interval="confidence",level = 0.95)

modsum = summary(lm1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
sen = modsum$coefficients[2,1]
rp = vector('expression',1)
sd1 = modsum$coefficients[2,2]

# rp[1] = substitute(expression(""~ delta == MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

# rp[1] = substitute(expression("Trend:"~ MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

rp[1] = substitute(expression("TC = (30,50) in NoP: "~ MYVALUE1 ~ "±" ~ MYVALUE2 ~ gC ~ m^{-2} ~ yr ^{-2}~""), 
		list(MYVALUE1 = format(sen, digits = 2), MYVALUE2 = format(sd1, digits = 2)))[2]

return(list(newx, conf_interval,lm1,rp))
}

x = dat.df0001$Year
y = dat.df0001$EnsemAnnualNoP_TreeM
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="#ef3b2c", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="#ef3b2c", lty=2)
abline(lm1[[3]], lwd = 3,col = "#ef3b2c")
legend(x = 1980, y = 10+45, legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "#ef3b2c")
# text(x = 2000, y = 0.2*365, labels = expression("Permafrost: 1.00 ± 0.10"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "blue")
	
par(new=TRUE)

plot(EnsemAnnualNoP_TreeL~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017,
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "#fcbba1",lwd = 3, 
								  xlab = "Year", ylab = expression("Net C Uptake (" ~ gC ~ m^{-2} ~ yr ^{-1}~ ")"),
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)
	
lmf = function(x,y){

lm1 = lm(y~x)

newx = 1980:2017
conf_interval <- predict(lm1, newdata=data.frame(Year=newx), interval="confidence",level = 0.95)

modsum = summary(lm1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
sen = modsum$coefficients[2,1]
rp = vector('expression',1)
sd1 = modsum$coefficients[2,2]

# rp[1] = substitute(expression(""~ delta == MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

# rp[1] = substitute(expression("Trend:"~ MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

rp[1] = substitute(expression("TC = (0,30) in NoP: "~ MYVALUE1 ~ "±" ~ MYVALUE2 ~ gC ~ m^{-2} ~ yr ^{-2}~""), 
		list(MYVALUE1 = format(sen, digits = 2), MYVALUE2 = format(sd1, digits = 2)))[2]

return(list(newx, conf_interval,lm1,rp))
}
	
x = dat.df0001$Year
y = dat.df0001$EnsemAnnualNoP_TreeL
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="#fcbba1", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="#fcbba1", lty=2)
abline(lm1[[3]], lwd = 3,col = "#fcbba1")
legend(x = 1980, y = 10+40, legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "#fcbba1")
# text(x = 2000, y = 0.2*365, labels = expression("Permafrost: 1.00 ± 0.10"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "blue")

## add statiscis 

	
# inset a subset to seasonsity
d1.df010 = d1.df01

par(new = TRUE)
par(fig = c(0.55, 1, 0.1, 0.5))
v1 = c()
v2 = c()
for(i in 1:36){

lm1 = lm(d1.df010[,i]~d1.df010[,37])
modsum = summary(lm1)
sen = modsum$coefficients[2,1]
v1 = c(v1, sen)

sd1 = modsum$coefficients[2,2]
v2 = c(v2, sd1)

}

v1.df = t(rbind(v1[c(1:9)],v1[c(10:18)],v1[c(19:27)]))

rownames(v1.df) = c("NoP_TreeH", "NoP_TreeM", "NoP_TreeL","DisconP_TreeH", "DisconP_TreeM", "DisconP_TreeL","P_TreeH", "P_TreeM", "P_TreeL")
colnames(v1.df) = c("EGS", "LGS", "Win")

# plot with sd
v2.df = t(rbind(v2[c(1:9)],v2[c(10:18)],v2[c(19:27)]))

rownames(v2.df) = c("NoP_TreeH", "NoP_TreeM", "NoP_TreeL","DisconP_TreeH", "DisconP_TreeM", "DisconP_TreeL","P_TreeH", "P_TreeM", "P_TreeL")
colnames(v2.df) = c("EGS", "LGS", "Win")

# us ggplot 2
# http://www.sthda.com/english/wiki/ggplot2-error-NDVIrs-quick-start-guide-r-software-and-data-visualization

v3.df = data.frame(mn = c(v1.df[,1], v1.df[,2], v1.df[,3]), 
                   sd1 = c(v2.df[,1], v2.df[,2], v2.df[,3]),
				   Season = c(rep("EGS",9),rep("LGS",9),rep("Win",9)),
				   Region = rep(c("NoP_TreeH", "NoP_TreeM", "NoP_TreeL","DisconP_TreeH", "DisconP_TreeM", "DisconP_TreeL","P_TreeH", "P_TreeM", "P_TreeL"),3))
				   

v3.df$Season <- factor(v3.df$Season, levels = c("EGS", "LGS", "Win"))
v3.df$Region <- factor(v3.df$Region, levels = c("NoP_TreeH", "NoP_TreeM", "NoP_TreeL","DisconP_TreeH", "DisconP_TreeM", "DisconP_TreeL","P_TreeH", "P_TreeM", "P_TreeL"))

require(ggplot2)
p <- ggplot(v3.df, aes(x=Season, y=mn, fill=Region)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mn-sd1, ymax=mn+sd1), width=.2,
                 position=position_dodge(.9)) + 
				     theme_bw() + 
  theme(legend.text = element_text(size = 18))+
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
 xlab("") + 
  ylab("") + 
 # ylab(expression(""~ gC ~ m^{-2} ~ yr ^{-2}~"")) + 
    theme(strip.text.x = element_text(size=18), strip.text.y = element_text(size=18)) + 
  theme(legend.position = "none") + 
  scale_fill_manual(values=c("#67000d","#ef3b2c","#fcbba1","#00441b","#41ab5d","#c7e9c0", "#08306b","#4292c6","#c6dbef"))

# print(p)


p = p +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

# overlay ggplot2 figure
require(grid)
par(new=TRUE)

print(p, vp=viewport(.78, .25, .45, .3)) #viewport, first two is the x,y coor; last 2 is the inset size

dev.off()




########################################
## different vegetation+P combination
## P = NoP vs P
## TC = TC < 50 vs TC > 50

## red for NoP
#11-NoP and TC > 50%::: 1
#12-NoP and TC [30-35] :::2
#13-NoP and TC < 30 :::2 

# brown for DisconP
#21-DisconP and TC > 50% :::3
#22-DisconP and TC [30-35] :::4
#23-DisconP and TC < 30 :::4

# blue for P
#31-P and TC > 50% :::3
#32-P and TC [30-35] :::4
#33-P and TC < 30 :::4

1 == HighTC + NoP
2 == LowTC + NoP
3 == HighTC + P
4 == LowTC + P


frost.grd1d2 = mask.1d
frost.grd1d2[] = NA 
frost.grd1d2[frost.grd00 == 11] = 1 
frost.grd1d2[frost.grd00 == 12] = 2 
frost.grd1d2[frost.grd00 == 13] = 2 

frost.grd1d2[frost.grd00 == 21] = 3 
frost.grd1d2[frost.grd00 == 22] = 4 
frost.grd1d2[frost.grd00 == 23] = 4 

frost.grd1d2[frost.grd00 == 31] = 3 
frost.grd1d2[frost.grd00 == 32] = 4 
frost.grd1d2[frost.grd00 == 33] = 4 

table(as.vector(as.matrix(frost.grd1d2)))


d1.df01 = data.frame(cbind(t(data.frame(zonal(SpringC.grd, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
                   # t(data.frame(zonal(SummerC.grd, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
				   t(data.frame(zonal(FallC.grd, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
				   t(data.frame(zonal(WinterC.grd, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
				   t(data.frame(zonal(AnnualC.grd, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365, 
                          c(1980:2017)))

names(d1.df01) = c("SpringNoP_TreeH", "SpringNoP_TreeL","SpringP_TreeH", "SpringP_TreeL", 
				   "FallNoP_TreeH", "FallNoP_TreeL","FallP_TreeH", "FallP_TreeL", 
				   "WinterNoP_TreeH", "WinterNoP_TreeL","WinterP_TreeH", "WinterP_TreeL", 
				   "AnnualNoP_TreeH", "AnnualNoP_TreeL","AnnualP_TreeH",  "AnnualP_TreeL", 
				   "Year")
		

ct.df = data.frame(t(zonal(ct.nee*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(ct.df) <- c("AnnualNoP_TreeH", "AnnualNoP_TreeL","AnnualP_TreeH", "AnnualP_TreeL")
ct.df$Year = as.numeric(as.character(substr(names(ct.nee), 2,5)))
ct.df = aggregate(cbind(AnnualNoP_TreeH, AnnualNoP_TreeL,AnnualP_TreeH, AnnualP_TreeL) ~ Year, FUN=mean, na.rm = TRUE,data = ct.df)

cte.df = data.frame(t(zonal(cte.nee*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(cte.df) <- c("AnnualNoP_TreeH", "AnnualNoP_TreeL","AnnualP_TreeH", "AnnualP_TreeL")
cte.df$Year = as.numeric(as.character(substr(names(cte.nee), 2,5)))
cte.df = aggregate(cbind(AnnualNoP_TreeH, AnnualNoP_TreeL,AnnualP_TreeH, AnnualP_TreeL) ~ Year, FUN=mean, na.rm = TRUE,data = cte.df)

jena2.df = data.frame(t(zonal(jena.nee2*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(jena2.df) <- c("AnnualNoP_TreeH", "AnnualNoP_TreeL","AnnualP_TreeH", "AnnualP_TreeL")
jena2.df$Year = as.numeric(as.character(substr(names(jena.nee2), 2,5)))
jena2.df = aggregate(cbind(AnnualNoP_TreeH, AnnualNoP_TreeL,AnnualP_TreeH, AnnualP_TreeL) ~ Year, FUN=mean, na.rm = TRUE,data = jena2.df)

jena.df = data.frame(t(zonal(jena.nee*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(jena.df) <- c("AnnualNoP_TreeH", "AnnualNoP_TreeL","AnnualP_TreeH", "AnnualP_TreeL")
jena.df$Year = as.numeric(as.character(substr(names(jena.nee), 2,5)))
jena.df = aggregate(cbind(AnnualNoP_TreeH, AnnualNoP_TreeL,AnnualP_TreeH, AnnualP_TreeL) ~ Year, FUN=mean, na.rm = TRUE,data = jena.df)

jams.df = data.frame(t(zonal(jams.nee*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(jams.df) <- c("AnnualNoP_TreeH", "AnnualNoP_TreeL","AnnualP_TreeH", "AnnualP_TreeL")
jams.df$Year = as.numeric(as.character(substr(names(jams.nee), 2,5)))
jams.df = aggregate(cbind(AnnualNoP_TreeH, AnnualNoP_TreeL,AnnualP_TreeH, AnnualP_TreeL) ~ Year, FUN=mean, na.rm = TRUE,data = jams.df)

cams.df = data.frame(t(zonal(cams.nee*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(cams.df) <- c("AnnualNoP_TreeH", "AnnualNoP_TreeL","AnnualP_TreeH", "AnnualP_TreeL")
cams.df$Year = as.numeric(as.character(substr(names(cams.nee), 2,5)))
cams.df = aggregate(cbind(AnnualNoP_TreeH, AnnualNoP_TreeL,AnnualP_TreeH, AnnualP_TreeL) ~ Year, FUN=mean, na.rm = TRUE,data = cams.df)

ensemble.df = data.frame(t(zonal(ensemble.nee*mask.1d, frost.grd1d2, fun='mean', digits=0, na.rm=TRUE)))[-1,]*365
names(ensemble.df) <- c("AnnualNoP_TreeH", "AnnualNoP_TreeL","AnnualP_TreeH", "AnnualP_TreeL")
ensemble.df$Year = as.numeric(as.character(substr(names(ensemble.nee), 2,5)))
ensemble.df = aggregate(cbind(AnnualNoP_TreeH, AnnualNoP_TreeL,AnnualP_TreeH, AnnualP_TreeL) ~ Year, FUN=mean, na.rm = TRUE,data = ensemble.df)


# plot permafrost and non-permrost
dummy.df = data.frame(dummy = NA, Year = 1980:2017)
dat.df0001 = merge(ct.df, dummy.df, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:5)] <- c("Year", "ctAnnualNoP_TreeH", "ctAnnualNoP_TreeL","ctAnnualP_TreeH", "ctAnnualP_TreeL")

dat.df0001 = merge(cte.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:5)] <- c("Year", "cteAnnualNoP_TreeH", "cteAnnualNoP_TreeL","cteAnnualP_TreeH", "cteAnnualP_TreeL")

dat.df0001 = merge(jena2.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:5)] <- c("Year", "jena2AnnualNoP_TreeH","jena2AnnualNoP_TreeL","jena2AnnualP_TreeH", "jena2AnnualP_TreeL")

dat.df0001 = merge(jena.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:5)] <- c("Year", "jenaAnnualNoP_TreeH", "jenaAnnualNoP_TreeL","jenaAnnualP_TreeH", "jenaAnnualP_TreeL")

dat.df0001 = merge(jams.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:5)] <- c("Year", "jamsAnnualNoP_TreeH", "jamsAnnualNoP_TreeL","jamsAnnualP_TreeH", "jamsAnnualP_TreeL")

dat.df0001 = merge(cams.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:5)] <- c("Year", "camsAnnualNoP_TreeH", "camsAnnualNoP_TreeL","camsAnnualP_TreeH", "camsAnnualP_TreeL")

dat.df0001 = merge(ensemble.df, dat.df0001, by = "Year", all.y = TRUE)
names(dat.df0001)[c(1:5)] <- c("Year", "EnsemAnnualNoP_TreeH", "EnsemAnnualNoP_TreeL","EnsemAnnualP_TreeH", "EnsemAnnualP_TreeL")

sde = function(x){
x1 = sd(x, na.rm = TRUE)
n = length(which(!is.na(x)))
return(x1/(sqrt(n)))
}

dat.df0001$EnsemAnnualNoP_TreeHsd = apply(dat.df0001[,c("ctAnnualNoP_TreeH", "cteAnnualNoP_TreeH", "jena2AnnualNoP_TreeH","jenaAnnualNoP_TreeH", "jamsAnnualNoP_TreeH", "camsAnnualNoP_TreeH")], 1, sd, na.rm = TRUE)
dat.df0001$EnsemAnnualNoP_TreeLsd = apply(dat.df0001[,c("ctAnnualNoP_TreeL", "cteAnnualNoP_TreeL", "jena2AnnualNoP_TreeL","jenaAnnualNoP_TreeL", "jamsAnnualNoP_TreeL", "camsAnnualNoP_TreeL")], 1, sd, na.rm = TRUE)
dat.df0001$EnsemAnnualP_TreeHsd = apply(dat.df0001[,c("ctAnnualP_TreeH", "cteAnnualP_TreeH", "jena2AnnualP_TreeH","jenaAnnualP_TreeH", "jamsAnnualP_TreeH", "camsAnnualP_TreeH")], 1, sd, na.rm = TRUE)
dat.df0001$EnsemAnnualP_TreeLsd = apply(dat.df0001[,c("ctAnnualP_TreeL", "cteAnnualP_TreeL", "jena2AnnualP_TreeL","jenaAnnualP_TreeL", "jamsAnnualP_TreeL", "camsAnnualP_TreeL")], 1, sd, na.rm = TRUE)


## plot 
png("R:/zhihua/L4C_NH_analysis_Revision3/Figures/CsinkTrend_3_regions_Permafrost+TC4.png",height = 3000, width = 3000, res = 300, units = "px", bg = "transparent")

par(mar=c(4,5,5,0)+.1)

# plot Permafrost region [deep blue = tree #08306b, mid blue = mix #2171b5, light blue = shrub #9ecae1]
plot(EnsemAnnualP_TreeH~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017, main = "",
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "#08306b",lwd = 3, 
								  xlab = "Year", ylab = expression("Net C Uptake (" ~ gC ~ m^{-2} ~ yr ^{-1}~ ")"),
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)


lmf = function(x,y){

lm1 = lm(y~x)

newx = 1980:2017
conf_interval <- predict(lm1, newdata=data.frame(Year=newx), interval="confidence",level = 0.95)

modsum = summary(lm1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
sen = modsum$coefficients[2,1]
rp = vector('expression',1)
sd1 = modsum$coefficients[2,2]

# rp[1] = substitute(expression(""~ delta == MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

# rp[1] = substitute(expression("Trend:"~ MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

rp[1] = substitute(expression("TC = (50,100) in P: "~ MYVALUE1 ~ "±" ~ MYVALUE2 ~ gC ~ m^{-2} ~ yr ^{-2}~""), 
		list(MYVALUE1 = format(sen, digits = 2), MYVALUE2 = format(sd1, digits = 2)))[2]

return(list(newx, conf_interval,lm1,rp))
}

x = dat.df0001$Year
y = dat.df0001$EnsemAnnualP_TreeH
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="#08306b", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="#08306b", lty=2)
abline(lm1[[3]], lwd = 3,col = "#08306b")
legend(x = 1980, y = 10+80, legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "#08306b")
# text(x = 2000, y = 0.2*365, labels = expression("Permafrost: 1.00 ± 0.10"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "blue")

par(new=TRUE)

plot(EnsemAnnualP_TreeL~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017,
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "#9ecae1",lwd = 3, 
								  xlab = "Year", ylab = expression("Net C Uptake (" ~ gC ~ m^{-2} ~ yr ^{-1}~ ")"),
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)


lmf = function(x,y){

lm1 = lm(y~x)

newx = 1980:2017
conf_interval <- predict(lm1, newdata=data.frame(Year=newx), interval="confidence",level = 0.95)

modsum = summary(lm1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
sen = modsum$coefficients[2,1]
rp = vector('expression',1)
sd1 = modsum$coefficients[2,2]

# rp[1] = substitute(expression(""~ delta == MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

# rp[1] = substitute(expression("Trend:"~ MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

rp[1] = substitute(expression("TC = (0,50) in P: "~ MYVALUE1 ~ "±" ~ MYVALUE2 ~ gC ~ m^{-2} ~ yr ^{-2}~""), 
		list(MYVALUE1 = format(sen, digits = 2), MYVALUE2 = format(sd1, digits = 2)))[2]

return(list(newx, conf_interval,lm1,rp))
}

x = dat.df0001$Year
y = dat.df0001$EnsemAnnualP_TreeL
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="#9ecae1", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="#9ecae1", lty=2)
abline(lm1[[3]], lwd = 3,col = "#9ecae1")
legend(x = 1980, y = 10+75, legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "#9ecae1")
# text(x = 2000, y = 0.2*365, labels = expression("Permafrost: 1.00 ± 0.10"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "blue")
			
# plot no Permafrost region [deep red = tree #67000d, mid  = mix #ef3b2c light  = shrub #fcbba1]
par(new=TRUE)

plot(EnsemAnnualNoP_TreeH~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017, main = "",
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "#67000d",lwd = 3, 
								  xlab = "Year", ylab = expression("Net C Uptake (" ~ gC ~ m^{-2} ~ yr ^{-1}~ ")"),
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)

lmf = function(x,y){

lm1 = lm(y~x)

newx = 1980:2017
conf_interval <- predict(lm1, newdata=data.frame(Year=newx), interval="confidence",level = 0.95)

modsum = summary(lm1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
sen = modsum$coefficients[2,1]
rp = vector('expression',1)
sd1 = modsum$coefficients[2,2]

# rp[1] = substitute(expression(""~ delta == MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

# rp[1] = substitute(expression("Trend:"~ MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

rp[1] = substitute(expression("TC = (50,100) in NoP: "~ MYVALUE1 ~ "±" ~ MYVALUE2 ~ gC ~ m^{-2} ~ yr ^{-2}~""), 
		list(MYVALUE1 = format(sen, digits = 2), MYVALUE2 = format(sd1, digits = 2)))[2]

return(list(newx, conf_interval,lm1,rp))
}

x = dat.df0001$Year
y = dat.df0001$EnsemAnnualNoP_TreeH
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="#67000d", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="#67000d", lty=2)
abline(lm1[[3]], lwd = 3,col = "#67000d")
legend(x = 1980, y = 10+70, legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "#67000d")
# text(x = 2000, y = 0.2*365, labels = expression("Permafrost: 1.00 ± 0.10"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "blue")
	
par(new=TRUE)

plot(EnsemAnnualNoP_TreeL~Year, type = "l", data = dat.df0001, subset = Year >= 1980 & Year <= 2017,
                                  ylim = c(-0.1,0.25)*365,
                                  xlim = c(1980,2017), 
								  cex = 2,  
								  col = "#fcbba1",lwd = 3, 
								  xlab = "Year", ylab = expression("Net C Uptake (" ~ gC ~ m^{-2} ~ yr ^{-1}~ ")"),
								  # xlab = expression(paste('Water Temperature, ',degree,'C',sep ='')),
								  cex.axis = 1.5, cex.lab = 1.6)
	
lmf = function(x,y){

lm1 = lm(y~x)

newx = 1980:2017
conf_interval <- predict(lm1, newdata=data.frame(Year=newx), interval="confidence",level = 0.95)

modsum = summary(lm1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
sen = modsum$coefficients[2,1]
rp = vector('expression',1)
sd1 = modsum$coefficients[2,2]

# rp[1] = substitute(expression(""~ delta == MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

# rp[1] = substitute(expression("Trend:"~ MYVALUE1 ~ gC ~ m^{-2} ~ d ^{-1}~""), 
		# list(MYVALUE1 = format(sen, digits = 2)))[2]

rp[1] = substitute(expression("TC = (0,50) in NoP: "~ MYVALUE1 ~ "±" ~ MYVALUE2 ~ gC ~ m^{-2} ~ yr ^{-2}~""), 
		list(MYVALUE1 = format(sen, digits = 2), MYVALUE2 = format(sd1, digits = 2)))[2]

return(list(newx, conf_interval,lm1,rp))
}
	
x = dat.df0001$Year
y = dat.df0001$EnsemAnnualNoP_TreeL
lm1 = lmf(x=x,y = y)
lines(lm1[[1]], lm1[[2]][,2], col="#fcbba1", lty=2)
lines(lm1[[1]], lm1[[2]][,3], col="#fcbba1", lty=2)
abline(lm1[[3]], lwd = 3,col = "#fcbba1")
legend(x = 1980, y = 10+65, legend = lm1[[4]], bty = 'n', cex = 1.5, text.col = "#fcbba1")
# text(x = 2000, y = 0.2*365, labels = expression("Permafrost: 1.00 ± 0.10"~ gC ~ m^{-2} ~ yr ^{-2}~""), cex = 1.5, col = "blue")

## add statiscis 

	
# inset a subset to seasonsity
d1.df010 = d1.df01

par(new = TRUE)
par(fig = c(0.55, 1, 0.1, 0.5))
v1 = c()
v2 = c()
for(i in 1:16){

lm1 = lm(d1.df010[,i]~d1.df010[,17])
modsum = summary(lm1)
sen = modsum$coefficients[2,1]
v1 = c(v1, sen)

sd1 = modsum$coefficients[2,2]
v2 = c(v2, sd1)

}

v1.df = t(rbind(v1[c(1:4)],v1[c(5:8)],v1[c(9:12)]))

rownames(v1.df) = c("NoP_TreeH", "NoP_TreeL","P_TreeH", "P_TreeL")
colnames(v1.df) = c("EGS", "LGS", "Win")

# plot with sd
v2.df = t(rbind(v2[c(1:4)],v2[c(5:8)],v2[c(9:12)]))

rownames(v2.df) = c("NoP_TreeH", "NoP_TreeL","P_TreeH", "P_TreeL")
colnames(v2.df) = c("EGS", "LGS", "Win")

# us ggplot 2
# http://www.sthda.com/english/wiki/ggplot2-error-NDVIrs-quick-start-guide-r-software-and-data-visualization

v3.df = data.frame(mn = c(v1.df[,1], v1.df[,2], v1.df[,3]), 
                   sd1 = c(v2.df[,1], v2.df[,2], v2.df[,3]),
				   Season = c(rep("EGS",4),rep("LGS",4),rep("Win",4)),
				   Region = rep(c("NoP_TreeH", "NoP_TreeL","P_TreeH", "P_TreeL"),3))
				   

v3.df$Season <- factor(v3.df$Season, levels = c("EGS", "LGS", "Win"))
v3.df$Region <- factor(v3.df$Region, levels = c("NoP_TreeH", "NoP_TreeL","P_TreeH", "P_TreeL"))

require(ggplot2)
p <- ggplot(v3.df, aes(x=Season, y=mn, fill=Region)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mn-sd1, ymax=mn+sd1), width=.2,
                 position=position_dodge(.9)) + 
				     theme_bw() + 
  theme(legend.text = element_text(size = 18))+
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
 xlab("") + 
  ylab("") + 
 # ylab(expression(""~ gC ~ m^{-2} ~ yr ^{-2}~"")) + 
    theme(strip.text.x = element_text(size=18), strip.text.y = element_text(size=18)) + 
  theme(legend.position = "none") + 
  scale_fill_manual(values=c("#67000d","#fcbba1","#08306b","#c6dbef"))

# print(p)


p = p +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

# overlay ggplot2 figure
require(grid)
par(new=TRUE)

print(p, vp=viewport(.78, .25, .45, .3)) #viewport, first two is the x,y coor; last 2 is the inset size

dev.off()




