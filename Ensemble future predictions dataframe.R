library(cmsaf)
library(ncdf4)
library(fields)
library(colorRamps)
library(sp)
library(rgdal)
library(raster)

my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))#colors to use in graphs

#-----------
#LOAD DATA

setwd("E:/Abiotic dataset/SST/future predictions/RCP 4.5")
dir()

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/Murilo/SST future predictions/ACESS1.3 RCP4.5/")
Acess_2016_2025=nc_open("tos_day_ACCESS1-3_rcp45_r1i1p1_20160101-20251231.nc")
print(Acess_2016_2025)

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/Murilo/SST future predictions/MRI-CGCM3 RCP4.5/")
mri_2016_2025=nc_open("tos_day_MRI-CGCM3_rcp45_r1i1p1_20160101-20251231.nc")
print(mri_2016_2025)

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/Murilo/SST future predictions/CNRM-CM5 RCP4.5/")
cnrm_2016_2025=nc_open("tos_day_CNRM-CM5_rcp45_r1i1p1_20160101-20251231.nc")
print(cnrm_2016_2025)

setwd("E:/Abiotic dataset/SST/future predictions/RCP 8.5")
dir()

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/Murilo/SST future predictions/ACESS1.3 RCP8.5/")
Acess_2016_2025_8.5=nc_open("tos_day_ACCESS1-3_rcp85_r1i1p1_20160101-20251231.nc")
print(Acess_2016_2025_8.5)

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/Murilo/SST future predictions/MRI-CGCM3 RCP8.5/")
mri_2016_2025_8.5=nc_open("tos_day_MRI-CGCM3_rcp85_r1i1p1_20160101-20251231.nc")
print(mri_2016_2025_8.5)

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/Murilo/SST future predictions/CNRM-CM5 RCP8.5/")
cnrm_2016_2025_8.5=nc_open("tos_day_CNRM-CM5_rcp85_r1i1p1_20160101-20251231.nc")
print(cnrm_2016_2025_8.5)

#------------
#Acess

#extract variables tos,lon,lat and time
sea_surface=ncvar_get(Acess_2016_2025,"tos")#get infor bout SST
print(sea_surface)
dim(sea_surface)
sea_surface = sea_surface - 273.15# convert units from kelvin to celcius
lons= ncvar_get(Acess_2016_2025,"lon") #get info about long
lons= lons[,1]
lats= ncvar_get(Acess_2016_2025,"lat") #get info about latitude
lats= lats[1,]
times= ncvar_get(Acess_2016_2025,"time") #get info about time
dim(times)
times1= as.Date(times, origin="0001-01-01")#convert to Julian date, since 0001-01-01
times1

#estimate years and days
years=as.numeric(format(times1,"%Y"))
doys= as.numeric(format(times1,"%j"))

#convert longituide
convert.lon= function(r0) ifelse(r0 > 180, -360 + r0, r0)
lons= convert.lon(lons)

#subset data to area of interest ##
lats.sel= which(lats> (-35) &lats<(30))#specific area of interest in atlantic ocean
lons.sel= which(lons>(-100)&lons<(-33))#specific area of interest in atlantic ocean
surface_temp=sea_surface[lons.sel,lats.sel,] #using the specified area
dim(surface_temp)

#---
#for each latitude, delete values before NA
for(lat.id in 1:125){
  na1=which.max(is.na(surface_temp[,lat.id,]))
  surface_temp[1:na1,lat.id,]<-NA
}

# selecting days from January 1st to March 31th and October 1st to December 31th of yeach year
times2=c(1:91,275:457,640:821,1005:1186,1370:1553,1737:1917,2101:2282,2466:2647,2831:3013,3197:3378,3562:3653)
#select seasons, start and end dates
#inds=cbind(c(275, 640, 1005, 1370, 1737, 2101, 2466, 2831, 3197, 3562),c(457,821,1186,1553,1917,2282,2647,3013,3378,3653))
inds=cbind(c(274, 639, 1004, 1370, 1735, 2100, 2465, 2831, 3196),c(455,820,1186,1551,1916,2281,2647,3012,3377))
years.inds= years[times2]
doys.inds= doys[times2]

#sst= surface_temp[,,1:91]

times3=c(275:366,1:91)

#update latitudes and longitudes
lons= lons[lons.sel]
lats= lats[lats.sel]

#plot
escala<-(c(16,18,20,22,24,26,28,30,32))
image.plot(lons,lats,surface_temp[,,1700],col=my.colors(1000), axis.args=list(at=escala, labels=names(escala)))


#----------
#MRI

#extract variables tos,lon,lat and time#need to fix
sea_surface.mri=ncvar_get(mri_2016_2025,"tos")#get infor bout SST
print(sea_surface.mri)
dim(sea_surface.mri)
sea_surface.mri=replace(sea_surface.mri,sea_surface.mri<1,NA)#replace 0 values from continent areas by NA
sea_surface.mri = sea_surface.mri - 273.15# convert units from kelvin to celcius
lons.mri= ncvar_get(mri_2016_2025,"lon") #get info about long
lons.mri= lons.mri[,1]
dim(lons.mri)
lats.mri= ncvar_get(mri_2016_2025,"lat") #get info about latitude
lats.mri= lats.mri[1,]
dim(lats.mri)
times.mri= ncvar_get(mri_2016_2025,"time") #get info about time
dim(times.mri)
times1.mri= as.Date(times.mri, origin="1850-01-01")#convert to Julian date, since 0001-01-01
times1.mri

#check dimensions world plot
escala<-(c(0,2,4,6,16,18,20,22,24,26,28,30,32))
image.plot(sea_surface.mri[,,17],col=my.colors(1000), axis.args=list(at=escala, labels=names(escala)))

# ROTATE GRID####need to fix
#Acess_stack=stack("tos_day_ACCESS1-3_rcp45_r1i1p1_20160101-20251231.nc")
#extent(Acess_stack)

#nsat<- stack ("tos_day_MRI-CGCM3_rcp45_r1i1p1_20160101-20251231.nc")
#extent(nsat)##check the extent## this will be in the form 0-360 degrees

#nsat1<-rotate(nsat)#change the coordinates
#extent(nsat1)#check result:##this should be in the format: -180/180

#convert longituide
lons.mri= convert.lon(lons.mri)

#estimate years and days
years.mri=as.numeric(format(times1.mri,"%Y"))
doys.mri= as.numeric(format(times1.mri,"%j"))

#----------
#FIND MRI VALUES THAT CORRESPOND TO ACESS GRID# (lats.mri>(-29)&lats.mri<(28)-lons.mri>(179)&lons.mri<(249))
sst= surface_temp

#make arrays for mri output
sst.mri= array(NA, dim= dim(surface_temp))

#match times
match.times2= match(times1, times1.mri)

#matching function
findit = function(x,vec){ 
  y = abs(vec - x) 
  if(all(is.na(y)))NA else which.min(y) }

#find closest lats and lons
match.lats2= sapply(lats,findit,lats.mri) 
match.lons2= sapply(lons,findit,lons.mri)

sst.mri= sea_surface.mri[match.lons2, match.lats2, match.times2]

#for each latitude, delete values before NA
for(lat.id in 1:125){
  na1=which.max(is.na(sst.mri[,lat.id,]))
  sst.mri[1:na1,lat.id,]<-NA
}

#plot
image.plot(lons,lats,sst.mri[,,17],col=my.colors(1000))

#----------
#CNRM

#extract variables tos,lon,lat and time
sea_surface.cnrm=ncvar_get(cnrm_2016_2025,"tos")#get infor bout SST
print(sea_surface.cnrm)
dim(sea_surface.cnrm)
sea_surface.cnrm = sea_surface.cnrm - 273.15# convert units from kelvin to celcius
lons.cnrm= ncvar_get(cnrm_2016_2025,"lon") #get info about long
dim(lons.cnrm)
lons.cnrm= lons.cnrm[,1]
lats.cnrm= ncvar_get(cnrm_2016_2025,"lat") #get info about latitude
dim(lats.cnrm)
lats.cnrm= lats.cnrm[1,]
times.cnrm= ncvar_get(cnrm_2016_2025,"time") #get info about time
dim(times.cnrm)
times1.cnrm= as.Date(times.cnrm, origin="2006-01-01")#convert to Julian date, since 0001-01-01
times1.cnrm

#convert longituide
lons.cnrm= convert.lon(lons.cnrm)

#estimate years and doys
years.cnrm=as.numeric(format(times1.cnrm,"%Y"))
doys.cnrm= as.numeric(format(times1.cnrm,"%j"))

#----------
#FIND CNRM VALUES THAT CORRESPOND TO ACESS GRID (lats.cnrm>(90)&lats.cnrm<(204)-lons.cnrm>(188)&lons.cnrm<(258))
sst= surface_temp

#make arrays for mri output
sst.cnrm= array(NA, dim= dim(surface_temp))

#match times
match.times= match(times1, times1.cnrm)

#find closest lats and lons
match.lats= sapply(lats,findit,lats.cnrm) 
match.lons= sapply(lons,findit,lons.cnrm) 

sst.cnrm= sea_surface.cnrm[match.lons, match.lats, match.times]

#for each latitude, delete values before NA
for(lat.id in 1:125){
  na1=which.max(is.na(sst.cnrm[,lat.id,]))
  sst.cnrm[1:na1,lat.id,]<-NA
}

#plot
image.plot(lons,lats,sst.cnrm[,,1700],col=my.colors(1000))


