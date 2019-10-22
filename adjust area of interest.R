
#libraries
library(ncdf4)
library(fields)
library(abind)
library(maptools)

#set directory
setwd ("E:/Abiotic dataset/SST/historic data/NOAA/Sea surface daily mean temp")
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/Murilo/Sea_surface_daily_temp/Sea surface daily mean temp/")

#load data
data.1982=nc_open("sst.day.mean.1982.nc")
data.1983=nc_open("sst.day.mean.1983.nc")

#extract sst, lon, lat, time fisrt year
temp=ncvar_get(data.1982,"sst")
lons= ncvar_get(data.1982,"lon") #get info about long
lats= ncvar_get(data.1982,"lat") #get info about latitude
times= ncvar_get(data.1982,"time") #julian date, calendar day, since 1800-01-01
times1= as.Date(times, origin="1800-01-01")#change to day of year
times1

#extract sst, lon, lat, time second year
temp2=ncvar_get(data.1983,"sst")
lons2= ncvar_get(data.1983,"lon") #get info about long
lats2= ncvar_get(data.1983,"lat") #get info about latitude
times2= ncvar_get(data.1983,"time") #julian date, calendar day, since 1800-01-01
times3= as.Date(times2, origin="1800-01-01")#change to day of year
times3

#convert longituide
convert.lon= function(r0) ifelse(r0 > 180, -360 + r0, r0)
lons= convert.lon(lons)
lons2= convert.lon(lons2)

#subset data to area of interest ##
lats.sel= which(lats> (-28) &lats<(29))#specific area of interest in atlantic ocean
lons.sel= which(lons>(-100)&lons<(-33))#specific area of interest in atlantic ocean
surface_temp=temp[lons.sel,lats.sel,] #using the specified area
dim(surface_temp)

lats.sel2= which(lats2> (-28) &lats2<(29))#specific area of interest in atlantic ocean
lons.sel2= which(lons2>(-100)&lons2<(-33))#specific area of interest in atlantic ocean
surface_temp2=temp2[lons.sel2,lats.sel2,] #using the specified area
dim(surface_temp2)

#subset
surface_temp= surface_temp[,,275:365]
surface_temp2= surface_temp2[,,1:121]

#abind merge the two files
temp.82_83= abind(surface_temp, surface_temp2, along=3)
dim(temp.82_83)

#temp.82_83_copy= temp.82_83
#temp.82_83= temp.82_83_copy

#for each latitude, after the last NA delete values after the 7?grid cell (approx. 200km)

nend= dim(temp.82_83)[1]

for(lat.id in 1:228){
  na1= min(which(is.na(temp.82_83[,lat.id,1])))
  #delete values before first NA
  temp.82_83[1:na1,lat.id,]<-NA 
  
  #find coast
  coast= min(which(!is.na(temp.82_83[,lat.id,1]))) 
  
  #find next island na
  w.isl=min(which(is.na(temp.82_83[coast:nend,lat.id,1])))  
  
  while( is.finite(w.isl) ){
  #keep 7 coastal cells on each side
  if ((coast+w.isl-7)>(coast+7)) temp.82_83[(coast+7):(coast+w.isl-7),lat.id,]<-NA
  #find next coast
  ind= coast+w.isl
  coast= min(which(!is.na(temp.82_83[ind:nend,lat.id,1]))) 
  coast= ind+coast
  #find next island
  w.isl=min(which(is.na(temp.82_83[coast:nend,lat.id,1])))
  }
  
  #last coastal grid cells
   if(coast+7 <nend) temp.82_83[(coast+7):nend,lat.id,]<-NA 
}

#PLOT
#update latitudes and longitudes
lons= lons[lons.sel]
lats= lats[lats.sel]

#plot
my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))#colors to use in graphs
e<-(c(15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45))
image.plot(lons,lats,temp.82_83[,,30], col=my.colors(34),legend.lab = "Superficial Sea Temperature",axis.args=list(at=e, labels=names(e)))
data(wrld_simpl)
plot(elide(wrld_simpl), add = TRUE)




