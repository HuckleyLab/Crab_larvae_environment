
#===================LIBRARIES===================
library(raster)
library(ncdf4)
library(fields)
library(colorRamps)
library(sdmpredictors)
library(leaflet)
library(abind)
library(csv)
library(googledrive)
library(maps)
library(maptools)
library(margins)
library(prediction)
library(ggplot2)
#===============================================

my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))#colors to use in graphs

setwd("C:/drive.google.com/drive/folders/1-EKhjU3WpZzlmgk5lD3D-f_NmFWUvRZZ")
setwd("C:/Users/Murilo/Desktop/Future predictions CMIP5")#set directory
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/Murilo/SST future predictions/ACESS1.3 RCP4.5/")
dir()

years=c(2016:2025,2026:2035,2036:2046,2046:2055,2056:2065,2066:2075,2076:2085,2086:2095,2096:2100)

Acess_2026_2035=nc_open("tos_day_ACCESS1-3_rcp45_r1i1p1_20260101-20351231.nc")
print(Acess_2026_2035)

#extract variables tos,lon,lat and time
sea_surface=ncvar_get(Acess_2026_2035,"tos")#get infor bout SST
print(sea_surface)
dim(sea_surface)
sea_surface = sea_surface - 273.15# convert units from kelvin to celcius
lons= ncvar_get(Acess_2026_2035,"i") #get info about long
dim(lons)
lats= ncvar_get(Acess_2026_2035,"j") #get info about latitude
dim(lats)
times= ncvar_get(Acess_2026_2035,"time") #get info about time
dim(times)
times1= as.Date(times, origin="0001-01-01")#convert to Julian date, since 0001-01-01
times1

#subset data to area of interest
lats.sel= which(lats>(80)&lats<(198))#specific area of interest in atlantic ocean
lons.sel= which(lons>(180)&lons<(260))#specific area of interest in atlantic ocean
surface_temp=sea_surface[lons.sel,lats.sel,] #using the specified area
dim(surface_temp)

#plot total area
escala<-(c(0,5,10,15,20,25,30,35))
image.plot(lons,lats,sea_surface[,,3000],col=my.colors(1000), axis.args=list(at=escala, labels=names(escala)))
data(wrld_simpl)
plot(elide(wrld_simpl,shift = c(360, 0)), add = TRUE)#NOT WORKING
map("world2",add=T,fill=TRUE,xlim=c(0,359),ylim=c(0,299), col="grey", bg="grey")#NOT WORKING

#plot area of interest (Tropical Atlantic)
escala<-(c(15,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35))
image.plot(lons[lons.sel],lats[lats.sel],surface_temp[,,3000],col=my.colors(1000), axis.args=list(at=escala, labels=names(escala)))

#Set species' parameters Ucides cordatus
To=12.37 ###lower developmental temperature
G=285.71 ### daily growing degree days 

survi.function=function(day,temp){
  surv= 1.161898  -1.289353e-01*day + 2.447864e-02*day^2 -9.780041e-04*day^3 + 8.937547e-06*day^4 -5.323727e-03*temp +4.684993e-03*day^1*temp  -1.010040e-03*day^2*temp  +3.879397e-05*day^3*temp  -3.207424e-07*day^4*temp
  #constraint from zero to 1
  surv[surv<0]=0
  surv[surv>1]=1
  return(surv)  
}
survi.function

#for each latitude, delete values before NA
for(lat.id in 1:117){
  na1=which.max(is.na(surface_temp[,lat.id,]))
  surface_temp[1:na1,lat.id,]<-NA
  print(lat.id)
}

#plot area of interest (Tropical Atlantic) without pacific
escala<-(c(15,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35))
image.plot(lons[lons.sel],lats[lats.sel],surface_temp[,,3000],col=my.colors(1000), axis.args=list(at=escala, labels=names(escala)))

# selecting days from January 1st to March 31th and October 1st to December 31th of yeach year
times2=c(1:91,275:457,640:821,1005:1186,1370:1553,1737:1917,2101:2282,2466:2647,2831:3013,3197:3378,3562:3653)

t= surface_temp[,,1:91]
t2= surface_temp[,,275:457]
t3= surface_temp[,,640:821]
t4= surface_temp[,,1005:1186]
t5= surface_temp[,,1370:1553]
t6= surface_temp[,,1737:1917]
t7= surface_temp[,,2101:2282]
t8= surface_temp[,,2466:2647]
t9= surface_temp[,,2831:3013]
t10= surface_temp[,,3197:3378]
t11= surface_temp[,,3562:3652]

#abind merge the two files
sst= abind(t,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11, along=3)
dim(sst)

times3=1:1823

#update latitudes and longitudes
lons= lons[lons.sel]
lats= lats[lats.sel]

#make arrays for output
#development
devel.out= array(NA, dim=c(length(lons),length(lats),length(times3)))
#survival
surv.out= array(NA, dim=c(length(lons),length(lats),length(times3)))

for(lon.k in 1:length(lons)){ #loop through longitude
  lats.sel2= which(!is.na(sst[lon.k,,1]))#find latitudes with data
  
  if(length(lats.sel2)>0){ #check there's cells with data
    
    for(lat.k in lats.sel2){ #loop through latitudes with data
      #estimate growing degree days
      GDDs= sst[lon.k,lat.k,]-To #how many growing degree days per day
      #cummulative sum of GDDs
      cumGDDs=cumsum(GDDs) 
      
      for(time.k in times3){ #estimate development time starting at different days
        cumGDDs= cumsum(GDDs[times3[time.k]:1823])####
        #substract off starting value
        cumGDDs= cumGDDs-cumGDDs[1]
        
        #find first date the exceed G
        devel.out[lon.k,lat.k,time.k]= which.max(cumGDDs>G)+time.k
        
        #estimate survival until reaching G
        day= which.max(cumGDDs>G)
        #mean temp
        ### constrain to end of year
        end.day= times3[time.k]+day
        if(end.day>1823)end.day=1823
        temps= mean(sst[lon.k,lat.k,times3[time.k]:end.day])
        #estimate survival
        surv.out[lon.k,lat.k,time.k]= survi.function(day, temps)
      } #end loop timing
    }#end latitude loop
  } #end check for data
} #end longitude loop

#plot out day to development assuming development starts on day 1 
devel.out=replace(devel.out,devel.out>38,NA)
e<-(c(20,21,22,23,24,25,26,27,28,29,30,31,32,32,33,34,35,36,37))
image.plot(lons,lats,devel.out[,,1], col=my.colors(1000),legend.lab = "days to megalopae phase",axis.args=list(at=e, labels=names(e)))
data(wrld_simpl)
plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE)#not working need to fix scale

#plot out corresponding survival
surv.out<-replace(surv.out, surv.out > 0.69 , 0)
s<-(c(0.3,0.35,0.4,0.45,0.5,0.55))
image.plot(lons,lats,surv.out[,,1230],col=my.colors(1000),legend.lab = "survival to megalopae phase",axis.args=list(at=s, labels=names(s)))
data(wrld_simpl)
plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE)

