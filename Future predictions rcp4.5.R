
#########libraries
library(cmsaf)
library(ncdf4)
library(fields)
library(colorRamps)
library(sp)
library(rgdal)
library(raster)
library(abind)
library(maptools)
##########

my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))#colors to use in graphs

#-----------
#LOAD DATA

setwd("E:/Abiotic dataset/SST/future predictions/RCP 4.5")
dir()

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/Murilo/SST future predictions/ACESS1.3 RCP4.5/")
Acess_2016_2025=nc_open("tos_day_ACCESS1-3_rcp45_r1i1p1_20160101-20251231.nc")

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
lats.sel= which(lats> (-28) &lats<(29))#specific area of interest in atlantic ocean
lons.sel= which(lons>(-100)&lons<(-33))#specific area of interest in atlantic ocean
surface_temp=sea_surface[lons.sel,lats.sel,] #using the specified area
dim(surface_temp)

#remove pacific ocean info
for(lat.id in 1:117){
  na1=which.max(is.na(surface_temp[,lat.id,]))
  surface_temp[1:na1,lat.id,]<-NA  
}

#delete +4 cell data from South American Atlantic coast 
nend= dim(surface_temp)[1]
for(lat.id in 1:89){
  na1= min(which(is.na(surface_temp[,lat.id,1])))
  #delete values before first NA
  surface_temp[1:na1,lat.id,]<-NA 
  #find coast
  coast= min(which(!is.na(surface_temp[,lat.id,1]))) 
  #find next island na
  w.isl=min(which(is.na(surface_temp[coast:nend,lat.id,1])))  
  while( is.finite(w.isl) ){
    #keep 4 coastal cells on each side
    if ((coast+w.isl-4)>(coast+4)) surface_temp[(coast+4):(coast+w.isl-4),lat.id,]<-NA
    #find next coast
    ind= coast+w.isl
    coast= min(which(!is.na(surface_temp[ind:nend,lat.id,1]))) 
    coast= ind+coast
    #find next island
    w.isl=min(which(is.na(surface_temp[coast:nend,lat.id,1])))
  }
  #last coastal grid cells
  if(coast+4 <nend) surface_temp[(coast+4):nend,lat.id,]<-NA 
}

#removing cells without islands
for(lat.id in 89){
  na2=which.min(is.na(surface_temp[,lat.id,]))  
  if((na2+96)<dim(surface_temp)[1]) surface_temp[(na2+96):dim(surface_temp)[1],lat.id,]<-NA  
} 
for(lat.id in 90:100){
  na2=which.min(is.na(surface_temp[,lat.id,]))  
  if((na2+26)<dim(surface_temp)[1]) surface_temp[(na2+26):dim(surface_temp)[1],lat.id,]<-NA  
} 
for(lat.id in 101){
  na2=which.min(is.na(surface_temp[,lat.id,]))  
  if((na2+24)<dim(surface_temp)[1]) surface_temp[(na2+24):dim(surface_temp)[1],lat.id,]<-NA  
} 
for(lat.id in 102){
  na2=which.min(is.na(surface_temp[,lat.id,]))  
  if((na2+37)<dim(surface_temp)[1]) surface_temp[(na2+37):dim(surface_temp)[1],lat.id,]<-NA  
}
for(lat.id in 103){
  na2=which.min(is.na(surface_temp[,lat.id,]))  
  if((na2+29)<dim(surface_temp)[1]) surface_temp[(na2+29):dim(surface_temp)[1],lat.id,]<-NA  
}
for(lat.id in 104:109){
  na2=which.min(is.na(surface_temp[,lat.id,]))  
  if((na2+30)<dim(surface_temp)[1]) surface_temp[(na2+30):dim(surface_temp)[1],lat.id,]<-NA  
}
for(lat.id in 110:112){
  na2=which.min(is.na(surface_temp[,lat.id,]))  
  if((na2+28)<dim(surface_temp)[1]) surface_temp[(na2+28):dim(surface_temp)[1],lat.id,]<-NA  
} 
for(lat.id in 113:114){
  na2=which.min(is.na(surface_temp[,lat.id,]))  
  if((na2+23)<dim(surface_temp)[1]) surface_temp[(na2+23):dim(surface_temp)[1],lat.id,]<-NA  
} 
for(lat.id in 115:117){
  na2=which.min(is.na(surface_temp[,lat.id,]))  
  if((na2+21)<dim(surface_temp)[1]) surface_temp[(na2+21):dim(surface_temp)[1],lat.id,]<-NA  
}

#PLOT
#update latitudes and longitudes
lons= lons[lons.sel]
lats= lats[lats.sel]

my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))#colors to use in graphs
e<-(c(15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45))
image.plot(lons,lats,surface_temp[,,1099], col=my.colors(34),legend.lab = "Superficial Sea Temperature",axis.args=list(at=e, labels=names(e)))
data(wrld_simpl)
plot(elide(wrld_simpl), add = TRUE)

#Set species' parameters Ucides cordatus
To=12.37 ###lower developmental temperature
G=285.71 ### daily growing degree days

#survival function
survi.function=function(day,temp){
  surv= 1.161898  -1.289353e-01*day + 2.447864e-02*day^2 -9.780041e-04*day^3 + 8.937547e-06*day^4 -5.323727e-03*temp +4.684993e-03*day^1*temp  -1.010040e-03*day^2*temp  +3.879397e-05*day^3*temp  -3.207424e-07*day^4*temp
  #constraint from zero to 1
  surv[surv<0]=0
  surv[surv>1]=1
  return(surv)  
}
survi.function

# selecting days from January 1st to March 31th and October 1st to December 31th of yeach year
times2=c(1:121,275:485,640:850,1005:1216,1370:1583,1737:1946,2101:2311,2466:2677,2831:3043,3197:3407,3562:3653)
#select seasons, start and end dates
#inds=cbind(c(275, 640, 1005, 1370, 1737, 2101, 2466, 2831, 3197, 3562),c(457,821,1186,1553,1917,2282,2647,3013,3378,3653))
inds=cbind(c(274, 639, 1004, 1370, 1735, 2100, 2465, 2831, 3196),c(455,820,1186,1551,1916,2281,2647,3012,3377))
years.inds= years[times2]
doys.inds= doys[times2]

times3=c(275:366,1:121)

#save surface_temp
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/Murilo/")
saveRDS(surface_temp, "surface_temp.rds")

#make arrays for output, add dimension for future periods
#development
devel.out= array(NA, dim=c(length(lons),length(lats),length(times3),nrow(inds) ))
#survival
surv.out= array(NA, dim=c(length(lons),length(lats),length(times3),nrow(inds) ))

for(lon.k in 1:length(lons)){ #loop through longitude
  lats.sel2= which(!is.na(surface_temp[lon.k,,1]))#find latitudes with data
  
  if(length(lats.sel2)>0){ #check there's cells with data
    
    for(lat.k in lats.sel2){ #loop through latitudes with data
      
      for(yr.k in 1:nrow(inds)){ #loop through seasons
        
        #extract start and end dates
        inds.seas= inds[yr.k,1]:inds[yr.k,2] 
        
        #estimate growing degree days
        GDDs= surface_temp[lon.k,lat.k,inds.seas]-To #how many growing degree days per day
        #cummulative sum of GDDs
        cumGDDs=cumsum(GDDs) 
        
        for(time.k in 1:length(inds.seas)){ #estimate development time starting at different days
          
          #substract off starting value
          cumGDDs.t= cumGDDs-cumGDDs[time.k]
          
          #find first date the exceed G
          dev.ind= which.max(cumGDDs.t>G) - time.k  #development duration
          if(dev.ind>0) devel.out[lon.k,lat.k,time.k,yr.k]= dev.ind
          #To sAVE: DOY days.inds[inds.seas[dev.ind]]
          
          #estimate survival until reaching G
          #mean temp
          if(dev.ind>1){#if complete development
            temps= mean(surface_temp[lon.k,lat.k,inds.seas[time.k]:inds.seas[dev.ind+time.k]])
            #estimate survival
            surv.out[lon.k,lat.k,time.k,yr.k]= survi.function(dev.ind, temps)
          } #end check complete development
        } #end loop timing
      } #end loop seasons
    }#end latitude loop
  } #end check for data
} #end longitude loop

#devel.out=replace(devel.out,devel.out<15,NA)#replace values below the minimum possible time for develop by NA
devel.out2=devel.out[,,1:183,]#removing extra days (not season days)
surv.out2=surv.out[,,1:183,]#removing extra days (not season days)
times4=c(275:366,1:91)#removing extra days (not season days)

#plot out day to development assuming development starts on day 1 
e<-(c(-120,-40,-5,15,17,19,21,23,25,27,29,31,33,35,37,39))
image.plot(lons,lats,devel.out2[,,160,1], col=my.colors(34),legend.lab = "days to megalopae phase",axis.args=list(at=e, labels=names(e)))
data(wrld_simpl)
plot(elide(wrld_simpl), add = TRUE)

#plot out corresponding survival
s<-(c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1))
image.plot(lons,lats,surv.out2[,,155,1],col=my.colors(1000),legend.lab = "survival to megalopae phase",axis.args=list(at=s, labels=names(s)))
data(wrld_simpl)
plot(elide(wrld_simpl), add = TRUE)

#average surv, devel and temp by column/row per year
surv.lat2= apply(surv.out2, MARGIN=c(1,2), FUN=mean, na.rm=TRUE)
s<-(c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1))
image.plot(lons,lats,surv.lat2,col=my.colors(1000),legend.lab = "survival to megalopae phase",axis.args=list(at=s, labels=names(s)))
data(wrld_simpl)
plot(elide(wrld_simpl), add = TRUE)

devel.lat2= apply(devel.out2, MARGIN=c(1,2), FUN=mean, na.rm=TRUE)
e<-(c(15,17,19,21,23,25,27,29,31,33,35,37,39))
image.plot(lons,lats,devel.lat2,col=my.colors(1000),legend.lab = "Devel time",axis.args=list(at=e, labels=names(s)))
data(wrld_simpl)
plot(elide(wrld_simpl), add = TRUE)

temp.out=array(surface_temp, dim=c(length(lons),length(lats),length(times4),nrow(inds)))
temp.lat2= apply(temp.out, MARGIN=c(1,2), FUN=mean, na.rm=TRUE)
e<-(c(22,24,26,28,30))
image.plot(lons,lats,temp.lat2,col=my.colors(1000),legend.lab = "Temp",axis.args=list(at=e, labels=names(s)))
data(wrld_simpl)
plot(elide(wrld_simpl), add = TRUE)