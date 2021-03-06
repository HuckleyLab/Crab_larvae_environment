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
##################atlantic crop area and nice colors for maps
atlantic.ext <- extent(-104,-28.58, -40,34.69)# atlantic crop for the area I will use - long lat
atlantic.ext2 <- extent(256,331.52, -40,34.69)#atlantic crop for the area I will use - geographic coordinates
my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))#color to use in graphs

##################load SST 

setwd ("C:/Users/Murilo/Desktop/Data/Abiotic_data_set/NOAA/Sea surface daily mean temp")
setwd("/Volumes/GoogleDrive/Team Drives/TrEnCh/Projects/Murilo/Sea_surface_daily_temp/Sea surface daily mean temp")
dir()

#years to analyze, correct years? - yes, it is correct.
years= 1982:2018

#construct file names for first year, Can eventually loop through years
year.ind=1
file.name= paste("sst.day.mean.", years[year.ind],".nc", sep="")

#load data for 1 year
sst<-nc_open(file.name)
print(sst)

#extract temperature data
temp=ncvar_get(sst,"sst")
dim(temp) #dimensions lon, lat, time
print(temp)

#extract lon, lat, time
lons= ncvar_get(sst,"lon") #get info about long
dim(lons)
lats= ncvar_get(sst,"lat") #get info about latitude
dim(lats)
times= ncvar_get(sst,"time") #julian date, calendar day, since 1800-01-01
times
#change to day of year
times1= as.Date(times, origin="1800-01-01")
times1
dim(times1)

#--------------
#subset data to area of interest
lats.sel= which(lats>(-40)&lats<(34.69))#specific area of interest in atlantic ocean
lons.sel= which(lons>(260)&lons<(331.52))#specific area of interest in atlantic ocean

temp.sub=temp[lons.sel,lats.sel,] #using the specified area
dim(temp.sub)
#plot out
image.plot(lons[lons.sel],lats[lats.sel],temp[lons.sel,lats.sel,365],col=my.colors(1000))
data(wrld_simpl)
plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE)
#get rid of Pacific values
temp.sub2<-temp.sub #lon, lat, time
#for each latitude, delete values before NA
for(lat.id in 1:280){
  na1=which.max(is.na(temp.sub2[,lat.id,1]))
  temp.sub2[1:na1,lat.id,]<-NA
  print(lat.id)
}
temp.sub2
dim(temp.sub2)

#update latitudes and longitudes
lons= lons[lons.sel]
lats= lats[lats.sel]
#update temp data to subset
temp=temp.sub2

#image.plot format longitude, latitude, matrix of data
image.plot(lons,lats,temp[,,1],col=my.colors(1000))
data(wrld_simpl)
plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE)
#--------------

#Set species' parameters 
######Menippe nodifrons
To=17.2 ###lower developmental temperature
G=204.08 ### daily growing degree days 

#survival function from below, needs checking and fixing
survi.function=function(day,temp){
  surv= 2.7931035  -3.3789256*day + 1.7749461*day^2 -0.2621086*day^3 + 0.0122824*day^4 -0.0720098*temp +0.1335313*day^1*temp  -0.0693260*day^2*temp  +0.0101202*day^3*temp  -0.0004729*day^4*temp
  #constraint from zero to 1
  surv[surv<0]=0
  surv[surv>1]=1
  return(surv)  
}

#define seasonal timing in terms of days of year
times2=1:90 # selecting days from January 1st to March 31th
times2.1=257:365 #selecting days from Setempber to December 31th    

#make arrays for output
#development
devel.out= array(NA, dim=c(length(lons),length(lats),length(times2)))
#survival
surv.out= array(NA, dim=c(length(lons),length(lats),length(times2)))

for(lon.k in 1:length(lons)){ #loop through longitude
  
  #find latitudes with data
  lats.sel2= which(!is.na(temp[lon.k,,1]))
  
  if(length(lats.sel2)>0){ #check there's cells with data
    for(lat.k in lats.sel2){ #loop through latitudes with data
      
      #estimate growing degree days
      GDDs= temp[lon.k,lat.k,]-To #how many growing degree days per day
      
      #cummulative sum of GDDs
      cumGDDs=cumsum(GDDs) 
      
      #estimate development time starting at different days
      for(time.k in times2){
        
        cumGDDs= cumsum( GDDs[times2[time.k]:365] )
        #substract off starting value
        cumGDDs= cumGDDs-cumGDDs[1]
        
        #find first date the exceed G
        
        devel.out[lon.k,lat.k,time.k]= which.max(cumGDDs>G)+time.k
        
        #estimate survival until reaching G
        day= which.max(cumGDDs>G)
        #mean temp
        ### constrain to end of year, NEED TO FIX
        end.day= times2[time.k]+day
        if(end.day>365)end.day=365
        temps= mean(temp[lon.k,lat.k,times2[time.k]:end.day])
        #estimate survival
        surv.out[lon.k,lat.k,time.k]= survi.function(day, temps)
      
      } #end loop timing
      
    }#end latitude loop
  } #end check for data
} #end longitude loop

devel.out[lon.k,lat.k,]

#plot out doy development assuming development starts on doy 50 
image.plot(lons,lats,devel.out[,,50],col=my.colors(1000))
data(wrld_simpl)
plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE)
#plot out corresponding survival
image.plot(lons,lats,surv.out[,,50],col=my.colors(1000))
data(wrld_simpl)
plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE)


#==========================
###survival data Menippe nodifrons
setwd ("C:/Users/Murilo/Desktop/Data/Abiotic_data_set/NOAA/Sea surface daily mean temp")
#to read from Google Drive
setwd("/Volumes/GoogleDrive/Team Drives/TrEnCh/Projects/Murilo/")
dir()
Survival_mydata<-read.table("survival_Menippe_mydata.txt",header=T)
Survival_Scotto=read.table("survival_Menippe_Scotto_1979.txt",header=T)

survi=Survival_mydata$surv
day=Survival_mydata$day
temperatura=Survival_mydata$temp

#combine data
names(Survival_Scotto)[3]<-"surv"
#turn Scotto data into percent
Survival_Scotto$surv= Survival_Scotto$surv/100
surv.dat= rbind(Survival_Scotto, Survival_mydata)

#plot data and check different polynomial degrees
ggplot(Survival_mydata, aes(x=day, y=surv, color=temp, group=temp)) + geom_point() + 
  geom_smooth(method="lm", formula=y~poly(x, 4))

ggplot(Survival_Scotto, aes(x=day, y=surv, color=temp, group=temp)) + geom_point() + 
  geom_smooth(method="lm", formula=y~poly(x, 4))

#plot combined
ggplot(surv.dat, aes(x=day, y=surv, color=temp, group=temp)) + geom_point() + 
  geom_smooth(method="lm", formula=y~poly(x, 4))

#-------------

mod1= lm(surv~poly(day,4, raw=TRUE)*temp, data=Survival_mydata)

mod1= lm(surv~ day + I(day^2) + I(day^3) + I(day^4) + temp +day:temp + I(day^2):temp + I(day^3):temp + I(day^4):temp, data=Survival_mydata)
#check model         
plot(day, survi)
points(day, predict(mod1), col="red")
#extract coefficients
coefs= coef(mod1)

#make function manually or coudl use the predict() function
survi.function=function(day,temp){
  surv= 2.7931035  -3.3789256*day + 1.7749461*day^2 -0.2621086*day^3 + 0.0122824*day^4 -0.0720098*temp +0.1335313*day^1*temp  -0.0693260*day^2*temp  +0.0101202*day^3*temp  -0.0004729*day^4*temp
  #constraint from zero to 1
  surv[surv<0]=0
  surv[surv>1]=1
  return(surv)  
}


#====================
summary(mod1)
par(mfrow=c(2,2))
plot(mod1)# vizualize the residuals effect and quality of the models R²=0.89
coef.temp=-0.1332
coef.days=9.0840
coef.days:temp=-0.4176

#function using my data 25°C,27°C,29°C until ZoeaIII
survi.function=function(day,temp){
  9.0840*day +(-0.1332)*temp + (-0.4176)*day*temp
}
survi.function
survi.function(2,29)# random values for days and temperature
names(mod1)
mod1$coefficients

survi2=Survival_Scotto$surv
day2=Survival_Scotto$day
temperatura2=Survival_Scotto$temp
mod2=lm(survi2~poly(day2)*temperatura2)
summary(mod2)
mod.coef=print(mod2)
par(mfrow=c(2,2))
plot(mod2)
coef.temp=0.621
coef.days=-63.663
coef.days:temp=-4.555

#function using data from Scotto 1979 25°C,30°C until megalopae
survi.function2=function(day2,temp2){
  (-63.663)*day2 +0.621*temp2 + (-4.555)*day2*temp2
}
survi.function2
survi.function2(2,29)# random values for days and temperature
names(mod2)
mod2$coefficients
