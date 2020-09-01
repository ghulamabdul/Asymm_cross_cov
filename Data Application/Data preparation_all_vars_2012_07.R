##### Preparing a trivariate data for PM2.5, WS1, WS2 #####
####################################################
############ Loading required libraries ############
####################################################
library(sp)
library(maps)
library(maptools)
library(geosphere)
library(fields)
library(MASS)
#library(scoringRules)
#library(doParallel)
library(rgdal)
#registerDoParallel(cores = 5)
#########################################################
################ Data Creation portion ##################
#########################################################

#setwd("/Users/qadirga/Documents/Project 2/Rcode/Data Analysis")
#setwd("/Users/qadirga/Documents/Project 3/Data Analysis (WS and PM 2.5)/Data preparation")
#setwd("/Users/qadirga/Documents/Project 3/Manuscript/Data Analysis/New day checking")
#setwd("D:/Projects 3/Data Analysis version 2/Data Analysis-gpawspm2.5/Summer data set")
setwd("/Users/qadirga/Documents/Project 3/Manuscript/Data Analysis Final Version")
dir()
#------ PM25 Data ------#

time_pm25_day_read <- system.time(pm25_day_2012 <- read.table("2012_pm25_daily_average.txt", header = TRUE, sep = ","))
time_pm25_day_read

str(pm25_day_2012); head(pm25_day_2012)
# tmp <- pm25_day_2006[pm25_day_2006$FIPS==1001020100,]

pm25_day_2012$Month <- as.numeric(substr(as.character(pm25_day_2012$Date),6,7))
unique(pm25_day_2012$Month)
time_pm25_mon <- system.time(pm25_mon_2012 <- aggregate(pm25_daily_average ~ Month + FIPS + Longitude + Latitude, data=pm25_day_2012, mean))
time_pm25_mon

str(pm25_mon_2012); head(pm25_mon_2012)
names(pm25_mon_2012)[names(pm25_mon_2012)=="pm25_daily_average"] <- "pm25_monthly_mean"
str(pm25_mon_2012); head(pm25_mon_2012) # 867396 obs

pm25_201207 <- subset(pm25_mon_2012, Month==7)
write.csv(pm25_201207, file = "my_pm25_201207_test.csv", row.names = FALSE)
pm25_201207<-read.csv("my_pm25_201207_test.csv")
par(mfrow=c(1,2))


#_______ Data over whole USA _______#
par(mfrow=c(1,1))
quilt.plot(pm25_201207$Longitude,pm25_201207$Latitude,pm25_201207$pm25_monthly_mean,nx=200,ny=200)


#-------------------------------------  WS DATA--------------------------------------------------------------#

# read raw data ####
grib_GDAL_201207 <- readGDAL("narrmon-a_221_20120701_0000_000.grb")


#==============================================================================================####
# transform LCC coordinate to Long/Lat ####
# str(grib_GDAL_200606)
grib_GDAL_201207@proj4string

# str(grib_GDAL_200606@data)
#dim(grib_GDAL_200606@data)

# longlat_200606 <- spTransform(SpatialPoints(coordinates(grib_GDAL_200606), proj4string=grib_GDAL_200606@proj4string), CRS("+proj=longlat +datum=WGS84"))
longlat <- spTransform(SpatialPoints(coordinates(grib_GDAL_201207), proj4string=grib_GDAL_201207@proj4string), CRS("+proj=longlat +datum=WGS84"))
dim(coordinates(grib_GDAL_201207))
dim(coordinates(longlat))

#==============================================================================================####
#wind speed u component at pressure 1000 =277
#wind speed v component at pressure 1000 =323
#WS =wind speed (magnitude)
#geopotential height=57 (pressure 1000) (GPH)
#Total Cloud cover TCDC =214
#Surface temperature : 231
# Temperature at pressure 1000 =232
# Total Precipitation (APCP) = 4
# Total Precipitation nearest gridpoint =5
# RH2m (Relative humidity) = 156
# RH hybrid = 157
# RH isotherm = 158

narr_201207 <- data.frame(Month=201207, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          WS1=grib_GDAL_201207@data[,277], WS2= grib_GDAL_201207@data[,323],
                          WS=sqrt(grib_GDAL_201207@data[,277]^2+grib_GDAL_201207@data[,323]^2),
                          GPH=grib_GDAL_201207@data[,57],
                          TC=grib_GDAL_201207@data[,214],
                          SrfTemp=grib_GDAL_201207@data[,231],
                          Temp_prs1000=grib_GDAL_201207@data[,232],
                          TPrp=grib_GDAL_201207@data[,4],
                          TPrpN=grib_GDAL_201207@data[,5],
                          RH2m=grib_GDAL_201207@data[,156],
                          RHhyb=grib_GDAL_201207@data[,157],
                          RHiso=grib_GDAL_201207@data[,158],
                          AirTemp=grib_GDAL_201207@data[,263])

write.csv(narr_201207, file = "my_narr_201207_test.csv", row.names = FALSE)


##### ______ Merging data files ________#######

#pm25_201201 <- read.csv("my_pm25_201201.csv")
#narr_201201 <- read.csv("my_narr_201201.csv")

#==============================================================================================####
long_pm25<- unique(pm25_201207$Longitude)
lat_pm25 <- unique(pm25_201207$Latitude)

data_201207 <- subset(narr_201207, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))

#==============================================================================================####

LL_pm25 <- cbind(pm25_201207$Longitude, pm25_201207$Latitude)
LL_narr <- cbind(data_201207$Longitude, data_201207$Latitude)
dim(LL_pm25); dim(LL_narr)
par(mfrow=c(3,4))

quilt.plot(x=pm25_201207$Longitude,y=pm25_201207$Latitude,z=pm25_201207$pm25_monthly_mean,nx=200,ny=200,main="PM2.5")
#quilt.plot(x=data_201206$Longitude,y=data_201206$Latitude,z=data_201206$WS1,main="U-component (ws)",nx=150,ny=150)
#quilt.plot(x=data_201206$Longitude,y=data_201206$Latitude,z=data_201206$WS2,main="V-Component (ws)",nx=150,ny=150)
quilt.plot(x=data_201207$Longitude,y=data_201207$Latitude,z=data_201207$WS,main="Magnitude Windspeed",nx=150,ny=150)
quilt.plot(x=data_201207$Longitude,y=data_201207$Latitude,z=data_201207$GPH,main="Geopotential Height",nx=150,ny=150)
quilt.plot(x=data_201207$Longitude,y=data_201207$Latitude,z=data_201207$TC,main="Total Cloud",nx=150,ny=150)
quilt.plot(x=data_201207$Longitude,y=data_201207$Latitude,z=data_201207$SrfTemp,main="Surface Temperature",nx=150,ny=150)
quilt.plot(x=data_201207$Longitude,y=data_201207$Latitude,z=data_201207$Temp_prs1000,main="Temperature at pressure 1000",nx=150,ny=150)
quilt.plot(x=data_201207$Longitude,y=data_201207$Latitude,z=data_201207$TPrp,main="Total Precipitation",nx=150,ny=150)
quilt.plot(x=data_201207$Longitude,y=data_201207$Latitude,z=data_201207$TPrpN,main="Total Precipitation (Nearest Grid Point)",nx=150,ny=150)
quilt.plot(x=data_201207$Longitude,y=data_201207$Latitude,z=data_201207$RH2m,main="Relative Humidity (2m)",nx=150,ny=150)
quilt.plot(x=data_201207$Longitude,y=data_201207$Latitude,z=data_201207$RHhyb,main="Relative Humidity hybrid",nx=150,ny=150)
quilt.plot(x=data_201207$Longitude,y=data_201207$Latitude,z=data_201207$RHiso,main="Relative Humidity isotherm",nx=150,ny=150)
quilt.plot(x=data_201207$Longitude,y=data_201207$Latitude,z=data_201207$AirTemp,main="Relative Humidity isotherm",nx=150,ny=150)

system.time(DMatrix <- distm(LL_pm25,LL_narr)); dim(DMatrix)    # 1398.215 sec
system.time(tmp <- apply(DMatrix, 1, min))                      # 20.734 sec
system.time(min_ind <- which(DMatrix==tmp, arr.ind=T))          # 3.275 sec
system.time(id_data <- min_ind[order(min_ind[,1]),2])           # 0.002 sec

summary(min_ind[,1]); summary(min_ind[,2])
length(unique(min_ind[,1])); length(unique(min_ind[,2]))

#==============================================================================================####
pm25_201207$id_data <- id_data
data_pm25_201207 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201207, mean)

#==============================================================================================####
n_narr <- dim(LL_narr)[1]

data_201207$id_data <- 1:n_narr
dir()
#==============================================================================================####
total_201207 <- merge(data_201207, data_pm25_201207, by="id_data"); str(total_201207)
write.csv(total_201207[,-1], file = "my_data_pm25_201207_testv2.csv", row.names = FALSE)
total_201207<-read.csv("my_data_pm25_201207_testv2.csv")
par(mfrow=c(3,4))
quilt.plot(x=total_201207$Longitude,y=total_201207$Latitude,z=total_201207$pm25_monthly_mean,nx=200,ny=200,main="PM2.5")
#quilt.plot(x=total_201206$Longitude,y=total_201206$Latitude,z=total_201206$WS1,main="U-component (ws)",nx=150,ny=150)
#quilt.plot(x=total_201206$Longitude,y=total_201206$Latitude,z=total_201206$WS2,main="V-Component (ws)",nx=150,ny=150)
quilt.plot(x=total_201207$Longitude,y=total_201207$Latitude,z=total_201207$WS,main="Magnitude Windspeed",nx=150,ny=150)
quilt.plot(x=total_201207$Longitude,y=total_201207$Latitude,z=total_201207$GPH,main="Geopotential Height",nx=150,ny=150)
quilt.plot(x=total_201207$Longitude,y=total_201207$Latitude,z=total_201207$TC,main="Total Cloud",nx=150,ny=150)
quilt.plot(x=total_201207$Longitude,y=total_201207$Latitude,z=total_201207$SrfTemp,main="Surface Temperature",nx=150,ny=150)
quilt.plot(x=total_201207$Longitude,y=total_201207$Latitude,z=total_201207$Temp_prs1000,main="Temperature at pressure 1000",nx=150,ny=150)
quilt.plot(x=total_201207$Longitude,y=total_201207$Latitude,z=total_201207$TPrp,main="Total Precipitation",nx=150,ny=150)
quilt.plot(x=total_201207$Longitude,y=total_201207$Latitude,z=total_201207$TPrpN,main="Total Precipitation (Nearest Grid Point)",nx=150,ny=150)
quilt.plot(x=total_201207$Longitude,y=total_201207$Latitude,z=total_201207$RH2m,main="Relative Humidity (2m)",nx=150,ny=150)
quilt.plot(x=total_201207$Longitude,y=total_201207$Latitude,z=total_201207$RHhyb,main="Relative Humidity hybrid",nx=150,ny=150)
quilt.plot(x=total_201207$Longitude,y=total_201207$Latitude,z=total_201207$RHiso,main="Relative Humidity isotherm",nx=150,ny=150)
#setwd("/Users/qadirga/Documents/Project 3/Manuscript/Data Analysis/New day checking/New region checking")

#part 2 
#####
#par(mfrow=c(1,2))
#quilt.plot(total_201206$Longitude,total_201206$Latitude,total_201206$WS1)
#quilt.plot(total_201206$Longitude,total_201206$Latitude,total_201206$pm25_monthly_mean)

##### -------- Assigning Climatic region --------########
latlong2state <- function(pointsDF) {
  # Prepare SpatialPolygons object with one SpatialPolygon
  # per state (plus DC, minus HI & AK)
  states <- map('state', fill=TRUE, col="transparent", plot=FALSE)
  IDs <- sapply(strsplit(states$names, ":"), function(x) x[1])
  states_sp <- map2SpatialPolygons(states, IDs=IDs,
                                   proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Convert pointsDF to a SpatialPoints object 
  pointsSP <- SpatialPoints(pointsDF, 
                            proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Use 'over' to get _indices_ of the Polygons object containing each point 
  indices <- over(pointsSP, states_sp)
  
  # Return the state names of the Polygons object containing each point
  stateNames <- sapply(states_sp@polygons, function(x) x@ID)
  stateNames[indices]
}


# Assign States
coords <- as.data.frame(cbind(total_201207$Longitude,total_201207$Latitude))
State <- latlong2state(coords)



#==============================================================================================####
unique(as.factor(State))
summary(as.factor(State))

id_NA <- (1:length(State))[is.na(State)]
id_ST <- (1:length(State))[!is.na(State)]
coords_NA <- coords[id_NA,]
coords_ST <- coords[id_ST,]
dim(coords_NA); dim(coords_ST)

DMatrix <- distm(coords_NA,coords_ST); dim(DMatrix)
tmp <- apply(DMatrix, 1, min)
min_ind <- which(DMatrix==tmp, arr.ind=T)
NA2ST <- id_ST[min_ind[order(min_ind[,1]),2]]

State[id_NA] <- State[NA2ST]
unique(as.factor(State))
summary(as.factor(State))
#==============================================================================================####
# Assing Climatic Regions (CR)
CR_NW <- c("washington", "oregon", "idaho")
CR_W <- c("california", "nevada")
CR_SW <- c("utah", "colorado", "arizona", "new mexico")
CR_WNC <- c("montana", "wyoming", "north dakota", "south dakota", "nebraska")
CR_ENC <- c("minnesota", "iowa", "wisconsin", "michigan")
CR_S <- c("kansas", "oklahoma", "texas", "arkansas", "louisiana", "mississippi")
CR_C <- c("illinois", "indiana", "ohio", "missouri", "kentucky", "west virginia", "tennessee")
CR_NE <- c("maine", "new hampshire", "vermont", "new york", "massachusetts", "connecticut", "rhode island",
           "pennsylvania", "new jersey", "delaware", "maryland")
CR_SE <- c("south carolina", "georgia", "alabama", "florida", "north carolina", "virginia")

CR <- rep(NA,length(State))
CR[State %in% CR_NW] <- "NW"
CR[State %in% CR_W]  <- "W"
CR[State %in% CR_SW] <- "SW"
CR[State %in% CR_WNC]<- "WNC"
CR[State %in% CR_ENC]<- "ENC"
CR[State %in% CR_S]  <- "S"
CR[State %in% CR_C]  <- "C"
CR[State %in% CR_NE] <- "NE"
CR[State %in% CR_SE] <- "SE"

summary(as.factor(CR))

#==============================================================================================####

total_201207$State<-State; total_201207$CR<-CR

coords_NW <- coords[CR == "NW",];  text_NW <- (apply(coords_NW,2,min) + apply(coords_NW,2,max))/2
coords_W  <- coords[CR == "W",];   text_W  <- (apply(coords_W,2,min) + apply(coords_W,2,max))/2
coords_SW <- coords[CR == "SW",];  text_SW <- (apply(coords_SW,2,min) + apply(coords_SW,2,max))/2
coords_WNC<- coords[CR == "WNC",]; text_WNC<- (apply(coords_WNC,2,min) + apply(coords_WNC,2,max))/2
coords_ENC<- coords[CR == "ENC",]; text_ENC<- (apply(coords_ENC,2,min) + apply(coords_ENC,2,max))/2
coords_S  <- coords[CR == "S",];   text_S  <- (apply(coords_S,2,min) + apply(coords_S,2,max))/2
coords_C  <- coords[CR == "C",];   text_C  <- (apply(coords_C,2,min) + apply(coords_C,2,max))/2
coords_NE <- coords[CR == "NE",];  text_NE <- (apply(coords_NE,2,min) + apply(coords_NE,2,max))/2
coords_SE <- coords[CR == "SE",];  text_SE <- (apply(coords_SE,2,min) + apply(coords_SE,2,max))/2

# text_NW <- apply(coords_NW,2,median)
text_W  <- apply(coords_W,2,median)
text_SW <- apply(coords_SW,2,median)
text_WNC<- (apply(coords_WNC,2,median) + text_WNC)/2
text_ENC<- apply(coords_ENC,2,median)
text_S  <- apply(coords_S,2,median)
text_C  <- apply(coords_C,2,median)
text_NE <- apply(coords_NE,2,median)
text_SE <- apply(coords_SE,2,median)



map('state'); 
points(total_201207$Longitude, total_201207$Latitude,  col = "red", cex = .1)
map('state', region = CR_NW, lwd=3, interior = FALSE, add = TRUE); text(text_NW[1], text_NW[2], "NW", cex=1.7)
map('state', region = CR_W, lwd=3, interior = FALSE, add = TRUE);  text(text_W[1]+1.5, text_W[2], "W", cex=1.7)
map('state', region = CR_SW, lwd=3, interior = FALSE, add = TRUE); text(text_SW[1], text_SW[2], "SW", cex=1.7)
map('state', region = CR_WNC, lwd=3, interior = FALSE, add = TRUE);text(text_WNC[1], text_WNC[2], "WNC", cex=1.7)
map('state', region = CR_ENC, lwd=3, interior = FALSE, add = TRUE);text(text_ENC[1], text_ENC[2], "ENC", cex=1.7)
map('state', region = CR_S, lwd=3, interior = FALSE, add = TRUE);  text(text_S[1], text_S[2], "S", cex=1.7)
map('state', region = CR_C, lwd=3, interior = FALSE, add = TRUE);  text(text_C[1], text_C[2], "C", cex=1.7)
map('state', region = CR_NE, lwd=3, interior = FALSE, add = TRUE); text(text_NE[1], text_NE[2], "NE", cex=1.7)
map('state', region = CR_SE, lwd=3, interior = FALSE, add = TRUE); text(text_SE[1]-2, text_SE[2], "SE", cex=1.7)
title(main="Climatic Regions", cex.main=2)



##########################################################
###----- Now plotting data of north east region only #####
##########################################################
NE_201207_total<-total_201207[total_201207$CR=="NE",]
NE_201207_total<-total_201207[total_201207$CR=="NE",]

par(mfrow=c(3,4))
quilt.plot(x=NE_201207_total$Longitude,y=NE_201207_total$Latitude,z=NE_201207_total$pm25_monthly_mean,main="PM2.5")
quilt.plot(x=NE_201207_total$Longitude,y=NE_201207_total$Latitude,z=NE_201207_total$WS1,main="U-component (ws)")
quilt.plot(x=NE_201207_total$Longitude,y=NE_201207_total$Latitude,z=NE_201207_total$WS2,main="V-Component (ws)")
quilt.plot(x=NE_201207_total$Longitude,y=NE_201207_total$Latitude,z=NE_201207_total$WS,main="Magnitude Windspeed")
quilt.plot(x=NE_201207_total$Longitude,y=NE_201207_total$Latitude,z=NE_201207_total$GPH,main="Geopotential Height")
quilt.plot(x=NE_201207_total$Longitude,y=NE_201207_total$Latitude,z=NE_201207_total$TC,main="Total Cloud")
quilt.plot(x=NE_201207_total$Longitude,y=NE_201207_total$Latitude,z=NE_201207_total$SrfTemp,main="Surface Temperature")
quilt.plot(x=NE_201207_total$Longitude,y=NE_201207_total$Latitude,z=NE_201207_total$Temp_prs1000,main="Temperature at pressure 1000")
quilt.plot(x=NE_201207_total$Longitude,y=NE_201207_total$Latitude,z=NE_201207_total$TPrp,main="Total Precipitation")
quilt.plot(x=NE_201207_total$Longitude,y=NE_201207_total$Latitude,z=NE_201207_total$TPrpN,main="Total Precipitation (Nearest Grid Point)")
quilt.plot(x=NE_201207_total$Longitude,y=NE_201207_total$Latitude,z=NE_201207_total$RH2m,main="Relative Humidity (2m)")
quilt.plot(x=NE_201207_total$Longitude,y=NE_201207_total$Latitude,z=NE_201207_total$RHhyb,main="Relative Humidity hybrid")
quilt.plot(x=NE_201207_total$Longitude,y=NE_201207_total$Latitude,z=NE_201207_total$RHiso,main="Relative Humidity isotherm")
quilt.plot(x=NE_201207_total$Longitude,y=NE_201207_total$Latitude,z=NE_201207_total$AirTemp,main="Air Temperature")



par(mfrow=c(1,1))
map('state')
map('state', region = CR_NE, lwd=3, interior = FALSE, add = TRUE); 
points(NE_201207_total$Longitude,NE_201207_total$Latitude,col=color.scale(NE_201207_total$pm25_monthly_mean,tim.colors(),zlim=c(min(NE_201207_total$pm25_monthly_mean)-0.01,max(NE_201207_total$pm25_monthly_mean)+0.01)),pch=19,cex=0.5)
image.plot(legend.only = T,horizontal = T,zlim=c(min(NE_201207_total$pm25_monthly_mean)-0.01,max(NE_201207_total$pm25_monthly_mean)+0.01))
title("PM25 Data")
map('usa')
map('state', region = c("maine", "new hampshire", "vermont", "new york", "massachusetts", "connecticut", "rhode island",
                        "pennsylvania", "new jersey", "delaware", "maryland"))
map('state')
map('state', region = CR_NE, lwd=3, interior = FALSE, add = TRUE); text(text_NE[1], text_NE[2], "NE", cex=1.7)
points(NE_201207_total$Longitude,NE_201207_total$Latitude,col=color.scale(NE_201207_total$WS1,tim.colors(),zlim=c(min(NE_201207_total$WS1)-0.01,max(NE_201207_total$WS1)+0.01)),pch=19,cex=0.5)
image.plot(legend.only = T,horizontal = T,zlim=c(min(NE_201207_total$WS1)-0.01,max(NE_201207_total$WS1)+0.01))
title("Wind speed 1")

map('state', region = CR_NE, lwd=3, interior = FALSE, add = TRUE); text(text_NE[1], text_NE[2], "NE", cex=1.7)
points(NE_201207_total$Longitude,NE_201207_total$Latitude,col=color.scale(NE_201207_total$WS2,tim.colors(),zlim=c(min(NE_201207_total$WS2)-0.01,max(NE_201207_total$WS2)+0.01)),pch=19,cex=0.5)
image.plot(legend.only = T,horizontal = T,zlim=c(min(NE_201207_total$WS2)-0.01,max(NE_201207_total$WS2)+0.01))
title("Wind speed 2")

rm.list<-ls()[ls()!="NE_201207_total"]
rm(list=rm.list)
my_work_data_july2012<-NE_201207_total[,-c(1,19)]
write.csv(NE_201207_total, file = "my_data_NE_201207_totalv2.csv", row.names = FALSE)
write.csv(my_work_data_july2012, file = "my_data_my_work_data_july2012v2.csv", row.names = FALSE)
my_work_data_july2012<-read.csv("my_data_my_work_data_july2012v2.csv")
#rm(rm.list,NE_201206_total,my_work_data)

par(mfrow=c(3,4))
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=my_work_data_july2012$pm25_monthly_mean,main="PM2.5")
#quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=my_work_data_july2012$WS1,main="U-component (ws)")
#quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=my_work_data_july2012$WS2,main="V-Component (ws)")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=my_work_data_july2012$WS,main="Magnitude Windspeed")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=my_work_data_july2012$GPH,main="Geopotential Height")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=my_work_data_july2012$TC,main="Total Cloud")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=my_work_data_july2012$SrfTemp,main="Surface Temperature")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=my_work_data_july2012$Temp_prs1000,main="Temperature at pressure 1000")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=my_work_data_july2012$TPrp,main="Total Precipitation")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=my_work_data_july2012$TPrpN,main="Total Precipitation (Nearest Grid Point)")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=my_work_data_july2012$RH2m,main="Relative Humidity (2m)")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=my_work_data_july2012$RHhyb,main="Relative Humidity hybrid")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=my_work_data_july2012$RHiso,main="Relative Humidity isotherm")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=my_work_data_july2012$AirTemp,main="Relative Humidity isotherm")

############################################################
######### Plotting data plots for the manuscript ###########
############################################################
library(viridis)
par(cex.axis=3, cex.lab=4, cex.main=1, cex.sub=1,mar=c(3,2.5,0,1)+.1)
par(mfrow=c(1,1))
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=my_work_data_july2012$pm25_monthly_mean,xlab="Longitude",ylab="Latitude")
map('state',region=c("maine", "new hampshire", "vermont", "new york", "massachusetts", "connecticut", "rhode island",
                     "pennsylvania", "new jersey", "delaware", "maryland"),add = T,lwd=2)


quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=my_work_data_july2012$RH2m,xlab="Longitude",ylab="Latitude")
map('state',region=c("maine", "new hampshire", "vermont", "new york", "massachusetts", "connecticut", "rhode island",
                     "pennsylvania", "new jersey", "delaware", "maryland"),add = T,lwd=2)

quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=my_work_data_july2012$WS,xlab="Longitude",ylab="Latitude")
map('state',region=c("maine", "new hampshire", "vermont", "new york", "massachusetts", "connecticut", "rhode island",
                     "pennsylvania", "new jersey", "delaware", "maryland"),add = T,lwd=2)



###### Computing residuals ####
PM25<-my_work_data_july2012$pm25_monthly_mean
WS<-my_work_data_july2012$WS
GPH<-my_work_data_july2012$GPH
TC<-my_work_data_july2012$TC
SrfTemp<-my_work_data_july2012$SrfTemp
Temp_prs1000<-my_work_data_july2012$Temp_prs1000
TPrp<-my_work_data_july2012$TPrp
TPrpN<-my_work_data_july2012$TPrpN
RH2m<-my_work_data_july2012$RH2m
RHhyb<-my_work_data_july2012$RHhyb
RHiso<-my_work_data_july2012$RHiso
AirTemp<-my_work_data_july2012$AirTemp

PM25reg<-lm(PM25~1+my_work_data_july2012$Longitude+my_work_data_july2012$Latitude)
WSreg<-lm(WS~1+my_work_data_july2012$Longitude+my_work_data_july2012$Latitude)
GPHreg<-lm(GPH~1+my_work_data_july2012$Longitude+my_work_data_july2012$Latitude)
TCreg<-lm(TC~1+my_work_data_july2012$Longitude+my_work_data_july2012$Latitude)
SrfTempreg<-lm(SrfTemp~1+my_work_data_july2012$Longitude+my_work_data_july2012$Latitude)
Temp_prs1000reg<-lm(Temp_prs1000~1+my_work_data_july2012$Longitude+my_work_data_july2012$Latitude)
TPrpreg<-lm(TPrp~1+my_work_data_july2012$Longitude+my_work_data_july2012$Latitude)
TPrpNreg<-lm(TPrpN~1+my_work_data_july2012$Longitude+my_work_data_july2012$Latitude)
RH2mreg<-lm(RH2m~1+my_work_data_july2012$Longitude+my_work_data_july2012$Latitude)
RHhybreg<-lm(RHhyb~1+my_work_data_july2012$Longitude+my_work_data_july2012$Latitude)
RHisoreg<-lm(RHiso~1+my_work_data_july2012$Longitude+my_work_data_july2012$Latitude)
AirTempreg<-lm(AirTemp~1+my_work_data_july2012$Longitude+my_work_data_july2012$Latitude)

par(mfrow=c(3,4))
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=PM25reg$residuals,main="PM2.5 (Residuals)")
#quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=my_work_data_july2012$WS1,main="U-component (ws)")
#quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=my_work_data_july2012$WS2,main="V-Component (ws)")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=WSreg$residuals,main="Magnitude Windspeed (Residuals)")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=GPHreg$residuals,main="Geopotential Height (Residuals)")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=TCreg$residuals,main="Total Cloud (Residuals)")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=SrfTempreg$residuals,main="Surface Temperature (Residuals)")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=Temp_prs1000reg$residuals,main="Temperature at pressure 1000 (Residuals)")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=TPrpreg$residuals,main="Total Precipitation (Residuals)")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=TPrpNreg$residuals,main="Total Precipitation (Nearest Grid Point) (Residuals)")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=RH2mreg$residuals,main="Relative Humidity (2m) (Residuals)")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=RHhybreg$residuals,main="Relative Humidity hybrid (Residuals)")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=RHisoreg$residuals,main="Relative Humidity isotherm (Residuals)")
quilt.plot(x=my_work_data_july2012$Longitude,y=my_work_data_july2012$Latitude,z=AirTempreg$residuals,main="Air Temp (Residuals)")

cor(PM25reg$residuals,RHhybreg$residuals)
cor(PM25reg$residuals,RHisoreg$residuals)
cor(WSreg$residuals,RHhybreg$residuals)
cor(WSreg$residuals,RHisoreg$residuals)



cor(PM25reg$residuals,Temp_prs1000reg$residuals)
cor(PM25reg$residuals,RHisoreg$residuals)
cor(WSreg$residuals,Temp_prs1000reg$residuals)
cor(WSreg$residuals,RHisoreg$residuals)
cor(WSreg$residuals,TCreg$residuals)
cor(PM25reg$residuals,TCreg$residuals)
cor(PM25reg$residuals,GPHreg$residuals)
cor(PM25reg$residuals,Temp_prs1000reg$residuals,TCreg$residuals)
res.data<-data.frame(PM25res=PM25reg$residuals,
                     WSres=WSreg$residuals,
                     GPHres=GPHreg$residuals,
                     TCres=TCreg$residuals,
                     SrfTempres=SrfTempreg$residuals,
                     Temp_prs1000res=Temp_prs1000reg$residuals,
                     TPrpres=TPrpreg$residuals,
                     TPrpNres=TPrpNreg$residuals,
                     RH2mres=RH2mreg$residuals,
                     RHhybres=RHhybreg$residuals,
                     RHisores=RHisoreg$residuals,
                     AirTempres=AirTempreg$residuals)
cor(res.data)
