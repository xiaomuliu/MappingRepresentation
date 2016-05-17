setwd("/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/ExploratoryAnalysis/Mapping/")
load("MatchedBurglaryData_portal.RData")

CrimeData <- BurglaryData
rm(BurglaryData)

CrimeData.sub <- subset(CrimeData,DATEOCC>=as.Date("2014-12-01") & DATEOCC<=as.Date("2014-12-31"))

library(ggplot2)
library(ggmap)
library(sp)
#library(raster)
library(plotKML)
library(rgeos)
library(rgdal)
library(maptools)

# crime.sp <- CrimeData.sub
# coordinates(crime.sp) <- ~X_COORD+Y_COORD
# prj <- paste("+proj=tmerc +lat_0=36.66666666666666 +lon_0=-88.33333333333333 +k=0.9999749999999999",
#              "+x_0=300000 +y_0=0 +datum=NAD83 +units=us-ft +no_defs +ellps=GRS80 +towgs84=0,0,0")
# proj4string(crime.sp) <- prj
# crime <- spTransform(crime.sp, CRS("+init=epsg:4326")) # reproject to WGS84 (LONG LAT)
# Crime <- fortify(crime) 

Path.census_tract <- "/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/GISData/Census_Tracts"
census_tract.shp <- readOGR(Path.census_tract,"Census_Tracts")
census_tract <- spTransform(census_tract.shp, CRS("+init=epsg:4326")) # reproject to WGS84 (LONG LAT)
censusTract <- fortify(census_tract,region="OBJECTID") 

Path.CTA_routes <- "/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/GISData/CTA_Routes"
CTA_routes.shp <- readOGR(Path.CTA_routes,"CTA_Routes")
CTA_routes <- spTransform(CTA_routes.shp, CRS("+init=epsg:4326")) # reproject to WGS84 (LONG LAT)
CTAroute <- fortify(CTA_routes,region="OBJECTID") #simplify to use with ggplot

# Path.park <- "/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/GISData/Parks_Aug2012"
# park.shp <- readOGR(Path.park,"Parks_Aug2012")
# parks <- spTransform(park.shp, CRS("+init=epsg:4326")) # reproject to WGS84 (LONG LAT)
# park <- fortify(parks,region="PARK_NO") #simplify to use with ggplot

Path.school <- "/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/GISData/School_Grounds"
school.shp <- readOGR(Path.school,"School_Grounds")
schools <- spTransform(school.shp, CRS("+init=epsg:4326")) # reproject to WGS84 (LONG LAT)
school <- fortify(schools,region="OBJECTID") #simplify to use with ggplot

ChiMap <- qmap(location='chicago', zoom=11, maptype='roadmap',source="google")

ChiMap+geom_path(aes(x=long, y=lat,group=group), data=CTAroute, color='blue', size=1, lineend="butt", linejoin="round", linemitre=1)+
  geom_polygon(aes(x=long, y=lat,group=group), data=censusTract, color='white', alpha=.5, size=.6)+
  geom_polygon(aes(x=long, y=lat,group=group), data=school, color='green',fill='green', alpha=1, size=.6)+
  geom_point(data=CrimeData.sub, aes(x=LONG, y=LAT), color="red", size=1, alpha=1)
