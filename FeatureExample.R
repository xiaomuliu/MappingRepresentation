setwd("/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/ExploratoryAnalysis/Mapping/")
load("MatchedBurglaryData_portal.RData")

library(sp)
library(rgeos)
library(rgdal)
library(raster)

Path.GIS <- "/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/GISData/"
Path.city <- paste0(Path.GIS,"City_Boundary")
city.shp <- readOGR(Path.city,"City_Boundary") 
cellsize.x <- 660
cellsize.y <- 660
X_range <- city.shp@bbox[1,]
Y_range <- city.shp@bbox[2,]
grd.full <- expand.grid(list(X_COORD=seq(X_range[1],X_range[2],by=cellsize.x),
                             Y_COORD=seq(Y_range[1],Y_range[2],by=cellsize.y)))
coordinates(grd.full) = ~X_COORD+Y_COORD # convert to SpatialPoints

prj <- paste("+proj=tmerc +lat_0=36.66666666666666 +lon_0=-88.33333333333333 +k=0.9999749999999999",
             "+x_0=300000 +y_0=0 +datum=NAD83 +units=us-ft +no_defs +ellps=GRS80 +towgs84=0,0,0")
proj4string(grd.full) <- prj

# rasterize the city spatial polygon to get a grid template
grd.full <- SpatialPixels(grd.full)
r <- raster(ncol=grd.full@grid@cells.dim[1],nrow=grd.full@grid@cells.dim[2],
            xmn=grd.full@bbox[1,1],xmx=grd.full@bbox[1,2],ymn=grd.full@bbox[2,1],ymx=grd.full@bbox[2,2],crs=CRS(prj))

city.raster <- rasterize(city.shp,r,0)
city.df_full <- as.data.frame(city.raster,xy=TRUE)
city.df_full <- city.df_full[,1:2]
names(city.df_full) <- c("X_COORD","Y_COORD")
RegGrd.full <- city.df_full

coordinates(city.df_full) <- c("X_COORD", "Y_COORD") 
proj4string(city.df_full) <- prj
BoundedOverFullGrd <- over(city.df_full, city.shp)
isInCity <- !is.na(BoundedOverFullGrd$OBJECTID)
RegGrd <- RegGrd.full[isInCity,]

load("SpatialFeatureMaps.RData")

isInZone <- NA
# ResiZone <- c(2,4,7,8,9)
# isInZone <- RegGrd$ZoneType %in% ResiZone
# RegGrd <- RegGrd[isInZone,]


RegGrdRaster.BldgDen <- rasterize(RegGrd[,c("X_COORD","Y_COORD")], r, RegGrd$BldgDen, fun=sum)
RegGrdRaster.StrDen <- rasterize(RegGrd[,c("X_COORD","Y_COORD")], r, RegGrd$StrDen, fun=sum)
RegGrdRaster.Dist2Street <- rasterize(RegGrd[,c("X_COORD","Y_COORD")], r, RegGrd$Dist2Street, fun=sum)

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
par(mfrow=c(1,2),mar=c(1,2,1,4),xaxs="i",yaxs="i",pty="s")
plot(RegGrdRaster.StrDen,col=jet.colors(256), xaxt="n", yaxt="n",xlab="",ylab="")
# plot(RegGrdRaster.Dist2Street,col=jet.colors(256),xaxt="n", yaxt="n",xlab="",ylab="")
plot(RegGrdRaster.BldgDen,col=jet.colors(256), xaxt="n", yaxt="n",xlab="",ylab="")