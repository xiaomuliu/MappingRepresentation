setwd("/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/ExploratoryAnalysis/Mapping/")
load("MatchedBurglaryData_portal.RData")

library(sp)
library(rgeos)
library(rgdal)
library(raster)

Path.GIS <- "/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/GISData/"
Path.city <- paste0(Path.GIS,"City_Boundary")
city.shp <- readOGR(Path.city,"City_Boundary") 
cellsize.x <- 1000
cellsize.y <- 1000
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

month <- 1
year <- 2013

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

# Point plot
CrimeSubset <- subset(BurglaryData,YEAR==year & MONTH==month,select = c("DATEOCC","BEAT","DISTRICT","X_COORD","Y_COORD","INC_CNT"))
CrimeSubset.sp <- SpatialPoints(coords=CrimeSubset[,c("X_COORD","Y_COORD")])

par(mar=c(2,1,1,1),xaxs="i",yaxs="i",cex.axis=0.8,cex.lab=0.8,pty="s")
plot(city.shp, border="black")
points(CrimeSubset.sp, pch=16, cex=.75,col="red")
box(which = "plot", lty = "solid")

par(mar=c(2,1,1,1),xaxs="i",yaxs="i",cex.axis=0.8,cex.lab=0.8,pty="s")
plot(city.shp, border="black")
points(CrimeSubset.sp, pch=16, cex=.75,col="red")
box(which = "plot", lty = "solid")


# Police beat plot
shapefilePath.new <- "/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/CPDShapeFiles/new/"
beat_new.rg <- readOGR(paste0(shapefilePath.new,"cpd_beats"), "cpd_beats")

CrimeBeat <- aggregate(INC_CNT~BEAT,data=CrimeSubset, FUN=sum, na.rm=TRUE)

Beat_template.spdf <- beat_new.rg
# remove some useless/redundant attributes
Beat_template.spdf@data$DISTRICT <- NULL
Beat_template.spdf@data$SECTOR <- NULL
Beat_template.spdf@data$BEAT <- NULL
Beat_template.spdf@data$BEAT_NUM <- NULL
# add an attribute INC_CNT
Beat_template.spdf@data$INC_CNT <- rep(0,nrow(Beat_template.spdf@data))

CrimeBeat.spdf <- Beat_template.spdf
BeatList <- sort(unique(BurglaryData$BEAT))
BeatList.df <- data.frame(BEAT=BeatList)
CrimeBeat <- merge(CrimeBeat,BeatList.df,all=TRUE)
CrimeBeat$INC_CNT[is.na(CrimeBeat$INC_CNT)] <- 0
for (j in BeatList){
  if (j == "3100"){next}
  CrimeBeat.spdf@data$INC_CNT[CrimeBeat.spdf@data$BEAT_NUMBE==j] <- CrimeBeat$INC_CNT[CrimeBeat$BEAT==j]
}

library(lattice)
# axis.line <- trellis.par.get("axis.line") 
# trellis.par.set(axis.line=list(col=NA)) 
# f <- spplot(CrimeBeat.spdf, zcol="INC_CNT", col.regions=jet.colors(256),colorkey=list(height=0.5,width=1),
#             par.settings =list(axis.line=list(col='transparent')))
f <- spplot(CrimeBeat.spdf, zcol="INC_CNT", col.regions=jet.colors(256),colorkey=list(height=0.5,width=1))
print(f)

Nx <- round(grd.full@grid@cells.dim[1]*3)
Ny <- round(grd.full@grid@cells.dim[2]*3)
r3 <- raster(ncol=Nx,nrow=Ny,xmn=grd.full@bbox[1,1],xmx=grd.full@bbox[1,2],ymn=grd.full@bbox[2,1],ymx=grd.full@bbox[2,2],crs=CRS(prj))
CrimeBeat.raster <- rasterize(CrimeBeat.spdf,r3,"INC_CNT",fun=sum)
plot(CrimeBeat.raster,col=jet.colors(256), xaxt="n", yaxt="n",xlab="",ylab="")
plot(beat_new.rg,add=TRUE)

# Grid plot
Nx <- round(grd.full@grid@cells.dim[1]/2.5)
Ny <- round(grd.full@grid@cells.dim[2]/2.5)
r2 <- raster(ncol=Nx,nrow=Ny,xmn=grd.full@bbox[1,1],xmx=grd.full@bbox[1,2],ymn=grd.full@bbox[2,1],ymx=grd.full@bbox[2,2],crs=CRS(prj))
CrimeSubset.raster <- rasterize(CrimeSubset[,c("X_COORD","Y_COORD")], r2, CrimeSubset$INC_CNT, fun=sum)
plot(CrimeSubset.raster, panel.first=grid(Nx,Ny,col="black", lty="dotted"),col=jet.colors(256), xaxt="n", yaxt="n",xlab="",ylab="")
plot(city.shp,add=TRUE)


# KDE  
source("ConvFunction.R")
sigma <- 2
window.x <- 3*sigma  # unit: raster cell size
window.y <- 3*sigma
kernel.Xdim <- 2*window.x+1
kernel.Ydim <- 2*window.y+1

kernel.Xgrd <- seq(-window.x,window.x,length.out=kernel.Xdim)
kernel.Ygrd <- seq(-window.y,window.y,length.out=kernel.Ydim)

kernel.x <- 1/(sqrt(2*pi)*sigma)*exp(-kernel.Xgrd^2/(2*sigma^2))
kernel.y <- 1/(sqrt(2*pi)*sigma)*exp(-kernel.Ygrd^2/(2*sigma^2))

KDE <- SpatialKernSmCrime(CrimeSubset,RegGrd.full,r,kernel.x,kernel.y,isInCity,prj)
KDE.df_inCity <- KDE$KernSm.df_inPoly
KDE.df_inCity$VALUE <- KDE.df_inCity$KS_VAL/sum(KDE.df_inCity$KS_VAL)

KDE.raster_inCity <- rasterize(KDE.df_inCity[,c("X_COORD","Y_COORD")], r, KDE.df_inCity$VALUE, fun=sum)
plot(KDE.raster_inCity,col=jet.colors(256), xaxt="n", yaxt="n",xlab="",ylab="")
plot(city.shp,add=TRUE)

# library(KernSmooth)
# h <- 900
# Nx <- grd.full@grid@cells.dim[1]
# Ny <- grd.full@grid@cells.dim[2]
# KDE.df <- expand.grid(list(X_COORD=seq(X_range[1],X_range[2],length.out=Nx),
#                            Y_COORD=seq(Y_range[1],Y_range[2],length.out=Ny)))
# kernSm <- bkde2D(data.matrix(CrimeSubset[,c("X_COORD","Y_COORD")]), bandwidth=c(h,h), 
#                    gridsize=c(Nx, Ny), range.x=list(X_range,Y_range))
#   
# KDE.df$VALUE <- as.vector(kernSm$fhat)     
# library(lattice)
# levelplot(VALUE~X_COORD*Y_COORD, data=KDE.df,col.regions=jet.colors(256),
#           xlab="X Coordinate", ylab="Y Coordinate", colorkey=list(width=0.5,height=0.8))