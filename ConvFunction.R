filterConv <- function(x,y){
  # Note that in R, the usual definition of convolution of two sequences x and y is given by convolve(x, rev(y), type = "open")
  # "filter" returns the middle sub-vector of "open". By zero-padding x, it can have the same effect of conv(x,y,'same') as in MATLAB
  return(zapsmall(convolve(c(rep(0,floor(length(y)/2)),x,rep(0,floor(length(y)/2))),rev(y),type="filter")))
}

truncConv <- function(x,y){
  # x is of length M, y is of length N, linear convolution returns a sequence z of length M+N-1
  # truncConv returns the first M samples of z
  z <- zapsmall(convolve(x,rev(y),type="open"))
  return(z[1:length(x)])
}

KernSmNearbyCrime <- function(CrimeData,QueryDay,Grid,raster,window.t,kernel.x,kernel.y,kernel.t,isInPoly,prj,isInZone=NA){
  Crime.NeighborArray <- array(data=NA,dim=c(raster@nrows,raster@ncols,window.t))
  NeighborArray_KS.xy <- array(data=NA,dim=c(raster@nrows,raster@ncols,window.t))
  KernSm.df_full <- Grid[,c("X_COORD","Y_COORD")]
  
  for (j in 1:window.t){
    CrimePts <- subset(CrimeData,DATEOCC==QueryDay-window.t+j-1,select=c("X_COORD","Y_COORD","INC_CNT"))
    if (nrow(CrimePts)==0){
      # no crime incident for a certain day
      Crime.NeighborArray[, , j] <- 0
    }else{
      CrimePts.raster <- rasterize(CrimePts[,c("X_COORD","Y_COORD")], raster, CrimePts$INC_CNT, fun=sum)
      Crime.NeighborArray[, , j] <- as.matrix(CrimePts.raster)
    }
  }
  Crime.NeighborArray[is.na(Crime.NeighborArray)] <- 0
  
  for (j in 1:dim(Crime.NeighborArray)[3]){
    NeighborArray_KS.x <- apply(Crime.NeighborArray[, , j],1,filterConv,kernel.x)
    NeighborArray_KS.y <- apply(NeighborArray_KS.x,1,filterConv,kernel.y)
    NeighborArray_KS.xy[, , j] <- NeighborArray_KS.y 
  }
  # We only need the last spatial layer. So use multiplication instead of convolution to speed up
  NeighborArray_KS.xyt <- apply(NeighborArray_KS.xy,c(1,2),FUN=function(x){x %*% rev(kernel.t)})
  
  KernSm.df_full$KS_VAL <- as.vector(t(NeighborArray_KS.xyt))
  
  KernSm.sp <- KernSm.df_full
  coordinates(KernSm.sp) <- c("X_COORD", "Y_COORD") 
  proj4string(KernSm.sp) <- prj
  KernSm.sp <- as(KernSm.sp,"SpatialPointsDataFrame")
  
  KernSm.sp_inPoly <- KernSm.sp
  KernSm.sp_inPoly@data <- data.frame(KS_VAL=KernSm.sp@data$KS_VAL[isInPoly])
  KernSm.sp_inPoly@coords <- KernSm.sp@coords[isInPoly,]
  
  KernSm.df_inPoly <- as.data.frame(KernSm.sp_inPoly)
  
  if (all(!is.na(isInZone))){
    KernSm.df_inPoly <- KernSm.df_inPoly[isInZone,]
    KernSm.sp_inPoly@data <- data.frame(KS_VAL=KernSm.sp@data$KS_VAL[isInZone])
    KernSm.sp_inPoly@coords <- KernSm.sp@coords[isInZone,]
  }
  
  return(list(KernSm.df_inPoly=KernSm.df_inPoly,KernSm.sp_inPoly=KernSm.sp_inPoly))
} 

KernSmNearbyCrime2 <- function(CrimeData,Query,Grid,raster,window.t,kernel.x,kernel.y,kernel.t,isInPoly,prj,isInZone=NA){
  Crime.NeighborArray <- array(data=NA,dim=c(raster@nrows,raster@ncols,window.t))
  NeighborArray_KS.xy <- array(data=NA,dim=c(raster@nrows,raster@ncols,window.t))
  KernSm.df_full <- Grid[,c("X_COORD","Y_COORD")]
  
  for (j in 1:window.t){
    CrimePts <- subset(CrimeData,GROUP==Query-window.t+j-1,select=c("X_COORD","Y_COORD","INC_CNT"))
    if (nrow(CrimePts)==0){
      # no crime incident for a certain day
      Crime.NeighborArray[, , j] <- 0
    }else{
      CrimePts.raster <- rasterize(CrimePts[,c("X_COORD","Y_COORD")], raster, CrimePts$INC_CNT, fun=sum)
      Crime.NeighborArray[, , j] <- as.matrix(CrimePts.raster)
    }
  }
  Crime.NeighborArray[is.na(Crime.NeighborArray)] <- 0
  
  for (j in 1:dim(Crime.NeighborArray)[3]){
    NeighborArray_KS.x <- apply(Crime.NeighborArray[, , j],1,filterConv,kernel.x)
    NeighborArray_KS.y <- apply(NeighborArray_KS.x,1,filterConv,kernel.y)
    NeighborArray_KS.xy[, , j] <- NeighborArray_KS.y 
  }
  # We only need the last spatial layer. So use multiplication instead of convolution to speed up
  NeighborArray_KS.xyt <- apply(NeighborArray_KS.xy,c(1,2),FUN=function(x){x %*% rev(kernel.t)})
  
  KernSm.df_full$KS_VAL <- as.vector(t(NeighborArray_KS.xyt))
  
  KernSm.sp <- KernSm.df_full
  coordinates(KernSm.sp) <- c("X_COORD", "Y_COORD") 
  proj4string(KernSm.sp) <- prj
  KernSm.sp <- as(KernSm.sp,"SpatialPointsDataFrame")
  
  KernSm.sp_inPoly <- KernSm.sp
  KernSm.sp_inPoly@data <- data.frame(KS_VAL=KernSm.sp@data$KS_VAL[isInPoly])
  KernSm.sp_inPoly@coords <- KernSm.sp@coords[isInPoly,]
  
  KernSm.df_inPoly <- as.data.frame(KernSm.sp_inPoly)
  
  if (all(!is.na(isInZone))){
    KernSm.df_inPoly <- KernSm.df_inPoly[isInZone,]
    KernSm.sp_inPoly@data <- data.frame(KS_VAL=KernSm.sp@data$KS_VAL[isInZone])
    KernSm.sp_inPoly@coords <- KernSm.sp@coords[isInZone,]
  }
  
  return(list(KernSm.df_inPoly=KernSm.df_inPoly,KernSm.sp_inPoly=KernSm.sp_inPoly))
} 

SpatialKernSmCrime <- function(CrimeData,Grid,raster,kernel.x,kernel.y,isInPoly,prj,isInZone=NA){
  Crime.NeighborArray <- matrix(NA,nrow=raster@nrows,ncol=raster@ncols)
  NeighborArray_KS.xy <- matrix(NA,nrow=raster@nrows,ncol=raster@ncols)
  KernSm.df_full <- Grid[,c("X_COORD","Y_COORD")]
  
  CrimePts <- subset(CrimeData,select=c("X_COORD","Y_COORD","INC_CNT"))
  CrimePts.raster <- rasterize(CrimePts[,c("X_COORD","Y_COORD")], raster, CrimePts$INC_CNT, fun=sum)
  Crime.NeighborArray <- as.matrix(CrimePts.raster)
  Crime.NeighborArray[is.na(Crime.NeighborArray)] <- 0
  
  NeighborArray_KS.x <- apply(Crime.NeighborArray,1,filterConv,kernel.x)
  NeighborArray_KS.xy <- apply(NeighborArray_KS.x,1,filterConv,kernel.y)
    
  KernSm.df_full$KS_VAL <- as.vector(t(NeighborArray_KS.xy))
  
  KernSm.sp <- KernSm.df_full
  coordinates(KernSm.sp) <- c("X_COORD", "Y_COORD") 
  proj4string(KernSm.sp) <- prj
  KernSm.sp <- as(KernSm.sp,"SpatialPointsDataFrame")
  
  KernSm.sp_inPoly <- KernSm.sp
  KernSm.sp_inPoly@data <- data.frame(KS_VAL=KernSm.sp@data$KS_VAL[isInPoly])
  KernSm.sp_inPoly@coords <- KernSm.sp@coords[isInPoly,]
  
  KernSm.df_inPoly <- as.data.frame(KernSm.sp_inPoly)
  
  if (all(!is.na(isInZone))){
    KernSm.df_inPoly <- KernSm.df_inPoly[isInZone,]
    KernSm.sp_inPoly@data <- data.frame(KS_VAL=KernSm.sp@data$KS_VAL[isInZone])
    KernSm.sp_inPoly@coords <- KernSm.sp@coords[isInZone,]
  }
  
  return(list(KernSm.df_inPoly=KernSm.df_inPoly,KernSm.sp_inPoly=KernSm.sp_inPoly))
} 

minmaxScale <- function(x,center=min(x),scale=max(x)-min(x)){
  x <- (x-center)/scale
  return(x)
}

MapCount <- function(Data,Attr,r,isInCity,isInZone=NA,fun=sum){
  Raster <- rasterize(Data[,c("X_COORD","Y_COORD")], r, Data[,Attr], fun=fun)
  Raster.df_inCity <- as.data.frame(Raster)[isInCity,]
  Raster.df_inCity[is.na(Raster.df_inCity)] <- 0
  if (!all(is.na(isInZone))){
    Raster.df_inCity <- Raster.df_inCity[isInZone]
  }
  return(Raster.df_inCity)
}