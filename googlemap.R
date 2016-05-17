setwd("/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/ExploratoryAnalysis/Mapping/")
load("MatchedBurglaryData_portal.RData")

CrimeData <- BurglaryData
rm(BurglaryData)

CrimeData.sub <- subset(CrimeData,DATEOCC>=as.Date("2014-12-01") & DATEOCC<=as.Date("2014-12-31"))


library(ggplot2)
library(ggmap)

#plot the  hybrid Google Maps basemap
map <- qmap(location='chicago', zoom = 12, maptype = 'hybrid')
#plot the crime points on top
map + geom_point(data = CrimeData.sub, aes(x = LONG, y = LAT), color="red", size=2, alpha=0.8)


map + stat_density2d(aes(x = LONG, y = LAT, fill = ..level.., alpha = ..level..*2),
                     size = 1, bins = 10, data = CrimeData.sub, geom = "polygon",h=c(0.03,0.03),show.legend =FALSE) +
  scale_fill_gradient(low = "black", high = "red")
