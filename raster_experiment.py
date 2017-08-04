#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
===================================
Explore raster layers and rasterizing data
===================================
Created on Tue Aug 30 12:25:25 2016

@author: xiaomuliu
"""

from osgeo import gdal, ogr, osr
import numpy as np

def ogr2raster(ogr_in, raster_out, nrows=None, ncols=None, pixel_size=None, NoDataVal=-9999, burnVal=[0]):
    """
    Take in a OGR file (e.g. shapefile) and create a new raster Tiff file based on the shapefile
    pixel_size: a tuple of (pixelHeight, pixelWidth)
    burnVal: default value when burning into raster layers
    """
         
    # Open the data source and read in the extent
    source_ds = ogr.Open(ogr_in)
    source_layer = source_ds.GetLayer()
    x_min, x_max, y_min, y_max = source_layer.GetExtent()
   
    # Define pixel_size and NoData value of new raster
    try: 
        if pixel_size is not None:
            pixelHeight, pixelWidth = pixel_size
            nrows = int(np.ceil((y_max-y_min)/float(pixelHeight)))
            ncols = int(np.ceil((x_max-x_min)/float(pixelWidth)))
        elif (nrows is not None) and (ncols is not None):
            pixelWidth = (x_max-x_min)/float(ncols)
            pixelHeight = (y_max-y_min)/float(nrows)
        else:
            raise ValueError('Either the raster dimensionn (rows and ncols) or the raster resolution (pixel_size) must be specified')
    except ValueError as errmsg:
        print errmsg  
    
    # Create the destination data source
    target_ds = gdal.GetDriverByName('GTiff').Create(raster_out, ncols, nrows, 1, gdal.GDT_Byte)
    target_ds.SetGeoTransform((x_min, pixelWidth, 0, y_max, 0, -pixelHeight))
    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(NoDataVal)
   
    # Set projection
    spatialRef = source_layer.GetSpatialRef()
    target_ds.SetProjection(spatialRef.ExportToProj4())
    
    # Rasterize
    gdal.RasterizeLayer(target_ds, [1], source_layer, burn_values=burnVal);

    # Flush data to disk
    band.FlushCache()
    
    # Close files
    source_ds = None
    target_ds = None
    return {'raster':target_ds,'band':target_ds.GetRasterBand(1)}


def raster2array(raster_in, xOffset=0, yOffset=0, xSize=None, ySize=None):
    """
    Get the (single) Band object by passing the band index (1-based) to the Dataset’s
    Read the data into a 2D Numeric array with ReadAsArray
    NOTE: For raster indexing, read individual pixels using [yoff, xoff]
          (math matrix notation is [row,col], not [x,y])
    """
    raster = gdal.Open(raster_in)
    band = raster.GetRasterBand(1)
    if xSize is None:
        xSize = raster.RasterXSize # ncols
    if ySize is None:
        ySize = raster.RasterYSize # nrows

    array = band.ReadAsArray(xOffset, yOffset, xSize, ySize)
    return array



def array2raster(raster_out, array, raster_in=None, layerExtent=(0,0,1,1), xOffset=0, yOffset=0, NoDataVal=None):
    """
    Write the numeric array into a single-band raster object (raster_out). If an existing raster object (raster_in) 
    is provided then its layer and spatial information will be used to create a new raster object.
    Otherwise the raster layer will be created according to layerExtent (a tuple in form of (xmin,ymin,xmax,ymax))
    """
    nrows,ncols = np.shape(array)
    if raster_in is None:
        x_min,y_min,x_max,y_max = layerExtent                       
        pixelWidth = (x_max-x_min)/float(ncols)
        pixelHeight = (y_max-y_min)/float(nrows)
        geotransform = (x_min, pixelWidth, 0, y_max, 0, -pixelHeight)
    else:  
        #**********************************************************#
        #adfGeoTransform[0] /* top left x */
        #adfGeoTransform[1] /* west-east pixel resolution */
        #adfGeoTransform[2] /* rotation, 0 if image is "north up" */
        #adfGeoTransform[3] /* top left y */
        #adfGeoTransform[4] /* rotation, 0 if image is "north up" */
        #adfGeoTransform[5] /* north-south pixel resolution */
        #**********************************************************#
        raster = gdal.Open(raster_in)
        geotransform = raster.GetGeoTransform()
        x_min = geotransform[0]
        y_max = geotransform[3]
        pixelWidth = geotransform[1]
        pixelHeight = geotransform[5]

    driver = gdal.GetDriverByName('GTiff')
    # Create(<filename>, <xsize>, <ysize>,[<bands>], [<GDALDataType>])
    outRaster = driver.Create(raster_out, ncols, nrows, 1, gdal.GDT_Byte)
    outRaster.SetGeoTransform((x_min, pixelWidth, 0, y_max, 0, -pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array, xOffset, yOffset)
    if NoDataVal is None:
        NoDataVal = raster.GetRasterBand(1).GetNoDataValue()
        if NoDataVal is None:
            # rasterfn does not have No Data value
            outband.SetNoDataValue(-9999)
    else:
       outband.SetNoDataValue(NoDataVal)
    
       
    # Georeference the new raster
    outRaster.SetGeoTransform(geotransform)
        
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromProj4(raster.GetProjectionRef())
    outRaster.SetProjection(outRasterSRS.ExportToProj4())
    
    # Flush data to disk
    outband.FlushCache()
    
    # close raster file
    outRaster = None
    return {'raster':outRaster,'band':outRaster.GetRasterBand(1)}
    
        
def pixelOffset2coord(raster, x_index=None, y_index=None, center=True):
    """
    Get X Y coordinate(s) indexing by x_index and y_index (can be scalar or vector).
    if x_index and y_index are None, return the coordinates of all the pixels
    If center is true, return the center coordinates of pixel cells,
    otherwise the top-left coordinates of pixel cells
    """
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]

    if x_index is None:
        ncols = raster.RasterXSize
        x_index = np.arange(ncols)
    if y_index is None:
        nrows = raster.RasterYSize
        y_index = np.arange(nrows) 
    # convert x_index and y_index to numpy arrays to enable vector operations
    x_index = np.array(x_index).astype(int)
    y_index = np.array(y_index).astype(int)
    if center==True:
       coordX = originX + pixelWidth*x_index + pixelWidth/float(2)
       coordY = originY + pixelHeight*y_index + pixelHeight/float(2)
    else:
       coordX = originX + pixelWidth*x_index
       coordY = originY + pixelHeight*y_index
    return coordX, coordY 

    
if __name__ == '__main__':   
    import matplotlib.pyplot as plt
    
    # Register all of the drivers
    gdal.AllRegister()
    
    #*****************************************************#
    # Take in a OGR file (e.g. shapefile) and create a new raster Tiff file based on the shapefile
    
    # Define pixel_size and NoData value of new raster
    pixel_size = (660,660)
    NoData_value = -9999
    
    # Filename of input OGR file
    shpfile_city = "/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/GISdata/City_Boundary/City_Boundary.shp"
    
    # Filename of the raster data that will be created
    raster_fn = 'raster_test.tif'
    
    # Open the data source and read in the extent
    source_ds = ogr.Open(shpfile_city)
    source_layer = source_ds.GetLayer()
    x_min, x_max, y_min, y_max = source_layer.GetExtent()
    
    spatialRef = source_layer.GetSpatialRef()
    spatialRef.ExportToProj4()
    # The above proj4. string should be similar to this one
    proj_str = "+proj=tmerc +lat_0=36.66666666666666 +lon_0=-88.33333333333333 +k=0.9999749999999999 " + \
              "+x_0=300000 +y_0=0 +datum=NAD83 +units=us-ft +no_defs +ellps=GRS80 +towgs84=0,0,0"   
    proj = osr.SpatialReference()
    proj.ImportFromProj4(proj_str)
    proj.ExportToProj4()
             
    # Create the destination data source
    # NOTE:
    # *Python 2.x floating point division
    # In order to be consistent with the mesh grid, the number of columns and rows are calculated as below
    # (Compare with the verification code below)  
    cols = int(np.ceil((int(x_max) - int(x_min)) / float(pixel_size[0])))
    rows = int(np.ceil((int(y_max) - int(y_min)) / float(pixel_size[1])))
    
    outData = np.random.uniform(low=0,high=100,size=(rows,cols)).astype(np.int32)
    
    target_ds = gdal.GetDriverByName('GTiff').Create(raster_fn, cols, rows, 1, gdal.GDT_Byte)
    target_ds.SetGeoTransform((x_min, pixel_size[0], 0, y_max, 0, -pixel_size[1]));
    #band = target_ds.GetRasterBand(1)
    #band.WriteArray(outData);
    #band.SetNoDataValue(NoData_value);
    #target_ds.SetProjection('LOCAL_CS["arbitrary"]');
    #band.FlushCache()
    target_ds.GetRasterBand(1).WriteArray(outData);
    # target_ds.SetProjection('LOCAL_CS["arbitrary"]');
    target_ds.SetProjection(proj_str)
    target_ds.FlushCache()
    # Close raster file
    # raster files must be closed otherwise the image values on disk are all zeros
    target_ds = None
        
    # Rasterize
    #gdal.RasterizeLayer(target_ds, [1], source_layer, burn_values=[0]);
    
    # plot
    bandArray = raster2array(raster_fn)
    plt.imshow(bandArray) 
    
    # Test functions
    raster_file1 = 'raster_test1.tif'
    raster_file2 = 'raster_test2.tif'
    outData = np.reshape(np.arange(rows*cols),(rows,cols))
    outData = (outData-np.min(outData))/float(np.max(outData)-np.min(outData))*255
    ogr2raster(ogr_in=shpfile_city, raster_out=raster_file1, pixel_size=pixel_size, NoDataVal=NoData_value)
    array2raster(raster_out=raster_file2, array=outData, raster_in=raster_file1, NoDataVal=NoData_value)
    bandArray = raster2array(raster_file2)
    plt.imshow(bandArray) 
    
    raster = gdal.Open(raster_file2)
    X_coord, Y_coord = pixelOffset2coord(raster,center=True)
  
    # Verify with mesh grid
    import geopandas as gpd
    city_shp = gpd.read_file(shpfile_city)
    range_x = city_shp.bounds[['minx','maxx']].values[0].astype(int)
    range_y = city_shp.bounds[['miny','maxy']].values[0].astype(int)    
    grd_x = np.arange(range_x[0], range_x[1], step=pixel_size[0], dtype=int)
    grd_y = np.arange(range_y[0], range_y[1], step=pixel_size[1], dtype=int)
    grd_rect = np.meshgrid(grd_x,grd_y)
    # NOTE: 
    # 1. The length of the np.arange result is ceil((stop - start)/step). 
    # Because of floating point overflow, this rule may result in the last element of out being greater than stop.
    # 2. Mesh grid gives the top-left coordinates of pixel cells
    # while pixel2Offset2coord (if center==True) gives the center coordinates of pixel cells
    # 3. Raster cells is in this order
    # Y direction: top (north) to bottom (south); X direction: from left (west) to right (east)
    # while mesh grid is in this order
    # Y direction: bottom (south) to top (north); X direction: from left (west) to right (east)