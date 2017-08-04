# -*- coding: utf-8 -*-
"""
==================================
Specify Spatial Grids
==================================
Created on Mon Aug 22 14:42:26 2016

@author: xiaomuliu
"""
      
import numpy as np
          
def vectorized_grid(shpfile,cellsize,center=True):
    cellsize_x, cellsize_y = cellsize
    range_x = shpfile.bounds[['minx','maxx']].values[0].astype(int)
    range_y = shpfile.bounds[['miny','maxy']].values[0].astype(int)
    # grids cover a rectangular area
    grd_x = np.arange(range_x[0], range_x[1], step=cellsize_x, dtype=int)
    grd_y = np.arange(range_y[0], range_y[1], step=cellsize_y, dtype=int)
    
    if center:
        # assume grid coordinates are the center of cells (for later calulating points in polygon)
        # shift by half cell size
        grd_x = grd_x+cellsize_x/2
        grd_y = grd_y+cellsize_x/2
    
    grd_rect = np.meshgrid(grd_x,grd_y)    
    # reshape the mesh grid of a vector of points of form [(x1,y1),(x2,y1),(x3,y1),...,(x1,y2),(x2,y2),...,(xn,ym)]
    grd_vec = np.dstack(grd_rect).reshape(len(grd_x)*len(grd_y),2)
    # equivalently,
    # grd_vec = np.vstack([grd_rect[0].ravel(), grd_rect[1].ravel()]).T
    
    return grd_vec, grd_x, grd_y

from shapely.geometry import Point  

def array_to_geoIm(array2d):
    """
    # ************************************************************ #
    # The input 2d  matrix are in the form of
    #         y_coord   181xxxx, ..., 195xxxx 
    # x_coord             
    # 109xxxx           val_00, ...,  val_0N
    # ...
    # 120xxxx           val_M0, ...,  val_MN
    #
    #
    # Return a re-arrange matrix with elements to be as the following
    #         x_coord   109xxxx, ..., 120xxxx 
    # y_coord             
    # 195xxxx           val_0N, ...,  val_MN
    # ...
    # 181xxxx           val_00, ...,  val_M0
    
    # ************************************************************ #
    """ 
    flipped_array = np.copy(np.fliplr(array2d).T) 
    return flipped_array


def grid_within_bndy(grid,shpfile,poly_index=0,im2d_mask=False,geoIm=False):
    """
    Returns grid cells within boundary which is defined by a polygon (indicated by poly_index) from the shpfile,
    as well as the corresponding mask.
    'grid' must be provided by a list with two vectors (or a 2d array) points coordinates
    (one for x coordinates, one for y coordinates)
    If im2d_mask=False, these coordinates can be of irregular grids and all grid coordinates must be given.
    If im2d_mask=True, the corrdinates can just be x-direction and y-direction coordinates.
    And the returned mask is in 2d (image) matrix form.
    X and y coordinates in grid are assumed to be in ascending order. 
    """
    if im2d_mask==True:
        grd_mesh = np.meshgrid(grid[0],grid[1])
        grid_vec = np.dstack(grd_mesh).reshape(len(grid[0])*len(grid[1]),2)
    else:
        grid_vec = grid
        
    mask_grdInPoly = [shpfile.geometry[poly_index].contains(Point(pt)) for pt in grid_vec]
    mask_grdInPoly = np.array(mask_grdInPoly)
    grd_inPoly = grid_vec[mask_grdInPoly,:]

    if im2d_mask==True:
        mask_grdInPoly = mask_grdInPoly.reshape((len(grid[0]),len(grid[1])),order='F') #column major 
        if geoIm==True:
            # please see the docstring of function 'array_to_geoIm' for the specifiction of x-y coordination layout
            mask_grdInPoly = array_to_geoIm(mask_grdInPoly)
            
    return {'mask':mask_grdInPoly, 'grid':grd_inPoly}


if __name__ == '__main__':
    ##############################################################################
    # Three equivalent ways of generating point-in-polyon mask
    #import time
    #from shapely.geometry import Point
    
    #start = time.time()
    #mask1 = [city_shp.geometry[0].contains(Point(pt)) for pt in grd_vec]
    #end = time.time()
    #print(end - start)
    #
    #start = time.time()
    #mask2 = [Point(pt).within(city_shp.geometry[0]) for pt in grd_vec]
    #end = time.time()
    #print(end - start)
    #
    # Pre-allocation doesn't improve performance
    #start = time.time()
    #mask3 = [False]*len(grd_vec)
    #for i, pt in enumerate(grd_vec):
    #    if city_shp.geometry[0].contains(Point(pt)):
    #        mask3[i] = True
    #end = time.time()
    #print(end - start)
    ##############################################################################

    import geopandas as gpd
    path_GIS = "/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/PredictiveMapping/Data/GISData/"
    shpfile_city = path_GIS + "City_Boundary/City_Boundary.shp"
    city_shp = gpd.read_file(shpfile_city)
        
    cellsize = (660,660)
    grd_vec, grd_x, grd_y = vectorized_grid(city_shp,cellsize)
    
    mask_grdInCity = grid_within_bndy(grd_vec,city_shp)['mask']
    grd_vec_inCity = grd_vec[mask_grdInCity,:]

    maskIm_grdInCity = grid_within_bndy([grd_x,grd_y],city_shp,im2d_mask=True,geoIm=True)['mask']
    
    import matplotlib.pyplot as plt
    
    fig = plt.figure()
    im = plt.imshow(maskIm_grdInCity, origin='upper',interpolation='none')


