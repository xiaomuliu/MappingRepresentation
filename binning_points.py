#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
===================================
Rasterize point data
===================================
Created on Thu Sep  1 16:45:51 2016

@author: xiaomuliu

One way of rasterizing points is through the process of creating raster and assigning points to the corresponding cells
However, GDAL does not provide handy tools for such task. As one has to do either of the followings
1. Convert point data to point shapefile with a field 'point count' and make it a spatial layer, then use gdal.RasterizeLayer to 
 rasterize this shapefile with point count field.
2. Create a raster then loop over all the raster cells and points to find their membership, which is inefficient.

This task can be simply done using 
1. numpy histogram and binning functionalities. Or
2. scipy binned statistics
Later the work of assigning numpy arrrays to rasters or even spatial data format such as shapefile is trivial.

"""

import numpy as np
from scipy.stats import binned_statistic_2d
import geopandas as gpd
from SetupGrid import vectorized_grid, grid_within_bndy

shpfile_city = "/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/GISdata/City_Boundary/City_Boundary.shp"
city_shp = gpd.read_file(shpfile_city)
    
cellsize = (660,660)
grd_vec, grd_x, grd_y = vectorized_grid(city_shp,cellsize)

mask_grdInCity = grid_within_bndy(grd_vec,city_shp)['mask']
grd_vec_inCity = grd_vec[mask_grdInCity,:]


# import crime (point) data
import cPickle as pickle

filePath = "/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/PredictiveMapping_py/CrimeData/"
pkl_file = open(filePath+'BURGLARY_01_14.pkl', 'rb')
CrimeData = pickle.load(pkl_file)
pkl_file.close()

x_coord = CrimeData['X_COORD'].values
y_coord = CrimeData['Y_COORD'].values

# bin
# Since grd_x grd_y correspond to the center of each grid cell, we add half cell size to each side of the grid coordinate to
# get the edge values
xedges = np.r_[grd_x-cellsize[0]/2, grd_x[-1]+cellsize[0]/2]
yedges = np.r_[grd_y-cellsize[1]/2, grd_y[-1]+cellsize[1]/2]

# counts (histogram)
Hist2d, xedges, yedges = np.histogram2d(x_coord, y_coord, bins=(xedges, yedges))
Hist2d_vec = Hist2d.reshape(-1,order='F') # Pay attention to 'C order' and 'Fortran order'
Hist2d_vec_inCity = Hist2d_vec[mask_grdInCity]

# scipy.stats.binned_statistic_2d gives a generalization of a histogram2d function allowing 
# the computation of the sum, mean, median, or other statistic of the values within each bin
Bin2d_stat, xedges2, yedges2, Bin_No = binned_statistic_2d(x_coord, y_coord, values=None, statistic='count', bins=(xedges, yedges))
np.all(Bin2d_stat==Hist2d)
np.all(xedges==xedges2) and np.all(yedges==yedges2)

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
    return np.fliplr(array2d).T 


Hist2d_flip = array_to_geoIm(Hist2d)


# A 2d mask
maskIm_grdInCity = grid_within_bndy([grd_x,grd_y],city_shp,im2d_mask=True,geoIm=True)['mask']
Hist2d_flip[maskIm_grdInCity==False] = np.nan

# plot
# ************************************************************ #
# numpy array shown as image
#                     geo_coord_x  109xx1  110xx5  110xx9 ...
#                     col_idx      0       1       2      ...
# geo_coord_y  row_idx                 
# 195xx1         0                 val_00  val_01  val_02 ...
# 194xx2         1                 val_10  val_11  val_12 ... 
# 193xx4         2                 val_20  val_21  val_22 ...
# ...           ...
# ************************************************************ #

import matplotlib.pyplot as plt

fig = plt.figure(figsize=(10, 3))
ax = fig.add_subplot(121)
im = plt.imshow(Hist2d, interpolation='nearest', origin='upper', cmap='jet',
                extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
ax.set_title('Binned burglary count (original array)')
ax = fig.add_subplot(122)
im = plt.imshow(Hist2d_flip, interpolation='nearest', origin='upper', cmap='jet',
                extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
ax.set_title('Binned burglary count (flipped array)')
plt.colorbar(im)