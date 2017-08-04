# -*- coding: utf-8 -*-
"""
================================
Explore spatial data
================================

Created on Wed Aug 17 11:52:25 2016

@author: xiaomuliu
"""

import fiona

shpfile_city = "/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/GISdata/City_Boundary/City_Boundary.shp"
chi_city = fiona.open(shpfile_city)
next(chi_city)
chi_city.close()
chi_city.closed

import pprint
from fiona.crs import to_string

with fiona.open(shpfile_city) as chi_city:
    pprint.pprint(chi_city[0])
    print(chi_city.driver)
    print(chi_city.crs)
    print(to_string(chi_city.crs))
    print(chi_city.bounds)

#GeoPandas extends the datatypes used by pandas to allow spatial operations on geometric types. 
#Geometric operations are performed by 'shapely'.
#Geopandas further depends on 'fiona' for file access and descartes and 'matplotlib' for plotting.

import geopandas as gpd
chi_city = gpd.read_file(shpfile_city)
chi_city.head()
chi_city.plot()

shpfile_beat = "/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/GISdata/cpd_beats/cpd_beats.shp"
beat = gpd.read_file(shpfile_beat)
beat.head()
beat.plot()

beat['centroid_column'] = beat.centroid
beat = beat.set_geometry('centroid_column')
beat.plot()

# GeoDataFrame keeps track of the active column by name, so if you rename the active geometry column,
# you must also reset the geometry:
beat = beat.rename(columns={'centroid_column': 'geometry'}).set_geometry('geometry')

# plot a subset of beat shapefile
district1 = beat[beat.DISTRICT=='01']
district1.plot(cmap='OrRd')

# overlay multiple layers
beat_ctr = beat.set_geometry('centroid_column')
beat_ctr.plot(marker='+', color='green', markersize=3)
beat_ctr = beat_ctr.to_crs(chi_city.crs) #ensure they share a common CRS 

# plotting method 1
base = chi_city.plot(color='white')
beat_ctr.plot(ax=base, marker='o', color='red', markersize=3)

# plotting method 2
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
# set aspect to equal. This is done automatically
# when using *geopandas* plot on it's own, but not when
# working with pyplot directly.
ax.set_aspect('equal')
chi_city.plot(ax=ax, color='white')
beat_ctr.plot(ax=ax, marker='o', color='red', markersize=3)

plt.show();