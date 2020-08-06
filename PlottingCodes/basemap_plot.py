#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 10:35:26 2020

@author: bcbirkel
"""

from obspy import read
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

#Ignore warnings due to python 2 and 3 conflict
import warnings
warnings.filterwarnings("ignore")


# make basemap figure
plt.figure(figsize=(10,5))

# setup mercator map projection.
m = Basemap(projection='merc',llcrnrlon=-120, llcrnrlat=33, urcrnrlon=-117, urcrnrlat=35,epsg=4269)

# add imagery layer
#http://server.arcgisonline.com/arcgis/rest/services -- EPSG Number for America is 4269 (see above)
m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 2000, verbose= True)

#load data in SAC format (with metadata in header, comes from IRIS like this)
path_dir = '/Users/bcbirkel/Documents/Research/LABasin/rawIRISdata/20200422/'
stream = read(path_dir + "C*")


# plot stations
for tr in stream:
    stlat = tr.stats.sac.stla; stlon = tr.stats.sac.stlo 
    #m.drawgreatcircle(stlon,stlat,evlon,evlat,linewidth=1,color='b')
    xx,yy = m(stlon,stlat)
    m.scatter(xx, yy, marker = "^" ,s=100, c="g" , edgecolors = "k", alpha = 1)
 
# plot the event
evlat = stream[0].stats.sac.evla; evlon = stream[0].stats.sac.evlo
xx,yy = m(evlon,evlat)
m.scatter(xx, yy, marker = "*" ,s=200, c="r" , edgecolors = "k", alpha = 1)    
    
plt.title("Event-station map")
plt.show()