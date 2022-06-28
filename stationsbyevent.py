#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 15:14:59 2021

@author: bcbirkel
"""

from obspy import read, read_inventory
from obspy.signal.filter import bandpass
import pandas as pd
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from obspy.signal.invsim import simulate_seismometer, corn_freq_2_paz
from obspy.geodetics import gps2dist_azimuth
from obspy.io.sac import sacpz
from obspy.io.sac import attach_paz
from obspy.io.gse2.paz import read_paz
import os
from obspy.core import Trace, Stream, AttribDict
from obspy.signal.rotate import rotate_ne_rt
import math
import obspy.realtime.signal as signal

#Ignore warnings due to python 2 and 3 conflict
import warnings
warnings.filterwarnings("ignore")

# %% Set parameters
#set event

net = []
code = []
lat = []
lon = []
text = []
    
for i in range(1):
    event_no = 4
    events = ['lahabra_2014', 'beverlyhills_2001', 'chatsworth_2007', 'chinohills_2008', 'inglewood_2009']
    event = events[event_no]
    
    
    os.chdir('/Users/bcbirkel/Documents/GitHub/LABasin/CompiledEvents/' + event + '/allData/fixed')
    # %% Set paths
    
    path_dir_data = '/Users/bcbirkel/Documents/GitHub/LABasin/CompiledEvents/' + event + '/allData/fixed/'
    # path_dir_syn = '/Users/bcbirkel/Documents/GitHub/LABasin/CompiledEvents/' + event + '/GravesSyn/CVM-S4/'

    if event_no == 0:
        event_title = 'lahabra'
        event_name = 'La Habra 2014'
        event_lat = 33.9325
        event_lon = -117.9158
    elif event_no == 3:
        event_title = 'chinohills'
        event_name = 'Chino Hills 2008'
        event_lat = 33.9465
        event_lon = -117.7667
    elif event_no == 4:
        event_title = 'inglewood'
        event_name ='Inglewood 2009'
        event_lat = 33.9377
        event_lon = -118.3357
    elif event_no == 2:
        event_title = 'chatsworth'
        event_name = 'Chatsworth 2007'
        event_lat = 34.2983
        event_lon = -118.6255
    elif event_no == 1:
        event_title = 'beverlyhills'
        event_name = 'Beverly Hills 2001'
        event_lat = 34.0541
        event_lon = -118.3929

    else: 
        print('unknown event file')

    stream_data = read(path_dir_data + "*.SAC")
    # stream_syn = read(path_dir_syn + "*.SAC")

    for tr in stream_data:
        net.append(tr.stats.network)
        code.append(tr.stats.station)
        lat.append(tr.stats.sac.stla)
        if tr.stats.sac.stlo > 0:
            tr.stats.sac.stlo = -tr.stats.sac.stlo
        lon.append(tr.stats.sac.stlo)


for i in range(len(code)):        
    text.append(net[i] + " " + str(code[i]) + " " + str(lat[i]) + " " + str(lon[i]))
text = set(text)
text = list(text)
text.sort()
    
file = open("/Users/bcbirkel/Documents/GitHub/LABasin/" + event_title + "_stations.txt", "w")
for i in range(len(text)):
    print(str(text[i]) + "\n")
    file.write(str(text[i]) + "\n")
        
file.close()


# %% MAP

#set corners for map
llcrnrlon=-119
llcrnrlat=33.25 #- 0.25
urcrnrlon=-117
urcrnrlat=34.75 #+ 0.25

fig = plt.figure()
ax = fig.add_subplot(111)
# setup mercator map projection.
m = Basemap(projection='merc',llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,epsg=4269) #http://server.arcgisonline.com/arcgis/rest/services; EPSG Number of America is 4269
m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 2000, verbose= True)

# plot stations
for i in range(len(code)):
    xx,yy = m(lon[i],lat[i])
    if net[i] == "CE":
        m.scatter(xx, yy, marker = "^" ,s=30, c='g', edgecolors = "k", alpha = 1)
    elif net[i] == "CI":
        m.scatter(xx, yy, marker = "^" ,s=30, c='b', edgecolors = "k", alpha = 1)
    else:
        m.scatter(xx, yy, marker = "^" ,s=30, c='y', edgecolors = "k", alpha = 1)

#Plot the event
xx,yy = m(event_lon,event_lat)
m.scatter(xx, yy, marker = "*" ,s=200, c="m" , edgecolors = "k", alpha = 1) 

plt.title(event_name + " Data Stations")
plt.xlabel("green=CE, blue=CI, yellow=other",fontsize=10)  
plt.show()
        
        
        