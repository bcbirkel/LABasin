#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 15:13:59 2021

@author: bcbirkel
"""


from obspy import read, read_inventory
from obspy.signal.filter import bandpass
# import pandas as pd
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
from mpl_toolkits.basemap import Basemap

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
    
for i in range(5):
    event_no = 0
    events = ['lahabra_2014', 'beverlyhills_2001', 'chatsworth_2007', 'chinohills_2008', 'inglewood_2009']
    event = events[event_no]
    
    
    os.chdir('/Users/bcbirkel/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/LABasin/CompiledEvents/' + event + '/allData/fixed')
    # %% Set paths
    
    path_dir_data = '/Users/bcbirkel/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/LABasin/CompiledEvents/' + event + '/allData/fixed/'
    # path_dir_syn = '/Users/bcbirkel/Documents/GitHub/LABasin/CompiledEvents/' + event + '/GravesSyn/CVM-S4/'

    if event_no == 0:
        event_title = 'lahabra'
        event_lat = 33.9325
        event_lon = -117.9158
    elif event_no == 3:
        event_title = 'chinohills'
        event_lat = 33.9465
        event_lon = -117.7667
    elif event_no == 4:
        event_title = 'inglewood'
        event_lat = 33.9377
        event_lon = -118.3357
    elif event_no == 2:
        event_title = 'chatsworth'
        event_lat = 34.2983
        event_lon = -118.6255
    elif event_no == 1:
        event_title = 'beverlyhills'
        event_lat = 34.0541
        event_lon = -118.3929

    else: 
        print('unknown event file')

    stream_data = read(path_dir_data + "*.SAC")
    # stream_syn = read(path_dir_syn + "*.SAC")

    for tr in stream_data:
        # if tr.stats.station[0].isdigit() == True:
        #       net.append("CE")
        net.append(tr.stats.network)
        code.append(tr.stats.station)
        lat.append(tr.stats.sac.stla)
        if tr.stats.sac.stlo > 0:
            tr.stats.sac.stlo = -tr.stats.sac.stlo
        lon.append(tr.stats.sac.stlo)
        
    # for tr in stream_syn:
    #     if tr.stats.station[0].isdigit() == True:
    #         net.append("CE")
    #         code.append(tr.stats.station)
    #         lat.append(tr.stats.sac.stla)
    #         if tr.stats.sac.stlo > 0:
    #             tr.stats.sac.stlo = -tr.stats.sac.stlo
    #         lon.append(tr.stats.sac.stlo)
  
#%% MAP    
  
# #set corners for map
# llcrnrlon=-119
# llcrnrlat=33.25
# urcrnrlon=-117
# urcrnrlat=34.75

# m = Basemap(projection='merc',llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,epsg=4269, resolution='l') #http://server.arcgisonline.com/arcgis/rest/services; EPSG Number of America is 4269
# m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 2000, verbose= True)

# m.scatter(lon, lat, c='b', marker = "o" , edgecolors = "k", alpha = 1, s=5)

#%%
for i in range(len(code)):        
    text.append(net[i] + " " + str(code[i]) + " " + str(lat[i]) + " " + str(lon[i]))
text = set(text)
text = list(text)
text.sort()

# write_text = []
# for i in range(len(text)):
#     write_text.append(str(text[i]) + "\n")
    
file = open("/Users/bcbirkel/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/LABasin/all_stationmaster.txt", "w")
for i in range(len(text)):
    print(str(text[i]) + "\n")
    file.write(str(text[i]) + "\n")
        
file.close()



        
        
        