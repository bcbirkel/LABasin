#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 09:24:25 2021

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
cha = []
lat = []
lon = []
elev = []
text = []
    
for i in range(1):
    event_no = 0
    events = ['lahabra_2014', 'beverlyhills_2001', 'chatsworth_2007', 'chinohills_2008', 'inglewood_2009']
    event = events[event_no]
    
    
    os.chdir('/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/currentSim/ObsData/Velocity')
    # %% Set paths
    
    path_dir_data = '/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/currentSim/ObsData/Velocity/'
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

    for tr in stream_data:
        net.append(tr.stats.network)
        code.append(tr.stats.station)
        cha.append(tr.stats.channel)
        lat.append(tr.stats.sac.stla)
        if tr.stats.sac.stlo > 0:
            tr.stats.sac.stlo = -tr.stats.sac.stlo
        lon.append(tr.stats.sac.stlo)
        if hasattr(tr.stats.sac, 'stel') == True:
            elev.append(tr.stats.sac.stel)
        else:
            elev.append(0)
        
for i in range(len(code)):        
    text.append(net[i] + " " + str(code[i]) + " " + str(cha[i]) + " "+ str(lat[i]) + " " + str(lon[i]) + " " + str(elev[i]))
text = set(text)
text = list(text)
text.sort()

file = open("/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/currentSim/Resources/Catalogs/StationList.txt", "w")
for i in range(len(text)):
    print(str(text[i]) + "\n")
    file.write(str(text[i]) + "\n")
        
file.close()



        
        
        