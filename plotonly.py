#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 13:50:59 2021

@author: bcbirkel
"""
from obspy import read, read_inventory
import obspy.signal.filter as flt
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
from obspy.core.utcdatetime import UTCDateTime
import shutil
import matplotlib.cm as cm
import matplotlib

#Ignore warnings due to python 2 and 3 conflict
import warnings
warnings.filterwarnings("ignore")

# TOGGLES

plotMap = False
plotTraces = True
plotBeams = False

#set corners for map
llcrnrlon=-119
llcrnrlat=33.25 #- 0.25
urcrnrlon=-117
urcrnrlat=34.75 #+ 0.25

rotate = False

# filetitle = "Beam_Grid_lahabra_CVM-S4_"
events = ['lahabra_2014', 'beverlyhills_2001', 'chatsworth_2007', 'chinohills_2008', 'inglewood_2009']
ev_name = ['lahabra', 'beverlyhills', 'chatsworth', 'chinohills', 'inglewood']
models = ["CVM-S4","CVM-H"]

ks = [0,1,2,3,4]
ms = [0,1]

ks = [1]
ms = [0]

# %% Set parameters
#set event
for event_no in ks:
    for j in ms:
        event = events[event_no]
        name = ev_name[event_no]
        mod = models[j]

        # %% Set paths
        
        dataPath = '/Users/bcbirkel/Documents/GitHub/LABasin/CompiledEvents/' + event + '/allData/fixed/'
        synPath = '/Users/bcbirkel/Documents/GitHub/LABasin/CompiledEvents/' + event + '/GravesSyn/' + mod + '/rotated/'
        beamPath = '/Users/bcbirkel/Documents/GitHub/LABasin/CompiledEvents/' + event + '/GravesSyn/' + mod + '/beamforming/'
        
        if event_no == 0:
            event_title = 'lahabra'
            event_name = 'La Habra'
            event_lat = 33.9325
            event_lon = -117.9158
            start = UTCDateTime("2014-03-29T04:09:42.994500Z")   
            
        elif event_no == 1:
            event_title = 'beverlyhills'
            event_name = 'Beverly Hill'
            event_lat = 34.0541
            event_lon = -118.3929
            start = UTCDateTime("2001-09-09T23:59:17.695Z")
            
        elif event_no == 2:
            event_title = 'chatsworth'
            event_name = 'Chatsworth 2007'
            event_lat = 34.2983
            event_lon = -118.6255
            start = UTCDateTime("2007-08-09T07:58:48.888")
            
        elif event_no == 3:
            event_title = 'chinohills'
            event_name = 'Chino Hills'
            event_lat = 33.9465
            event_lon = -117.7667
            start = UTCDateTime("2008-07-29T18:42:15.960")
    
        elif event_no == 4:
            event_title = 'inglewood'
            event_name = 'Inglewood 2009'
            event_lat = 33.9377
            event_lon = -118.3357
            start = UTCDateTime("2009-05-18T03:39:36.126")
            
        else: 
            print('unknown event file')
            
        st_data = read(dataPath + "*.SAC")
        st_syn = read(synPath + "*.SAC")
        st_beam = read(beamPath + "*.SAC")
        
        # %% PLOT MAP
        
        if plotMap == True:
            m = Basemap(projection='merc',llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,epsg=4269) #http://server.arcgisonline.com/arcgis/rest/services; EPSG Number of America is 4269
            m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 2000, verbose= True)
        
            for tr in st_data:
                if tr.stats.sac.stlo > 0:
                    tr.stats.sac.stlo = -tr.stats.sac.stlo
                if tr.stats.sac.stla < 0:
                    tr.stats.sac.stla = -tr.stats.sac.stla
                xx,yy = m(tr.stats.sac.stlo,tr.stats.sac.stla)
                m.scatter(xx,yy, marker = "o" ,s=2, color='r', alpha = 1)
                # print(tr.stats.sac.stlo,tr.stats.sac.stla)
                for sytr in st_syn:
                    if tr.stats.station == sytr.stats.station:
                        xx,yy = m(tr.stats.sac.stlo,tr.stats.sac.stla)
                        m.scatter(xx,yy, marker = "o" ,s=0.1, color='b', alpha = 1)
                
            # for tr in st_syn:
            #     xx,yy = m(tr.stats.sac.stlo,tr.stats.sac.stla)
            #     m.scatter(xx,yy, marker = "o" ,s=0.1, color='b', alpha = 1)
                
            for tr in st_beam:
                xx,yy = m(tr.stats.sac.stlo,tr.stats.sac.stla)
                m.scatter(xx,yy, marker = "o" ,s=0.01, color='k', alpha = 1)
                    
            xx,yy = m(event_lon,event_lat)
            m.scatter(xx,yy, marker = "*" ,s=100, color='y', edgecolors = "k", alpha = 1)
            # plt.text(event_lon+0.03, event_lat-0.07, name)
        
            plt.title(name + " " + mod + ': Data (red), Syn (blue), Beam (black)')
            plt.show()
        
        # %% PLOT TRACES
        
        if plotTraces == True:
        
            figN = plt.figure(0, figsize=[15,5])
            aN = figN.add_subplot(131)
            sN = figN.add_subplot(132)
            bN = figN.add_subplot(133)
            
            figE = plt.figure(1, figsize=[15,5])
            aE = figE.add_subplot(131)
            sE = figE.add_subplot(132)
            bE = figE.add_subplot(133)
            
            figR = plt.figure(2, figsize=[15,5])
            aR = figR.add_subplot(131)
            sR = figR.add_subplot(132)
            bR = figR.add_subplot(133)
            
            figT = plt.figure(3, figsize=[15,5])
            aT = figT.add_subplot(131)
            sT = figT.add_subplot(132)
            bT = figT.add_subplot(133)
            
            figZ = plt.figure(4, figsize=[15,5])
            aZ = figZ.add_subplot(131)
            sZ = figZ.add_subplot(132)
            bZ = figZ.add_subplot(133)
            
            figs = [figN,figE,figR,figT,figZ]
            axes = [aN,aE,aR,aT,aZ]
            syax = [sN,sE,sR,sT,sZ]
            bax = [bN,bE,bR,bT,bZ]
            channels = ["BHN","BHE","BHR","BHT","BHZ"]
            
            for i in range(5):
                fig = figs[i]
                a = axes[i]
                s = syax[i]
                b = bax[i]
                cha = channels[i]
                
                a.set_ylabel("Distance from source (km)")    
                a.set_xlabel("Time (s)")
                a.set_xlim(-10,60)
                # a.set_ylim(0,50)
                a.set_title(name + " " + mod + " " + cha + " - data")
                    
                s.set_xlabel("Time (s)")
                s.set_xlim(-10,60)
                # s.set_ylim(0,50)
                s.set_title(name + " " + mod + " " +  cha + " - synthetics")
                
                b.set_xlabel("Time (s)")
                b.set_xlim(-10,60)
                # b.set_ylim(0,50)
                b.set_title(name + " " + mod + " " + cha + " - beamforming")
                
                # stalist = ['LGB','08x19']
                
                for tr in st_data:
                    if tr.stats.channel == cha: # and tr.stats.station in stalist:
                        if tr.stats.network == 'CI':    
                            a.plot(tr.times()-10, tr.data+tr.stats.sac.dist, c="k",linewidth=0.5)
                        else:
                            a.plot(tr.times(), tr.data+tr.stats.sac.dist, c="b",linewidth=0.5)
                        # a.plot(tr.times(), tr.data, c="k",linewidth=0.5)
                        # a.plot(tr.data+tr.stats.sac.dist, c="k",linewidth=0.5)
                        a.text(50,tr.stats.sac.dist+0.2,tr.stats.station)
                        
                for tr in st_syn:
                    if tr.stats.channel == cha: # and tr.stats.station in stalist:
                        s.plot(tr.times(), tr.data+tr.stats.sac.dist, c="b",linewidth=0.5)
                        # s.plot(tr.data+tr.stats.sac.dist, c="b",linewidth=0.5)
                        s.text(50,tr.stats.sac.dist+0.2,tr.stats.station)
                        
                for tr in st_beam:
                    if tr.stats.channel == cha: # and tr.stats.station in stalist:
                        b.plot(tr.times(), tr.data+tr.stats.sac.dist, c="b",linewidth=0.5)
                        # b.plot(tr.data+tr.stats.sac.dist, c="b",linewidth=0.5)
                        b.text(50,tr.stats.sac.dist+0.2,tr.stats.station)
                        
            plt.show()
        
        # %% PLOT BEAMFORMING ONLY
        
        if plotBeams == True: 
            figBeam = plt.figure(5, figsize=[25,15])
            bN = figBeam.add_subplot(231)
            bE = figBeam.add_subplot(232)
            bR = figBeam.add_subplot(233)
            bT = figBeam.add_subplot(234)
            bZ = figBeam.add_subplot(235)
            
            bax = [bN,bE,bR,bT,bZ]
            channels = ["BHN","BHE","BHR","BHT","BHZ"]
            
            for i in range(5):
                b = bax[i]
                cha = channels[i]
                
                b.set_ylabel("Distance from source (km)") 
                b.set_xlabel("Time (s)")
                b.set_xlim(-10,60)
                # b.set_ylim(0,50)
                b.set_title(name + " " + mod + " " + cha + " - beamforming")
                
                for tr in st_beam:
                    if tr.stats.channel == cha:
                        # tr.data = np.pad(tr.data,(int(20*10),0),'constant',constant_values=(0))
                        b.plot(tr.times(), tr.data+tr.stats.sac.dist, c="k",linewidth=0.5)
                        
            plt.show()
                        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
            