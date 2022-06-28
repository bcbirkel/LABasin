#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 12:11:35 2022

@author: bcbirkel
"""

# basic data plotting

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

plotTraces = True

#set corners for map
llcrnrlon=-119
llcrnrlat=33.25 #- 0.25
urcrnrlon=-117
urcrnrlat=34.75 #+ 0.25

events = ['lahabra_2014', 'beverlyhills_2001', 'chatsworth_2007', 'chinohills_2008', 'inglewood_2009']
ev_name = ['lahabra', 'beverlyhills', 'chatsworth', 'chinohills', 'inglewood']
models = ["CVM-S4","CVM-H"]

ks = [0,1,2,3,4]
ms = [0,1]

ks = [0]
ms = [0]

# %% Set parameters
#set event
for event_no in ks:
    for j in ms:
        event = events[event_no]
        name = ev_name[event_no]
        mod = models[j]

        # %% Set paths
        
        dataPath = './CompiledEvents/' + event + '/allData/rotated/'
        synPath = './CompiledEvents/' + event + '/GravesSyn/' + mod + '/' #'rotated/'
        
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
            
        st_data = read(dataPath + "*BH*.SAC")
        st_syn = read(synPath + "*.SAC")
        
        # # %% PLOT MAP
        
        # m = Basemap(projection='merc',llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,epsg=4269) #http://server.arcgisonline.com/arcgis/rest/services; EPSG Number of America is 4269
        # m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 2000, verbose= True)
    
        # for tr in st_data:
        #     if tr.stats.sac.stlo > 0:
        #         tr.stats.sac.stlo = -tr.stats.sac.stlo
        #     if tr.stats.sac.stla < 0:
        #         tr.stats.sac.stla = -tr.stats.sac.stla
        #     xx,yy = m(tr.stats.sac.stlo,tr.stats.sac.stla)
        #     m.scatter(xx,yy, marker = "o" ,s=2, color='k', alpha = 1)
        #     # print(tr.stats.sac.stlo,tr.stats.sac.stla)
        #     for sytr in st_syn:
        #         if tr.stats.station == sytr.stats.station:
        #             xx,yy = m(tr.stats.sac.stlo,tr.stats.sac.stla)
        #             m.scatter(xx,yy, marker = "o" ,s=0.1, color='b', alpha = 1)

        # xx,yy = m(event_lon,event_lat)
        # m.scatter(xx,yy, marker = "*" ,s=100, color='y', edgecolors = "k", alpha = 1)
    
        # plt.title(name + " " + mod + ': Data (black), Syn (blue)')
        # plt.show()
        # %% ###
        
        # read in station file
        stationFile = "./all_stationmaster.txt"
        stations = []
        savedpts = []
        stalist = []
        tmpsta = []
        
        # strip station file
        staCoord = open(stationFile, 'r')
        lines = staCoord.readlines()
        for line in lines:
            if not line.startswith("#"): # ignore first line
                split_line = line.split()
                net = split_line[0]
                code = split_line[1]
                lat = float(split_line[2])
                lon = float(split_line[3])
                stations.append([net,code,lat,lon])
                tmpsta.append(code)

                    
        # for sta in stations:           
        #     for tr in st_syn:
        #         if sta[1] == tr.stats.station:
        #             [dist,az,baz] = gps2dist_azimuth(event_lat,event_lon,sta[2],sta[3])
        #             distance = dist/1000
        #             tr.stats.distance = distance
        #             tr.stats.sac.dist = distance
        #             print(str(tr.stats.distance) + " " + tr.stats.station)
       
        # for tr in st_data:
        #     if tr.stats.sac.stlo > 0:
        #        tr.stats.sac.stlo = -tr.stats.sac.stlo 
        #     [dist,az,baz] = gps2dist_azimuth(event_lat,event_lon,tr.stats.sac.stla,tr.stats.sac.stlo)
        #     distance = dist/1000
        #     tr.stats.sac.dist = distance
        #     tr.stats.distance = distance
    
        # for station in tmpsta:
        #     for tr in st_data:
        #         if tr.stats.station == station:
        #             if str(tr)[0:2] == 'CI':
        #                 stalist.append(station)
        #             if str(tr)[0:2] == 'CE':
        #                 stalist.append(station)    
        #                 print(tr.stats.sac.dist)
                    
        stalist = set(tmpsta)

        
        # %% FIX
        
        # for tr in st_data: 
        #     if max(tr.data) > 100:
        #         print(str(tr) + " " + str(max(tr.data)))
    
        # %% PLOT TRACES
        gain = 1
        
        if plotTraces == True:
        
            figN = plt.figure(0, figsize=[15,5])
            aN = figN.add_subplot(111)
            
            figE = plt.figure(1, figsize=[15,5])
            aE = figE.add_subplot(111)
            
            # figR = plt.figure(2, figsize=[15,5])
            # aR = figR.add_subplot(111)
            
            # figT = plt.figure(3, figsize=[15,5])
            # aT = figT.add_subplot(111)
            
            figZ = plt.figure(4, figsize=[15,5])
            aZ = figZ.add_subplot(111)
            
            figs = [figN,figE,figZ]
            axes = [aN,aE,aZ]
            channels =["BHN","BHE","BHZ"]
            
            ymax = 5
            xmax = 50
            lowfreq = 1/5
            
            # for sytr in st_syn:
            #     sytr.stats.delta = 0.05
            #     sytr.stats.sampling_rate = 20
                
            for tr in st_data:
                if tr.stats.sac.dist > 500:
                    tr.stats.sac.dist = tr.stats.sac.dist/1000
            
            # st_data.filter('lowpass',freq=lowfreq)
            # st_syn.filter('lowpass',freq=lowfreq)
            
            for i in range(3):
                fig = figs[i]
                a = axes[i]
                cha = channels[i]
                
                a.set_ylabel("Distance from source (km)")    
                a.set_xlabel("Time (s)")
                a.set_xlim(-5,xmax)
                a.set_ylim(0,ymax)
                a.set_title(name + " " + mod + " " + cha + " - data")
                
                for tr in st_data:
                    if tr.stats.channel == cha and tr.stats.station in stalist and tr.stats.sac.dist < ymax:
                        
                        # if str(tr)[0:2] == 'CI': # and tr.stats.sac.dist < 50:
                        #     tr.filter("lowpass", freq=lowfreq)
                        #     a.plot(tr.times(reftime=start), tr.data*gain+tr.stats.sac.dist, c="r",linewidth=0.5)
                        #     a.text(xmax-5,tr.stats.sac.dist+0.2,tr.stats.station)
                        #     print(str(tr))
                        # if str(tr)[0:2] == 'CE':
                        #     tr.filter("lowpass", freq=lowfreq)
                            # a.plot(tr.times(reftime=start), tr.data*gain+tr.stats.sac.dist, c="y",linewidth=0.5)
                            # a.text(xmax-5,tr.stats.sac.dist+0.2,tr.stats.station)
                        a.plot(tr.times(reftime=start), tr.data*gain+tr.stats.sac.dist, c="k",linewidth=0.5)
                        a.text(xmax-5,tr.stats.sac.dist+0.2,tr.stats.station)
                        
                for sytr in st_syn:
                    if sytr.stats.channel == cha and sytr.stats.station in stalist and sytr.stats.sac.dist < ymax:
                        # if str(tr)[0:2] == 'CI':
                        #     tr.filter("lowpass", freq=lowfreq)
                        
                        a.plot(sytr.times(reftime=start), sytr.data*gain+sytr.stats.sac.dist, c="b",linewidth=0.5)
                            # a.text(xmax - 5,tr.stats.dist+0.2,tr.stats.station)
                            # print(str(tr))
                    
            plt.show()
    
    
    
        # %% PLOT TRACES
        
        # f = plt.figure(figsize=[30,15])
        
        # data_plot = Stream()
        # syn_plot = Stream()
        
        # for tr in st_data:
        #     if tr.stats.distance < 500:
        #         tr.stats.distance = tr.stats.distance*1000
        #     if tr.stats.channel.startswith('B') and tr.stats.channel.endswith('N'):
        #         data_plot.append(tr)
        # for tr in st_syn:
        #     if tr.stats.distance < 500:
        #         tr.stats.distance = tr.stats.distance*1000
        #     if tr.stats.channel.endswith('N'):
        #         tr.stats.starttime = start
        #         syn_plot.append(tr)
        
        # data_plot.filter('lowpass',freq=0.2)       
        # syn_plot.filter('lowpass',freq=0.2)          
        
        # data_plot.plot(type='section',fig=f,reftime=start,linewidth=1,scale=2,orientation='horizontal',recordlength=30,offset_max=10000)
        # syn_plot.plot(type='section',fig=f,color='r',reftime=start,linewidth=1,scale=2,orientation='horizontal',recordlength=30,offset_max=10000)
    
    
    
    
    
    
    
    
    
    
    
    
            