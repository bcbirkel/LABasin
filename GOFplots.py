#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 12:11:48 2022

@author: bcbirkel
"""

from obspy import read
from obspy.core import Trace, Stream
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import os
from obspy.geodetics import gps2dist_azimuth
from numpy import meshgrid
import matplotlib.mlab as mlab
import scipy
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from collections import Counter

# %% Set variables

figPath = '../../Research/currentFigs/contourPlotFigs/GOF/'
stationsZ = []
stationsR = []
stationsT = []
rstationsZ = []
rstationsR = []
rstationsT = []


for ii in range(1):  
    event_no = ii
    events = ['lahabra_2014', 'chatsworth_2007', 'chinohills_2008', 'inglewood_2009', 'beverlyhills_2001']
    event = events[event_no]
    
    dataPath = './CompiledEvents/' + event + '/allData/fixed/'
    rgPath = './CompiledEvents/' + event + '/GravesSyn/synthetics_stations/CVM-S4/'
    
    low_freq = 1/20
    high_freq = 1/5
    
    if event_no == 0:
        event_title = 'lahabra'
        event_lat = 33.9325
        event_lon = -117.9158
        gain = 2/high_freq
    elif event_no == 2:
        event_title = 'chinohills'
        event_lat = 33.9465
        event_lon = -117.7667
        gain = 100*high_freq
    elif event_no == 3:
        event_title = 'inglewood'
        event_lat = 33.9377
        event_lon = -118.3357
        gain = 500*high_freq
    elif event_no == 1:
        event_title = 'chatsworth'
        event_lat = 34.2983
        event_lon = -118.6255
        gain = 500*high_freq
    # elif event_no == 1:
    #     event_title = 'beverlyhills'
    #     event_lat = 34.0541
    #     event_lon = -118.3929
    else: 
        print('unknown event file')
    
    
    #set corners for map
    llcrnrlon=-119
    llcrnrlat=33.25
    urcrnrlon=-117
    urcrnrlat=34.75
    

    

    # %% Import data
    
    stream_dataZ = read(dataPath + "*Z_vel.SAC") 
    for tr in stream_dataZ:
        tr.stats.channel = 'Z'
    stream_dataR = read(dataPath + "*R_vel.SAC") 
    for tr in stream_dataR:
        tr.stats.channel = 'R'
    stream_dataT = read(dataPath + "*T_vel.SAC") 
    for tr in stream_dataT:
        tr.stats.channel = 'T'   ,
        
    
    stream_rgZ = read(rgPath + "*Z_vel.SAC") 
    for tr in stream_rgZ:
        tr.stats.channel = 'Z'
    stream_rgR = read(rgPath + "*R_vel.SAC") 
    for tr in stream_rgR:
        tr.stats.channel = 'R'
    stream_rgT = read(rgPath + "*T_vel.SAC") 
    for tr in stream_rgT:
        tr.stats.channel = 'T'
        
    # %% Stations
    stationFile = "./all_stationmaster.txt"
    # stationFile = "./stationLists/stationsByEvent/lahabra_stations.txt"
    
    # Load station coords into arrays, many more stations than used
    st_netw  = []
    st_name  = []
    st_dist  = []
    st_az    = []
    st_baz   = []
    st_lat   = []
    st_lon   = []
    
    
    # strip station file
    staCoord = open(stationFile, 'r')
    lines = staCoord.readlines()
    for line in lines:
        split_line = line.split()
        st_netw.append(split_line[0])
        st_name.append(split_line[1])
        st_lat.append(float(split_line[2]))
        st_lon.append(float(split_line[3]))
        [distance,az,baz] = gps2dist_azimuth(event_lat, event_lon, float(split_line[2]), float(split_line[3])) # Get traveltime and azimuth
        print(distance/1000.,az,baz)
        st_dist.append(distance/1000.) # distance
        st_az.append(az) # azimuth
        st_baz.append(baz) # back-azimuth
 
    # %% set up lists for averages, streams + bandpass
        
    stream_data = stream_dataR + stream_dataT + stream_dataZ
    stream_rg = stream_rgZ + stream_rgR + stream_rgT
    
    stream_data.filter("bandpass", freqmin=low_freq, freqmax=high_freq)
    stream_rg.filter("bandpass", freqmin=low_freq, freqmax=high_freq)
    
    channels = ["Z","R","T"]
    
    
    for j in range(len(channels)):
        for tr in stream_data:
            if tr.stats.channel == channels[j]:
                if j == 0:
                    stationsZ.append(tr.stats.station)
                if j == 1:
                    stationsR.append(tr.stats.station)
                if j == 2:
                    stationsT.append(tr.stats.station)
        for tr in stream_rg:
            if tr.stats.channel == channels[j]:
                if j == 0:
                    rstationsZ.append(tr.stats.station)
                if j == 1:
                    rstationsR.append(tr.stats.station)
                if j == 2:
                    rstationsT.append(tr.stats.station)            
        
    if ii == 3:
        z_stacount = Counter(stationsZ)
        r_stacount = Counter(stationsR)
        t_stacount = Counter(stationsT)
        
        rz_stacount = Counter(rstationsZ)
        rr_stacount = Counter(rstationsR)
        rt_stacount = Counter(rstationsT)
            
    dataall = [[0 for j in range(len(st_name))] for i in range(len(channels))]
    synall = [[0 for j in range(len(st_name))] for i in range(len(channels))]


    da_st = [] 
    rg_st = [] 
    for tr in stream_data:
        da_st.append(tr.stats.station)
    for tr in stream_rg:
        rg_st.append(tr.stats.station)    
        
    st_ov = set(da_st).intersection(set(rg_st))
    
    st_overlap = list(st_ov)
    
    stov_lat =  [0 for j in range(len(st_overlap))]
    stov_lon =  [0 for j in range(len(st_overlap))]
    
    for i in range(len(st_name)):
        for j in range(len(st_overlap)):
            if st_name[i] == st_overlap[j]:
                stov_lat[j] = st_lat[i]
                stov_lon[j] = st_lon[i]

    # %% GOF contour plots
    
    for i in range(len(channels)):
        cha = channels[i]

        fig = plt.figure()
        ax = fig.add_subplot(111)
        xx = []
        yy = []
        datamax = []
        points = ([])
        all_d = []

        for dtr in stream_data:
            all_d.append(max(abs(dtr.data)))
                 
        # set bounds for outliers 
        upper = (np.mean(all_d)+2*np.std(all_d)) 
        lower = (np.mean(all_d)-2*np.std(all_d))
    
        synmax = []
        points = ([])
        
        vel_gof = []
                
        for i in range(len(st_overlap)):
            for dtr in stream_data:
                if dtr.stats.station == st_overlap[i] and dtr.stats.channel == cha:
                    if (max(abs(dtr.data))) < upper and (max(abs(dtr.data))) > lower:
                        for rgtr in stream_rg:
                            if rgtr.stats.station == st_overlap[i] and rgtr.stats.channel == cha:   
                                if (max(abs(rgtr.data))) < upper and (max(abs(rgtr.data))) > lower:
                                    yy.append(stov_lat[i])
                                    xx.append(stov_lon[i])
                                    datamax.append((max(abs(dtr.data))))
                                    synmax.append((max(abs(rgtr.data))))

        for i in range(len(datamax)):
            vel_gof.append(10*np.exp(-((datamax[i]-synmax[i])/(np.minimum(datamax[i],synmax[i])))**2))
                        
        # setup mercator map projection.
        m = Basemap(projection='merc',llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,epsg=4269, resolution='l') #http://server.arcgisonline.com/arcgis/rest/services; EPSG Number of America is 4269
        m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 2000, verbose= True)
        # m.drawcoastlines()
        
        # linear interpolation between stations at grid points
        x = (xx)
        y = (yy)
        yi = np.linspace(llcrnrlat, urcrnrlat, 100)
        xi = np.linspace(llcrnrlon, urcrnrlon, 100)
        xi, yi = np.mgrid[-120:-115:100j, 32:35:100j]
        
        zid = scipy.interpolate.griddata((x,y), vel_gof, (xi,yi), method='linear')
        
        # create contour plot
        plt.contourf(xi, yi, zid, alpha=.5, antialiased=True, cmap='jet') #, vmin=0, vmax=5)
        plt.colorbar(orientation="horizontal",shrink=2/3)  # draw colorbar  
        
        
        m.scatter(x, y, c=vel_gof, cmap='jet', marker = "^" , edgecolors = "k", alpha = 1, s=40)
    
        # plot the event
        xx,yy = m(event_lon,event_lat)
        m.scatter(xx, yy, marker = "*" ,s=50, c="m" , edgecolors = "k", alpha = 1) 
        # plt.title(event_name + " Transverse Data/Synthetic Vel Ratio - " + str(1/high_freq) + "-" + str(1/low_freq) + " sec")
        # plt.title("La Habra - Data Peak Velocity for " + cha + "-comp")
        plt.title(event_title + " - Peak Velocity GOF " + cha + " - " + str(1/high_freq) + "-" + str(1/low_freq) + " sec")
        # plt.show()
        # plt.savefig(figPath + event_title + " " + cha + " - peak velocity GOF.png")
        
        
        # %% plot traces
        ymax = 50 
        xmax = 120
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        for tr in stream_data:
            if tr.stats.channel == cha and tr.stats.sac.dist < ymax and tr.stats.station in st_overlap:
                # if (max(abs(dtr.data))) < upper and (max(abs(dtr.data))) > lower:
                    ax.plot(tr.times(), tr.data*gain+tr.stats.sac.dist,color='k',lw=0.3)
                    ax.text(xmax+20,tr.stats.sac.dist,tr.stats.station)
        for tr in stream_rg:
            if tr.stats.channel == cha and tr.stats.sac.dist < ymax and tr.stats.station in st_overlap:
                # if (max(abs(dtr.data))) < upper and (max(abs(dtr.data))) > lower:
                    ax.plot(tr.times(), tr.data*gain+tr.stats.sac.dist,color='b',lw=0.3)
        plt.ylim(0,ymax)
        plt.xlim(0,xmax)
        plt.title(event_title + " - Data v Synthetic " + cha)
       
