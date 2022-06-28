#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 11:36:33 2022

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

figPath = '../../Research/currentFigs/contourPlotFigs/'
stationsZ = []
stationsR = []
stationsT = []
rstationsZ = []
rstationsR = []
rstationsT = []


for ii in range(4):  
    event_no = ii
    events = ['lahabra_2014', 'chatsworth_2007', 'chinohills_2008', 'inglewood_2009', 'beverlyhills_2001']
    event = events[event_no]
    
    dataPath = './CompiledEvents/' + event + '/allData/fixed/'
    rgPath = './CompiledEvents/' + event + '/GravesSyn/synthetics_stations/CVM-S4/'
    
    if event_no == 0:
        event_title = 'lahabra'
        event_lat = 33.9325
        event_lon = -117.9158
    elif event_no == 2:
        event_title = 'chinohills'
        event_lat = 33.9465
        event_lon = -117.7667
    elif event_no == 3:
        event_title = 'inglewood'
        event_lat = 33.9377
        event_lon = -118.3357
    elif event_no == 1:
        event_title = 'chatsworth'
        event_lat = 34.2983
        event_lon = -118.6255
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
    
    low_freq = 1/5
    high_freq = 1/1
    

    # %% Import data
    
    stream_dataZ = read(dataPath + "*Z_vel.SAC") 
    for tr in stream_dataZ:
        tr.stats.channel = 'Z'
    stream_dataR = read(dataPath + "*R_vel.SAC") 
    for tr in stream_dataR:
        tr.stats.channel = 'R'
    stream_dataT = read(dataPath + "*T_vel.SAC") 
    for tr in stream_dataT:
        tr.stats.channel = 'T'   
        
    
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


    # %% data contour plots
    
    for i in range(len(channels)):
        cha = channels[i]
        
        # plot traces by distance
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for tr in stream_data:
            if tr.stats.channel == cha and tr.stats.sac.dist < 100:
                ax.plot(tr.times(), tr.data+tr.stats.sac.dist,color='k',lw=0.1)
                ax.text(220,tr.stats.sac.dist,tr.stats.station)
        plt.ylim(0,100)
        plt.xlim(0,200)
        plt.title(event_title + " - Data " + cha)
        # plt.show()
        plt.savefig(figPath + event_title + " " + cha + " - all data.png")
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        xx = []
        yy = []
        datamax = []
        points = ([])
        all_d = []
        doutliers = []
        doxx = []
        doyy = []
        
        for dtr in stream_data:
            all_d.append(max(abs(dtr.data)))
          
            
        # set bounds for outliers 
        upper = (np.mean(all_d)+5*np.std(all_d)) 
        lower = (np.mean(all_d)-5*np.std(all_d))
        
        # setup mercator map projection.
        m = Basemap(projection='merc',llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,epsg=4269, resolution='l') #http://server.arcgisonline.com/arcgis/rest/services; EPSG Number of America is 4269
        m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 2000, verbose= True)
        # m.drawcoastlines()
    
        # attach data with stations in seperate lists
        for i in range(len(st_name)):
            for dtr in stream_data:
                if dtr.stats.station == st_name[i] and dtr.stats.channel == cha:
                    if (max(abs(dtr.data))) < upper and (max(abs(dtr.data))) > lower:
                        yy.append(st_lat[i])
                        xx.append(st_lon[i])
                        datamax.append((max(abs(dtr.data))))
                    if (max(abs(dtr.data))) > upper or (max(abs(dtr.data))) < lower:
                        doyy.append(st_lat[i])
                        doxx.append(st_lon[i])
                        doutliers.append(dtr)
        
        # linear interpolation between stations at grid points
        x = (xx)
        y = (yy)
        yi = np.linspace(llcrnrlat, urcrnrlat, 100)
        xi = np.linspace(llcrnrlon, urcrnrlon, 100)
        xi, yi = np.mgrid[-120:-115:100j, 32:35:100j]
        
        zid = scipy.interpolate.griddata((x,y), datamax, (xi,yi), method='linear')
        
        # create contour plot
        plt.contourf(xi, yi, zid) #, vmin=0, vmax=5)
        plt.colorbar(orientation="horizontal",shrink=2/3)  # draw colorbar  
    
        # plot the event
        xx,yy = m(event_lon,event_lat)
        m.scatter(xx, yy, marker = "*" ,s=50, c="m" , edgecolors = "k", alpha = 1) 
        # plt.title(event_name + " Transverse Data/Synthetic Vel Ratio - " + str(1/high_freq) + "-" + str(1/low_freq) + " sec")
        # plt.title("La Habra - Data Peak Velocity for " + cha + "-comp")
        plt.title(event_title + " - Data Peak Velocity for " + cha + " - " + str(1/high_freq) + "-" + str(1/low_freq) + " sec")
        # plt.show()
        plt.savefig(figPath + event_title + " " + cha + " - data contour.png")
        
        # add to list for average plot, created later
        dataspreadall = []
        dataresidualall = []
        for i in range(len(st_name)):
            for dtr in stream_data:
                if dtr.stats.station == st_name[i] and dtr.stats.channel == cha:
                    if (max(abs(dtr.data))) < upper and (max(abs(dtr.data))) > lower:
                        # effectively take out variation due to event magnitude:
                        datanorm = (max(abs(dtr.data)))/(max(datamax)) #normalize by event maxima
                        # take out geometric spreading:
                        # dataspread = (datanorm*dtr.stats.sac.dist)
                        dataspreadall.append(datanorm*dtr.stats.sac.dist)
                        
                    
        for i in range(len(st_name)):
            for dtr in stream_data:
                if dtr.stats.station == st_name[i] and dtr.stats.channel == cha:
                    if (max(abs(dtr.data))) < upper and (max(abs(dtr.data))) > lower:
                        
                        
                        # effectively take out variation due to event magnitude:
                        datanorm = (max(abs(dtr.data)))/(max(datamax)) #normalize by event maxima
                        # take out geometric spreading:
                        dataspread = datanorm*dtr.stats.sac.dist
                        # take out baseline average:
                        dataresidual = dataspread-(sum(dataspreadall)/len(dataspreadall))

                        dataall[channels.index(cha)][i] += dataresidual
                        
            if ii == 3 and z_stacount[st_name[i]] != 0 and r_stacount[st_name[i]] != 0 and t_stacount[st_name[i]] != 0:
                dataall[0][i] = (dataall[0][i])/(z_stacount[st_name[i]])
                dataall[1][i] = (dataall[1][i])/(r_stacount[st_name[i]])
                dataall[2][i] = (dataall[2][i])/(t_stacount[st_name[i]])
        
        if len(doutliers) > 0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            for tr in doutliers:
                ax.plot(tr.times(), tr.data+tr.stats.sac.dist,color='k',lw=1)
                ax.text(max(tr.times())-2,tr.stats.sac.dist,tr.stats.station + ", " + str(max(tr.data)))
            plt.title(event_title + " - Data Outliers " + cha + " - Mean: " + str(np.mean(all_d)) + ", StDev: " + str(np.std(all_d)))
            # plt.show()
                
            mo = Basemap(projection='merc',llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,epsg=4269, resolution='l') #http://server.arcgisonline.com/arcgis/rest/services; EPSG Number of America is 4269
            mo.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 2000, verbose= True)
            mo.scatter(doxx, doyy, marker = "^" ,s=50, c="r" , edgecolors = "k", alpha = 1) 
            plt.title(event_title + " - Data Outlier Stations for " + cha)
            # plt.show()
            
# %%% synthetic plots
    
    # for i in range(len(channels)):
        
        # cha = channels[i]
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for tr in stream_rg:
            if tr.stats.channel == cha and tr.stats.sac.dist < 100:
                ax.plot(tr.times(), tr.data+tr.stats.sac.dist,color='k',lw=0.1)
                ax.text(220,tr.stats.sac.dist,tr.stats.station)
        plt.ylim(0,100)
        plt.xlim(0,200)
        plt.title(event_title + " - Synthetic " + cha)
        # plt.show()
        plt.savefig(figPath + event_title + " " + cha + " - all synthetic.png")
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        xx = []
        yy = []
        synmax = []
        points = ([])
        all_rg = []
        soutliers = []
        soxx = []
        soyy = []
        
        for rgtr in stream_rg:
            all_rg.append(max(abs(rgtr.data)))
            
        upper = (np.mean(all_rg)+10*np.std(all_rg)) 
        lower = (np.mean(all_rg)-10*np.std(all_rg))
        
        
        
        # setup mercator map projection.
        m2 = Basemap(projection='merc',llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,epsg=4269, resolution='l') #http://server.arcgisonline.com/arcgis/rest/services; EPSG Number of America is 4269
        m2.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 2000, verbose= True)
        # m.drawcoastlines()
        # m.etopo(service='ESRI_Imagery_World_2D', xpixels = 2000, verbose= True)
        
        for i in range(len(st_name)):
            for rgtr in stream_rg:
                if rgtr.stats.station == st_name[i] and rgtr.stats.channel == cha: # and tr.stats.network == 'CI':
                            # for datr in stream_data:
                            #     if datr.stats.channel == cha and datatr.stats.station == st_name[i]:
                    if (max(abs(rgtr.data))) < upper and (max(abs(rgtr.data))) > lower:
                        yy.append(st_lat[i])
                        xx.append(st_lon[i])
                        synmax.append((max(abs(rgtr.data))))
                    if (max(abs(rgtr.data))) > upper or (max(abs(rgtr.data))) < lower:
                        soyy.append(st_lat[i])
                        soxx.append(st_lon[i])
                        soutliers.append(rgtr)
        
        x = (xx)
        y = (yy)
        yi = np.linspace(llcrnrlat, urcrnrlat, 100)
        xi = np.linspace(llcrnrlon, urcrnrlon, 100)
        xi, yi = np.mgrid[-120:-115:100j, 32:35:100j]
        
        zi = scipy.interpolate.griddata((x,y), synmax, (xi,yi), method='linear')
        # levels = MaxNLocator(nbins=15).tick_values(zi.min(), zi.max())
        # levels = [0:0.1:3]
        
        plt.contour(xi, yi, zi) #, vmin=0, vmax=2)
        # plt.contourf(xi, yi, zi) #, norm=plt.Normalize(vmax=1, vmin=0))
        plt.colorbar(orientation="horizontal",shrink=2/3)  # draw colorbar  
    
        # #Plot the event
        xx,yy = m(event_lon,event_lat)
        m2.scatter(xx, yy, marker = "*" ,s=50, c="m" , edgecolors = "k", alpha = 1) 
        plt.title(event_title + " - Synthetic Peak Velocity for " + cha + " - " + str(1/high_freq) + "-" + str(1/low_freq) + " sec")
        # plt.show()
        plt.savefig(figPath + event_title + " " + cha + " - synthetic contour.png")
        
        # add to list for average plot, created later
        synspreadall = []
        for i in range(len(st_name)):
            for rtr in stream_data:
                if rtr.stats.station == st_name[i] and rtr.stats.channel == cha:
                    if (max(abs(rtr.data))) < upper and (max(abs(rtr.data))) > lower:
                        # effectively take out variation due to event magnitude:
                        synnorm = (max(abs(rtr.data)))/(max(synmax)) #normalize by event maxima
                        # take out geometric spreading:
                        synspreadall.append(synnorm*rtr.stats.sac.dist)
          
        for i in range(len(st_name)):
            for rtr in stream_rg:
                if rtr.stats.station == st_name[i] and rtr.stats.channel == cha:
                    if (max(abs(rtr.data))) < upper and (max(abs(rtr.data))) > lower:
                        # effectively take out variation due to event magnitude:
                        synnorm = (max(abs(rtr.data)))/(max(synmax)) #normalize by event maxima
                        # take out geometric spreading:
                        synspread = synnorm*rtr.stats.sac.dist
                        # take out baseline average:
                        synresidual = synspread-(sum(synspreadall)/len(synspreadall))
                        synall[channels.index(cha)][i] += synresidual
            if ii == 3 and rz_stacount[st_name[i]] != 0 and rr_stacount[st_name[i]] != 0 and rt_stacount[st_name[i]] != 0:
                synall[0][i] = (synall[0][i])/(rz_stacount[st_name[i]])
                synall[1][i] = (synall[1][i])/(rr_stacount[st_name[i]])
                synall[2][i] = (synall[2][i])/(rt_stacount[st_name[i]])

#%%% averages - data
    
    if ii == 3: 
        
        for i in range(len(channels)):
            xx = []
            yy = []
            
            for jj in range(len(st_name)):
                yy.append(st_lat[jj])
                xx.append(st_lon[jj])
                            
                            

            fig = plt.figure()
            ax = fig.add_subplot(111)
            cha = channels[i]

            # setup mercator map projection.
            m = Basemap(projection='merc',llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,epsg=4269, resolution='l') #http://server.arcgisonline.com/arcgis/rest/services; EPSG Number of America is 4269
            m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 2000, verbose= True)
        
            x = (xx)
            y = (yy)
            yi = np.linspace(llcrnrlat, urcrnrlat, 100)
            xi = np.linspace(llcrnrlon, urcrnrlon, 100)
            xi, yi = np.mgrid[-120:-115:100j, 32:35:100j]
            
            zi = scipy.interpolate.griddata((x,y), dataall[i], (xi,yi), method='linear')
            
            
            plt.contour(xi, yi, zi)
            plt.colorbar(orientation="horizontal",shrink=2/3)  # draw colorbar  

            plt.title("Average all events  - Data Peak Velocity for " + cha + " - " + str(1/high_freq) + "-" + str(1/low_freq) + " sec")
            plt.savefig(figPath + "Data - average all events " + cha + ".png")
            
#%%% averages - synthetics
        for i in range(len(channels)):
            xx = []
            yy = []
            
            for jj in range(len(st_name)):
                yy.append(st_lat[jj])
                xx.append(st_lon[jj])
        
            fig = plt.figure()
            ax = fig.add_subplot(111)
            cha = channels[i]
            
            # setup mercator map projection.
            m = Basemap(projection='merc',llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,epsg=4269, resolution='l') #http://server.arcgisonline.com/arcgis/rest/services; EPSG Number of America is 4269
            m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 2000, verbose= True)
            m.drawcoastlines()
            
            
            x = (xx)
            y = (yy)
            yi = np.linspace(llcrnrlat, urcrnrlat, 100)
            xi = np.linspace(llcrnrlon, urcrnrlon, 100)
            xi, yi = np.mgrid[-120:-115:100j, 32:35:100j]
            
            zi = scipy.interpolate.griddata((x,y), synall[i], (xi,yi), method='linear')
            
            
            plt.contour(xi, yi, zi) #, vmin=0, vmax=5)
            plt.colorbar(orientation="horizontal",shrink=2/3)  # draw colorbar  
            
            plt.title("Average all events  - Synthetic Peak Velocity for " + cha + " - " + str(1/high_freq) + "-" + str(1/low_freq) + " sec")
            # plt.show()
            plt.savefig(figPath + "Synthetic - average all events " + cha + ".png")
            
            
            