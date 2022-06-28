#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 16:33:36 2022

@author: bcbirkel
"""

import struct
from obspy import read
import numpy as np
import matplotlib.pyplot as plt
import obspy.io.sac 
from obspy import Trace
from obspy import Stream
from obspy.geodetics.base import gps2dist_azimuth
import os
from obspy.core.utcdatetime import UTCDateTime
import shutil
from mpl_toolkits.basemap import Basemap
import time
import math
from scipy.signal import argrelextrema
from obspy.signal.rotate import rotate_ne_rt
import csv

# %% ##### SET UP PATHS #####
write_st_syn = True
extract = True

try: event_chk
except NameError: event_chk = None

try: mod_chk
except NameError: mod_chk = None

lowf = 1/10
highf = 1/2

# iis = [0,1,2,3,4]
iis = [0]
ms = [0]

for kk in iis:
    for j in ms:
        event_no = kk
        events = ['lahabra', 'beverlyhills', 'chatsworth', 'chinohills', 'inglewood']
        eventnames = ['lahabra_2014', 'chatsworth_2007', 'chinohills_2008', 'inglewood_2009', 'beverlyhills_2001']
        event = events[event_no]
        eventname = eventnames[event_no]
        model_no = j
        models = ['CVM-S4', 'CVM-H']
        model = models[model_no]
        
        dataPath = './CompiledEvents/' + eventname + '/allData/fixed/'
        
        if event_no == 0:
            starttime = UTCDateTime("2014-03-29T04:09:42.994500Z")
            eventfol = "lahabra_2014/GravesSyn/" + model + "/"    
            eventlat = 33.9325
            eventlon = -117.9158
            if model_no == 0:
                fileName = "./timeslices/GravesSim_fullgrid_tsfiles/epw_102_m5.09-3.5x3.5-s266318098_cvmsi-vs500_sc01-h0.100/OutBin/sc01_xyts.e3d"
            if model_no == 1:
                fileName = "./timeslices/GravesSim_fullgrid_tsfiles/epw_102_m5.09-3.5x3.5-s266318098_cvmh-vs500_sc01-h0.100/OutBin/sc01_xyts.e3d"
        if event_no == 1:
            starttime = UTCDateTime("2001-09-09T23:59:17.695")
            eventfol = "beverlyhills_2001/GravesSyn/" + model + "/"
            eventlat = 34.0590
            eventlon = -118.3885 
            if model_no == 0:
                fileName = "./timeslices/GravesSim_fullgrid_tsfiles/e1036_m4.24-1.3x1.3-s1098915986_cvmsi-vs500_sc01-h0.100/OutBin/sc01_xyts.e3d"
            if model_no == 1:
                fileName = "./timeslices/GravesSim_fullgrid_tsfiles/e1036_m4.24-1.3x1.3-s1098915986_cvmh-vs500_sc01-h0.100/OutBin/sc01_xyts.e3d"
        if event_no == 2:
            starttime = UTCDateTime("2007-08-09T07:58:48.888 ")
            eventfol = "chatsworth_2007/GravesSyn/" + model + "/"
            eventlat = 34.2995 
            eventlon = -118.6195 
            if model_no == 0:
                fileName = "./timeslices/GravesSim_fullgrid_tsfiles/e1019_m4.66-2.1x2.1-s60148050_cvmsi-vs500_sc01-h0.100/OutBin/sc01_xyts.e3d"
            if model_no == 1:
                fileName = "./timeslices/GravesSim_fullgrid_tsfiles/e1019_m4.66-2.1x2.1-s60148050_cvmh-vs500_sc01-h0.100/OutBin/sc01_xyts.e3d"
        if event_no == 3:
            starttime = UTCDateTime("2008-07-29T18:42:15.960")
            eventfol = "chinohills_2008/GravesSyn/" + model + "/"
            eventlat = 33.9530
            eventlon = -117.7613 
            if model_no == 0:
                fileName = "./timeslices/GravesSim_fullgrid_tsfiles/e1002_m5.39-5.0x5.0-s650146834_cvmsi-vs500_sc01-h0.100/OutBin/sc01_xyts.e3d"
            if model_no == 1:
                fileName = "./timeslices/GravesSim_fullgrid_tsfiles/e1002_m5.39-5.0x5.0-s650146834_cvmh-vs500_sc01-h0.100/OutBin/sc01_xyts.e3d"
        if event_no == 4:
            starttime = UTCDateTime("2009-05-18T03:39:36.126")
            eventfol = "inglewood_2009/GravesSyn/" + model + "/"
            eventlat = 33.9377
            eventlon = -118.3357 
            if model_no == 0:
                fileName = "./timeslices/GravesSim_fullgrid_tsfiles/e1011_m4.70-2.2x2.2-s1247256210_cvmsi-vs500_sc01-h0.100/OutBin/sc01_xyts.e3d"
            if model_no == 1:
                fileName = "./timeslices/GravesSim_fullgrid_tsfiles/e1011_m4.70-2.2x2.2-s1247256210_cvmh-vs500_sc01-h0.100/OutBin/sc01_xyts.e3d"
        
        
        #set corners for map
        llcrnrlon=-119
        llcrnrlat=33.25 #- 0.25
        urcrnrlon=-117
        urcrnrlat=34.75 #+ 0.25
        
        
        # %% ##### GET DATA #####
        stream_dataZ = read(dataPath + "*Z_vel.SAC") 
        for tr in stream_dataZ:
            tr.stats.channel = 'Z'
        stream_dataR = read(dataPath + "*R_vel.SAC") 
        for tr in stream_dataR:
            tr.stats.channel = 'R'
        stream_dataT = read(dataPath + "*T_vel.SAC") 
        for tr in stream_dataT:
            tr.stats.channel = 'T'   
        
            # %% #### GET HEADER INFO #####
        if event_chk == event and mod_chk == model:
            extract = False
        
        if extract == True:
            with open(fileName, mode='rb') as file: # rb-> read binary
                fileContent = file.read(60)
                tshead = struct.unpack("iiiiiiiifffffff", fileContent[:60])
               
            # pull out header information    
            ix0 = tshead[0]            #/* starting x grid location for output */
            iy0 = tshead[1]            #/* starting y grid location for output */
            iz0 = tshead[2]            #/* starting z grid location for output */
            it0 = tshead[3]            #/* starting time step for output */
            nx = tshead[4]             #/* number of x points                                */
            ny = tshead[5]             #/* number of y points                                */
            nz = tshead[6]             #/* number of z points                                */
            nt = tshead[7]             #/* number of time points                                */
            dx = tshead[8]             #/* X direction spacing btw adjacent points        */
            dy = tshead[9]             #/* Y direction spacing btw adjacent points        */
            dz = tshead[10]            #/* Z direction spacing btw adjacent points        */
            dt = tshead[11]            #/* time step *
            modelrot = tshead[12]      #/* rotation of y-axis from south (clockwise positive)   */
            modellat = tshead[13]      #/* latitude of model origin                             */
            modellon = tshead[14]      #/* longitude of model origin                            */
            
    
                    
            # %% ##### MORE EFFICIENT UNPACK MOTIONS FROM SIMULATION BINARY FILE (IN TIMESLICES) #####
            
            # define variables
            fsize = nx*ny*3*nt # number of points
            spapts = nx*ny*3
            fpts = nx*ny # number of grid points
            
            t0 = time.time()
            ts_array = np.zeros((spapts,nt), dtype=float)
            count = 0
            
            def read_file_in_chunks(file_object, chunk_size=int(fsize/nt)):
                while True:
                    data = []
                    # Just read chunk_size size data.
                    f = file_object.read(chunk_size*4)
                    
                    if not f:
                        # Break the loop.
                        break
                    
                    data = struct.unpack("f"*chunk_size, f[0:chunk_size*4])
    
                    yield data
                    
            # Open the big data file.
            with open(fileName, mode='rb+') as f:
                # Invoke the above lazy file read data function.
                header = f.read(60)
                tshead = struct.unpack("iiiiiiiifffffff", fileContent[:60])
                print("Header: " + str(tshead))
                for tslice in read_file_in_chunks(f):
                    # Process the piece of data such as write the data to another file.
                    ts_array[:, count] = tslice
                    count += 1
                    if count % 1000 == 0:
                        t1 = time.time()
                        print("Elapsed time on data read: " + str(t1-t0) + " seconds. " + str(count) + " of 8000 timeslices done.")
                    
        
            # %% ##### PULL IN TRANSFORMATION FROM X,Y TO LAT,LON #####
                
            latlonmod = []
            xymod = []
            model_trans = []
            
            # pull in latlon to xy coordinate transformation info
            modelcoord = './model_coords_sc01-h0.100'
            coord = open(modelcoord, 'r')
            pairs = coord.readlines()
            for line in pairs:
                split_line = line.split()
                lon = float(split_line[0])
                lat = float(split_line[1])
                xmod = int(split_line[2])
                ymod = int(split_line[3])
                latlonmod.append([lat,lon])
                xymod.append([xmod,ymod])
                model_trans.append([xmod,ymod,lat,lon])
                
            print("coordinates stripped") #check
            
            # %% ##### SORT LAT,LON TO FIT ORDERING OF DATA FROM SIMULATION #####
            
            xymod_sorted = []
            
            # skips to every 10th point in coord transform file to match sim grid
            model_trans_skip = model_trans[::10]
            
            #pulls out only useful points from coord transform file
            for pt in model_trans_skip:
                if pt[0]%10 == 0 and pt[1]%10 == 0:
                    xymod_sorted.append(pt)
              
            print("lat lon indices sorted")
            
            # %% ##### PUT DATA IN TRACES W/ HEADER INFO #### 
            
            st = Stream()
            stN = Stream() 
            stE = Stream() 
            stZ = Stream() 
            
            dist_list = []
            # start = 15              # was 61 ?
            for i in range(3*fpts):
                # data_tr = data[start+i::3*fpts] # pull out data by spatial point (stored as timeslice in binary file)
                data_tr = ts_array[i,:]
                tr = Trace()
                if i < fpts:
                    cmp = 'N'
                    head = xymod_sorted[i]
                    tr.stats.xcoord = head[0]
                    tr.stats.ycoord = head[1]
                    tr.stats.lat = head[2]
                    tr.stats.lon = head[3]
                    
                if i >= fpts and i < 2*fpts:
                    cmp = 'E'
                    head = xymod_sorted[i-fpts]
                    tr.stats.xcoord = head[0]
                    tr.stats.ycoord = head[1]
                    tr.stats.lat = head[2]
                    tr.stats.lon = head[3]
                    
                if i >= 2*fpts and i < 3*fpts:
                    cmp = 'Z'
                    head = xymod_sorted[i-2*fpts]
                    tr.stats.xcoord = head[0]
                    tr.stats.ycoord = head[1]
                    tr.stats.lat = head[2]
                    tr.stats.lon = head[3]
                    
                [dist,az,baz] = gps2dist_azimuth(tr.stats.lat,tr.stats.lon,eventlat,eventlon)
                tr.stats.dist = dist/1000   
                tr.stats.baz = baz
                dist_list.append(dist)
                tr.stats.cmp = cmp
                tr.stats.dt = dt
                tr.stats.delta = 0.05
                tr.stats.starttime = starttime
                tr.data = np.asarray(data_tr) 
                # st.append(tr)
                if cmp == 'N':
                    stN.append(tr)
                if cmp == 'E':
                    stE.append(tr)
                if cmp == 'Z':
                    stZ.append(tr)
            print("data sorted into traces, added to stream")
            

            # %% ##### REMOVE POINTS OUT OF BOUNDS #####
            
            lats=[]; lons=[]
            for pt in xymod_sorted:
                lats.append(pt[2])
                lons.append(pt[3])
            maxlat = max(lats)
            minlat = min(lats)
            maxlon = max(lons)
            minlon = min(lons)
            
            # %% PLOT CHECK 
            
            fig = plt.figure(0, figsize=[8,5])
            a = fig.add_subplot(1, 1, 1)
                    
            # count = 0
            for tr in stN:
                if tr.stats.dist < 80:
                # if tr.stats.station == "13878":
                    # if count%300==0:
                    # plt.text(55,tr.stats.dist,tr.stats.station)
                    a.plot(tr.times(),tr.data+tr.stats.dist,color='k',lw=0.01)
                        
                    # print(tr.stats.station + " " + str(tr.stats.dist))
                # count = count+1
            plt.xlim(0,60)
            plt.ylim(0,80)
            plt.show()
        
            print("check plots")
            
            # %% DON'T EXTRACT AGAIN IF JUST EXTRACTED
            event_chk = event
            mod_chk = model
        
        # %% ##### IMPORT DATA

        stream_dataZ = read(dataPath + "*Z_vel.SAC") 
        for tr in stream_dataZ:
            tr.stats.channel = 'Z'
        stream_dataR = read(dataPath + "*R_vel.SAC") 
        for tr in stream_dataR:
            tr.stats.channel = 'R'
        stream_dataT = read(dataPath + "*T_vel.SAC") 
        for tr in stream_dataT:
            tr.stats.channel = 'T'  
        
        # %% ##### FIND SIMULATION POINTS CLOSEST TO ACTUAL STATIONS #####
        
        # # read in station file
        # stationFile = "./all_stationmaster.txt"
        # stations = []
        # savedpts = []
        
        # # strip station file
        # staCoord = open(stationFile, 'r')
        # lines = staCoord.readlines()
        # for line in lines:
        #     if not line.startswith("#"): # ignore first line
        #         split_line = line.split()
        #         net = split_line[0]
        #         code = split_line[1]
        #         lat = float(split_line[2])
        #         lon = float(split_line[3])
        #         if minlat < lat < maxlat and minlon < lon < maxlon:
        #             stations.append([net,code,lat,lon])
        
        # for sta in stations: # all stations
        #     leng = []
        #     for pt in xymod_sorted:    
        #         [dist,az,baz] = gps2dist_azimuth(pt[2], pt[3], sta[2], sta[3])
        #         leng.append(dist)
                
        #     ind_min = np.argmin(leng)
        #     savedpts.append([xymod_sorted[ind_min],sta])
            
            
        # %% ##### GET LINES FOR TRACING 
        
        linepts = []
        savedlnpts = []
        count = 0
        user_def = False
        
        if user_def == False:
            lat1 = []
            lon1 = []
            lat2 = []
            lon2 = []
            name = []
            dist_from_source = []
            with open("LABasin_BukaLines.csv", 'r') as file:
                csvreader = csv.reader(file)
                header = next(csvreader)
                for row in csvreader:
                    lat1.append(float(row[1]))
                    lon1.append(float(row[0]))
                    lat2.append(float(row[3]))
                    lon2.append(float(row[2]))
                    name.append(str(row[4]))
                    [distd,azd,bazd] = gps2dist_azimuth(eventlat, eventlon, np.mean([float(row[1]),float(row[3])]), np.mean([float(row[0]),float(row[2])]))
                    dist_from_source.append(distd)
        
        
        if user_def == True:
            lat1 = float(input("Latitude of beginning point: "))
            lon1 = float(input("Longitude of beginning point: "))
            lat2 = float(input("Latitude of end point: "))
            lon2 = float(input("Longitude of end point: "))
        
        
        
        # %% ##### PULL OUT POINTS CLOSEST TO LINE(S)
        
        t0 = time.time()
        
        for ii in range(len(lat1)):
            savedlnpts = []
            [dist,az,baz] = gps2dist_azimuth(lat1[ii],lon1[ii],lat2[ii],lon2[ii])
            
            def intermediates(p1, p2, nb_points=int(dist/50)):
                # If we have 8 intermediate points, we have 8+1=9 spaces
                # between p1 and p2
                x_spacing = (p2[0] - p1[0]) / (nb_points + 1)
                y_spacing = (p2[1] - p1[1]) / (nb_points + 1)
            
                return [[p1[0] + i * x_spacing, p1[1] +  i * y_spacing] 
                        for i in range(1, nb_points+1)]
            
            linepts = intermediates([lat1[ii], lon1[ii]], [lat2[ii], lon2[ii]], nb_points=int(dist/500))

            
            for lnpt in linepts: # points along line
                len_arr = np.zeros(len(xymod_sorted), dtype=float)
                ind = 0
                for pt in xymod_sorted:    
                    len_arr[ind] = np.sqrt(np.sum(np.square(np.array((pt[2], pt[3])) - np.array((lnpt[0], lnpt[1])))))
                    ind += 1
                    
                ind_min = np.argmin(len_arr)
                # ind_min = np.argmin(leng)
                savedlnpts.append([xymod_sorted[ind_min]])
                count += 1
                if count % int(10) == 0:
                    t1 = time.time()
                    print("Elapsed time for finding closest line points: " + str(t1-t0) + " seconds. Currently on coordinate set: " + str(ii) + ".")
                
                
            ##### REMOVE DUPLICATES    
            tmp = sorted(savedlnpts)
            savedlnpts = [tmp[i] for i in range(len(tmp)) if i == 0 or tmp[i] != tmp[i-1]]    
            
            
            #####  PLOT STATIONS AND SAVED POINTS FROM GRID ####
            fig = plt.figure()
            ax = fig.add_subplot(111)
            
            for j in range(len(savedlnpts)):
                ax.scatter(savedlnpts[j][0][3],savedlnpts[j][0][2], color='b')    
                
            
           ##### SET UP SYNTHETIC STREAMS FOR ROTATED COMPONENTS
            
            stR = Stream()
            stT = Stream()
            
            for i in range(len(stN)):
                    if stN[i].stats.xcoord == stE[i].stats.xcoord and stN[i].stats.ycoord == stE[i].stats.ycoord:
                        for pt in savedlnpts:
                                if stN[i].stats.xcoord == pt[0][0] and stN[i].stats.ycoord == pt[0][1]:
                                    trR = stN[i]
                                    trT = stE[i]
                                    trR.stats.cmp = 'R'
                                    trT.stats.cmp = 'T'
                                    trR.data, trT.data = rotate_ne_rt(stN[i].data,stE[i].data,stN[i].stats.baz)
                                    trR.filter("bandpass", freqmin=lowf, freqmax=highf)
                                    trT.filter("bandpass", freqmin=lowf, freqmax=highf)
                                    stR.append(trR)
                                    stT.append(trT)
    
                                    # print("data rotated")
            # %% ##### FIND DATA STATIONS IN REGION OF INTEREST in R     
            saved_dataptsR = []
            
            for lnpt in linepts: # points along line
                len_arr = np.zeros(len(stream_dataR), dtype=float)
                ind_data = 0
                for tr in stream_dataR:    
                    len_arr[ind_data] = np.sqrt(np.sum(np.square(np.array((tr.stats.sac.stlo, tr.stats.sac.stla)) - np.array((lnpt[0], lnpt[1])))))
                    ind_data += 1
            
                    
                ind_min = np.argmin(len_arr)
                saved_dataptsR.append([stream_dataR[ind_min]])
                count += 1
                if count % int(10) == 0:
                    t1 = time.time()
                    print("Elapsed time for finding closest line points: " + str(t1-t0) + " seconds. Currently on coordinate set: " + str(ii) + ".")
                     
            # %% ##### FIND DATA STATIONS IN REGION OF INTEREST in T    
            saved_dataptsT = []
            
            for lnpt in linepts: # points along line
                len_arr = np.zeros(len(stream_dataR), dtype=float)
                ind_data = 0
                for tr in stream_dataT:    
                    len_arr[ind_data] = np.sqrt(np.sum(np.square(np.array((tr.stats.sac.stlo, tr.stats.sac.stla)) - np.array((lnpt[0], lnpt[1])))))
                    ind_data += 1
            
                    
                ind_min = np.argmin(len_arr)
                saved_dataptsT.append([stream_dataR[ind_min]])
                count += 1
                if count % int(10) == 0:
                    t1 = time.time()
                    print("Elapsed time for finding closest line points: " + str(t1-t0) + " seconds. Currently on coordinate set: " + str(ii) + ".")
                     
                         
            # %% ##### PLOT                  
            
            ylim_min = np.min(dist_from_source)/1000 - 5
            ylim_max = np.max(dist_from_source)/1000 + 5
            
            fig = plt.figure(0, figsize=[8,5])
            a = fig.add_subplot(1, 1, 1)
                                    
                        
            for tr in stR:
                a.plot(tr.times(),tr.data+tr.stats.dist,color='k',lw=0.2)
                
            for kk in range(len(saved_dataptsR)):
                tr = saved_dataptsR[kk][0]
                a.plot(tr.times(),tr.data+tr.stats.sac.dist,color='b',lw=0.2)
            
            plt.xlim(0,100)
            plt.ylim(ylim_min,ylim_max)
            plt.xlabel("Time (seconds)")
            plt.ylabel("Distance from starting point of line (km)")
            plt.title("Radial - " + name[ii] + " - " + event + "," + model)
            plt.show()
            
            fig = plt.figure(0, figsize=[8,5])
            a = fig.add_subplot(1, 1, 1)
            
            # for tr in stT:
            #     linedist = gps2dist_azimuth(lat1, lon1, tr.stats.lat, tr.stats.lon)
            #     a.plot(tr.times(),tr.data+linedist/1000,color='k',lw=0.2)
            
            for tr in stT:
                a.plot(tr.times(),tr.data+tr.stats.dist,color='k',lw=0.2)
                
            for kk in range(len(saved_dataptsT)):
                tr = saved_dataptsT[kk][0]
                a.plot(tr.times(),tr.data+tr.stats.sac.dist,color='b',lw=0.2)
    
            plt.xlim(0,100)
            plt.ylim(ylim_min,ylim_max)
            plt.xlabel("Time (seconds)")
            plt.ylabel("Distance from starting point of line (km)")
            plt.title("Tranverse - " + name[ii] + " - " + event + "," + model)
            plt.show()
            
            print("check plots")
        
        # %% ##### WRITE EXTRACTED DATA CLOSEST TO STATION AS SAC FILE TO COMPARE TO OTHER STATION SYNTHETICS #####
        
        # # for all lat/lon in station file, search through latlon variable,
        # # find (x,y) point with min dist using gps2dist_azimuth, return index and assign
        # # that seismogram to that station. save.
        # checklist = []
        
        # os.chdir('./CompiledEvents/' + eventfol)
        # for i in range(len(savedpts)):
        #     for tr in st:
        #         if tr.stats.xcoord == savedpts[i][0][0] and tr.stats.ycoord == savedpts[i][0][1]:
        #             tr.stats.network = savedpts[i][1][0]
        #             tr.stats.station = savedpts[i][1][1]
        #             tr.stats.channel = tr.stats.cmp
        #             [dist,az,baz] = gps2dist_azimuth(eventlat, eventlon, tr.stats.lat, tr.stats.lon)
        #             tr.stats.dist = dist/1000
        #             tr.stats.starttime = tr.stats.starttime-5
        #             np.pad(tr.data,(int(20*4.2),0),'constant',constant_values=(0))
        #             if tr.stats.channel == 'N':
        #                 checklist.append([tr.stats.network, tr.stats.station, tr.stats.dist])
                    
        #             # print(dist)
        #             if write_st_syn == True:
        #                 tr.write('Extracted_RGsyn_' + event + "_" + model + "_" + tr.stats.network + '_' + tr.stats.station + '_' + tr.stats.cmp + '.SAC',format='SAC')
        # checklist.sort(key = lambda x: x[2])
        # print(checklist)
        
        # # %% ##### PULL OUT BEAMFORMING SQUARE(S) ##### 
        
        # beamfol = './CompiledEvents/' + eventfol + 'beamforming/'
        
        # if os.path.exists(beamfol):
        #     shutil.rmtree(beamfol)
        # os.mkdir(beamfol)
                
        # for tr in st:
        #     if not tr.stats.station == "":
        #         # print(tr.stats.station)
        #         if tr.stats.station == "WTT2" or tr.stats.station == "WTT":
        #             x_orig = tr.stats.xcoord - 100
        #             y_orig = tr.stats.ycoord - 100
        #             print("WTT2 station location: " + str(tr.stats.xcoord) + ", " + str(tr.stats.ycoord))
        
                
        # square = []
        # for i in range(20):
        #     for j in range(20):
        #         square.append([x_orig+j*10,y_orig+i*10])
            
        # beamst = st.copy()
        
        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        # fig2 = plt.figure()
        # ax2 = fig2.add_subplot(111)
        # m = Basemap(projection='merc',llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,epsg=4269) #http://server.arcgisonline.com/arcgis/rest/services; EPSG Number of America is 4269
        # m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 2000, verbose= True)
        
        # for pt in square:
        #     for tr in beamst:
        #         if tr.stats.xcoord == pt[0] and tr.stats.ycoord == pt[1]:
        #             [dist,az,baz] = gps2dist_azimuth(eventlat, eventlon, tr.stats.lat, tr.stats.lon)
        #             tr.stats.dist = dist/1000
        #             tr.stats.network = 'BEA'
        #             tr.stats.station = str(pt[0]-x_orig) + str(pt[1]-y_orig)
        #             tr.stats.sac = obspy.core.util.attribdict.AttribDict()
        #             tr.stats.sac.dist = dist/1000
        #             tr.stats.sac.az = az
        #             tr.stats.sac.baz = baz
        #             tr.stats.sac.stla = tr.stats.lat
        #             tr.stats.sac.stlo = tr.stats.lon
        #             tr.stats.sac.evla = eventlat
        #             tr.stats.sac.evlo = eventlon
        #             tr.stats.sac.user0 = tr.stats.xcoord
        #             tr.stats.sac.user1 = tr.stats.ycoord
        #             tr.stats.sac.user2 = tr.stats.xcoord-x_orig
        #             tr.stats.sac.user3 = tr.stats.ycoord-y_orig
        #             # print("BEA trace stat check: " + tr.stats)
        #             ax.plot(tr.times(),tr.data*10+tr.stats.dist,color='k',lw=0.2) 
        #             ax.set_xlim(0,60)
        #             xx,yy = m(tr.stats.lon,tr.stats.lat)
        #             m.scatter(xx,yy, marker = "o" ,s=0.5, color='k', alpha = 1)
        #             tr.stats.starttime = starttime
        #             # np.pad(tr.data,(int(20*4.2),0),'constant',constant_values=(0)) ## PADDING OCCURS DURING ROTATION, NOT HERE.
        #             tr.write('Beam_Grid_' + event + "_" + model + "_" + str(pt[0]-x_orig) + "," + str(pt[1]-y_orig) + '_' + tr.stats.cmp + '.SAC',format='SAC')
                    
        # xx,yy = m(eventlon,eventlat)
        # m.scatter(xx,yy, marker = "*" ,s=100, color='m', edgecolors = "k", alpha = 1)
        # plt.title(event + ' Beamforming Grid')
        # plt.show()
           
        # # for tr in st:
        # #     if xymod_sorted[i][0] == x_orig and savedpts[i][1] == y_orig:
        # #         tr.stats.channel = tr.stats.cmp
        # #         [dist,az,baz] = gps2dist_azimuth(eventlat, eventlon, tr.stats.lat, tr.stats.lon)
        # #         tr.stats.dist = dist/1000
        # #         tr.stats.starttime = tr.stats.starttime-0.6
        
                
        # #         # print(dist)
        # #         tr.write('Extracted_RGsyn_' + event + "_" + model + "_" + tr.stats.network + '_' + tr.stats.station + '_' + tr.stats.cmp + '.SAC',format='SAC')
        # # checklist.sort(key = lambda x: x[2])
        # # print(checklist)
        
        
        # # %% PLOT CHECK
        # fig = plt.figure()
        # a = fig.add_subplot(111)
        # for i in range(len(savedpts)):
        #     for tr in st:
        #         # if tr.stats.xcoord == int(round(savedpts[i][0][0],-1)) and tr.stats.ycoord == int(round(savedpts[i][0][1],-1)) and tr.stats.channel == 'N':
        #         if tr.stats.xcoord == savedpts[i][0][0] and tr.stats.ycoord == savedpts[i][0][1] and tr.stats.channel == 'N':
        #             if tr.stats.dist < 50:
        #             # if count%300==0:
        #                 a.text(55,tr.stats.dist,tr.stats.station)
        #                 a.plot(tr.times(),tr.data*10+tr.stats.dist,color='k',lw=0.2)    
        #             # print(tr.stats.station + " " + str(tr.stats.dist))
        # plt.xlim(0,60)
        # plt.ylim(0,50)
        # plt.show()
            
            
            
            
            
            