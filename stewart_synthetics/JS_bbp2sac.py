#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 13:41:16 2022

@author: bcbirkel
"""


# Plot ASCII files to see Jon Stewart's synthetics
# units are cm/s
# Brianna Birkel - 11/2019
#     from __future__ import print_function

 
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from obspy import read_events, UTCDateTime, read, Trace, Stream
from obspy.geodetics import gps2dist_azimuth
from obspy.core.util.attribdict import AttribDict
from obspy.core.trace import Stats
from obspy.io.sac import SACTrace
import pandas as pd
import time
from scipy.spatial import distance
from pykdtree.kdtree import KDTree
import pickle
from obspy.io.sac import SACTrace
from obspy.geodetics.base import gps2dist_azimuth
import shutil


# %% VARIABLES
#    freqmin = 0.1
#    freqmax = 0.5
#    start_buff = 0
#    end_buff = 150
#    min_plot_dist = 0
#    max_plot_dist = 100
#    scale_fact = 1

vmodel = 'S'
eventfile = 'pw_102'
event_title = 'lahabra'

  

count = 0
# crap for file naming to work right
if vmodel == 'S' or vmodel == 'cvms426-223':
    mod = 'si'
    mod_title = 'S4'
elif vmodel == 'H' or vmodel == 'cvmhy':
    mod = 'h'
    mod_title = 'H'

if eventfile == 'pw_102':
    event = 'lahabra_2014'
    event_lat = 33.9325
    event_lon = -117.9158
    start = UTCDateTime("2014-03-29T04:09:42.994500Z") 
    
elif eventfile == '1002':
    event = 'chinohills_2008'
    event_lat = 33.9465
    event_lon = -117.7667
elif eventfile == '1011':
    event = 'inglewood_2009'
    event_lat = 33.9377
    event_lon = -118.3357
elif eventfile == '1019':
    event = 'chatsworth_2007'
    event_lat = 34.2983
    event_lon = -118.6255
elif eventfile == '1036':
    event = 'beverlyhills_2001'
    event_lat = 34.0541
    event_lon = -118.3929

# %% SETUP 
dir = '/Users/bcbirkel/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/LABasin/stewart_synthetics/CVM-' + mod_title + '/' + eventfile + '_' + event_title + '/Vel/'

print(dir)
os.chdir(dir)

event_files = os.listdir(dir)
num_files = len(event_files)
print(num_files) 
print(event_files)
#stastr = str(stations)
#print(stastr)

stT_all = Stream()
stE_all = Stream()
stN_all = Stream()
stZ_all = Stream()


#%% stations

# read in station file
stationFile = "../../../../all_stationmaster.txt"
stations = []
savedpts = []
stalist = []
tmpsta = []
ngasts = []

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
        stalist.append([net,code,lat,lon])
        # stations.append(code)


#Stewart syn - pull stations into list, mapped with station code and 

nga_stfile = '../../../Table 2 - List of Socal NGA-West2 Station Locations.txt'

staCoord = open(nga_stfile, 'r')
lines = staCoord.readlines()
for line in lines:
    if not line.startswith("Station"): # ignore first line
        split_line = line.split()
        seqno = split_line[0]
        lat = float(split_line[1])
        lon = float(split_line[2])
        ngasts.append([seqno,lat,lon])
        # stations.append(code)
        
        # %%## match up stations
        t0 = time.time()
        
        
        overlap = []
        XX = []  
        YY = []
        
        for sta in stalist: # all stations
            XX.append([sta[2], sta[3]])
        for pt in ngasts:  
            YY.append([pt[1], pt[2]])

        X = np.array(XX)
        Y = np.array(YY)
        
        tree = KDTree(Y)
        neighbor_dists, neighbor_indices = tree.query(X)
        
        for i in range(len(neighbor_indices)):
            if neighbor_dists[i] < 0.01:
                overlap.append([neighbor_dists[i],ngasts[neighbor_indices[i]],stalist[i]])
                
        t1 = time.time()
        
        print(overlap)
        print("Elapsed time finding closest points: " + str(t1-t0) + " seconds. Station: " + str(sta))
# %% create streams

                
stN = Stream()
stE = Stream()
stZ = Stream()   

for ii in range(num_files):
    fname = event_files[ii]
    
    for jj in range(len(overlap)):
        
        stastr="s" + str(overlap[jj][1][0]) + ".bbp"

        if stastr in fname:
            with open(fname, 'r') as file:
                time = []
                comp_N = []
                comp_E = []
                comp_Z = []
                all_text = '' 
                # Remove comments, put txt file into lists for time/comps
                num_lines = range(len(lines))
                #print('check2')
            
            
                # put data into list
                for line in open(fname):
                    if not line.startswith("#"):
                        all_text = all_text + line
                        line = line.split('\t')
                        #print(line)
#                            for num in range(0,4):
#                                line[num] = float(line[num])
                        time.append(line[0])
                        comp_N.append(line[1])
                        comp_E.append(line[2])
                        comp_Z.append(line[3])
                    if "longitude" in line:
                        line = line.split()
                        sta_lon = line[2]
                    if "latitude" in line:
                        line = line.split()
                        sta_lat = line[2]

# %% set time/component lists as NumPy arrays
                time = np.asarray(time)
                N_data = np.asarray(comp_N,dtype=float)
                E_data = np.asarray(comp_E,dtype=float)
                Z_data = np.asarray(comp_Z,dtype=float)
                
# %% fill header            
 
                # Fill header attributes
                print('station code going into header is ' + overlap[jj][2][1])
                
               # initialize traces
                trtmp = Trace()
                tr2 = SACTrace.from_obspy_trace(trtmp)
                tr = tr2.to_obspy_trace()
                trN = tr.copy()
                trE = tr.copy()
                trZ = tr.copy()
                
                streams = [stN,stE,stZ]
                traces = [trN,trE,trZ]
                data_arrs = [N_data,E_data,Z_data]
                cha = ["BHN","BHE","BHZ"]
                count = 0
                
                for tr in traces:    
                    tr.data = data_arrs[count]
                    tr.stats.delta = float(time[1])-float(time[0])
                    tr.stats.starttime = start
                    tr.stats.npts = len(time)
                    tr.stats.sampling_rate = 1/(float(time[1])-float(time[0]))
                    tr.stats.network = overlap[jj][2][0]
                    tr.stats.station = overlap[jj][2][1]
                    tr.stats.channel = cha[count]
                    tr.stats.sac.stla = overlap[jj][2][2]
                    tr.stats.sac.stlo = overlap[jj][2][3]
                    tr.stats.sac.evla = event_lat
                    tr.stats.sac.evla = event_lon
                    st = streams[count]
                    st.append(tr)
                    
                    [dist_m,az,baz] = gps2dist_azimuth(event_lat,event_lon,tr.stats.sac.stla,tr.stats.sac.stlo)
                    tr.stats.distance = dist_m
                    tr.stats.sac.dist = dist_m/1000
                    tr.stats.sac.az = az
                    tr.stats.sac.baz = baz
                
                    count = count + 1
                    
#%%### PLOT

stN.plot(type='section',recordlength=60,orientation='horizontal',offset_max=50000)
stE.plot(type='section',recordlength=60,orientation='horizontal',offset_max=50000)
stZ.plot(type='section',recordlength=60,orientation='horizontal',offset_max=50000)
    # %% convert to SAC and save files
outpath = ('/Users/bcbirkel/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/LABasin/CompiledEvents/' + event + "/JS_Syn/")
if os.path.exists(outpath):
    shutil.rmtree(outpath)
os.mkdir(outpath)
os.chdir(outpath)

st = stN + stE + stZ

for tr in st:
    tr.write("JS_Syn_" + event_title + "_" + tr.stats.network + "_" + tr.stats.station+ "_" + tr.stats.channel + ".SAC")



# print(os.getcwd())
# print("for loop ran " + str(count) + " times")
# print(stE_all)
# print(stE)
# print("check here")
# print(eventfile)
# #print(stN.stats.code)
# os.chdir('/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/compare_syn')
# #stT_all.write("time-" + eventfile + ".mseed", format='MSEED', reclen=256, ) 
# #stN_all.write("N_bp-" + eventfile +  ".mseed", format='MSEED', reclen=256)
# #stE_all.write("E_bp-" + eventfile +  ".mseed", format='MSEED', reclen=256)
# #stZ_all.write("Z_bp-" + eventfile +  ".mseed", format='MSEED', reclen=256)
        

#
