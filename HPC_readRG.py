#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 16:41:40 2020
@author: bcbirkel

Script for reading in Rob Grave's synthetics (FULL GRID)
Sorts binary file into spatial grid points (binary file is in timeslices),
saves points nearest existing stations to SAC files
"""

import struct
import numpy as np
import matplotlib.pyplot as plt
# import obspy.io.sac 
from obspy import Trace
from obspy import Stream
from obspy.geodetics.base import gps2dist_azimuth
import os
from obspy.core.utcdatetime import UTCDateTime

# read in binary file
hpc_fol = "/home1/birkel/"
#fileName = hpc_fol + "GravesSim_fullgrid_tsfiles/epw_102_m5.09-3.5x3.5-s266318098_cvmsi-vs500_sc01-h0.100/OutBin/sc01_xyts.e3d"
fileName = hpc_fol + "GravesSim_fullgrid_tsfiles/e1002_m5.39-5.0x5.0-s650146834_cvmh-vs500_sc01-h0.100/OutBin/sc01_xyts.e3d"
model = "CVM-H"
event = "chinohills"
#starttime = UTCDateTime("2014-03-29T04:09:42.994500Z")
starttime=UTCDateTime("2008-07-29T18:42:15.960Z")

# %% #### GET HEADER INFO #####

with open(fileName, mode='rb') as file: # rb-> read binary
    fileContent = file.read()
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


# %% ##### UNPACK MOTIONS FROM SIMULATION (IN TIMESLICES) #####

# define variables
fsize = nx*ny*3*nt # number of poimts
fpts = nx*ny # number of grid points

# unpack binary file
with open(fileName, mode='rb+') as f: # rb-> read binary
    fileContent = f.read()
    data = struct.unpack("f"*fsize, fileContent[0:fsize*4])
print("data read in")

##### SET UP GRID POINTS #####

# initialize variables
x = []
y = []
xy = []

# create array for grid points in x and y
a = np.arange(0,nx)
b = np.arange(0,ny)

# set up array with all grid points in order
# SQUARE portion of grid
for i in a:
    jloop = np.arange(i,-1,-1)
    for j in jloop:
        pta = [i*10, (i-j)*10]
        xy.append(pta)
        if j != 0:
            ptb = [(i-j)*10, i*10]
            xy.append(ptb)
            
# piece of grid that is NOT square:
rect = [i for i in b if i not in a]
for k in rect:
    lloop = np.arange(0,nx,1)
    for l in lloop:
        pta = [l*10, k*10]
        xy.append(pta)
print("xy built") #check
    
# %% ##### PULL IN TRANSFORMATION FROM X,Y TO LAT,LON #####
    
latlonmod = []
xymod = []
model_trans_fin = []

# pull in latlon to xy coordinate transformation info
modelcoord = hpc_fol + 'model_coords_sc01-h0.100'
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

xymod_fin = []
model_trans_fin = []
xymod_sorted = []

# skips to every 10th point in coord transform file to match sim grid
# xymod_skip = xymod[::10]
model_trans_skip = model_trans[::10]

#pulls out only useful points from coord transform file
for pt in model_trans_skip:
    if pt[0]%10 == 0 and pt[1]%10 == 0:
        model_trans_fin.append(pt)

#sorts points from coord transform file to match data ordering
for xypair in xy: #sorted correctly
    for pt in model_trans_fin:
        if xypair[0] == pt[0] and xypair[1] == pt[1]:
            xymod_sorted.append(pt)
            #print(pt)    
print("lat lon indices sorted")



# %% ##### PUT DATA IN TRACES W/ HEADER INFO #### 

st = Stream()
start = 15 
for i in range((fpts*3)):
    data_tr = data[start+i::fpts*3] # pull out data by spatial point (stored as timeslice in binary file)
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
    tr.stats.cmp = cmp
    tr.stats.dt = dt
    tr.stats.delta = 0.05
    tr.stats.starttime=starttime
    tr.data = np.asarray(data_tr) 
    st.append(tr)
print("data sorted into traces, added to stream")


# %% ##### FIND SIMULATION POINTS CLOSEST TO ACTUAL STATIONS #####

# read in station file
stationFile = hpc_fol + 'all_stationmaster.txt'
stations = []
savedpts = []

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

for sta in stations: # John's stations
    leng = []
    codept = []
    synpt = []
    for pt in model_trans_fin: #coords from transform file
        [dist,az,baz] = gps2dist_azimuth(pt[2], pt[3], sta[2], sta[3])
        leng.append(dist)
        # codept.append(sta[1])
        # synpt.append(pt)
    ind_min = np.argmin(leng)
    # print(ind_min)
    savedpts.append([model_trans_fin[ind_min],sta])

print("paired simulation points to stations")
    
# %% #####  PLOT STATIONS AND SAVED POINTS FROM GRID ####
# fig = plt.figure()
# ax = fig.add_subplot(111)

# for i in range(len(stations)):
#     ax.scatter(stations[i][2],stations[i][3],color='k')
# for j in range(len(savedpts)):
#     ax.scatter(savedpts[j][0][2],savedpts[j][0][3],color='r')    

# %% ##### WRITE EXTRACTED DATA CLOSEST TO STATION AS SAC FILE TO COMPARE TO OTHER STATION SYNTHETICS #####

# for all lat/lon in station file, search through latlon variable,
# find (x,y) point with min dist using gps2dist_azimuth, return index and assign
# that seismogram to that station. save.

if os.path.isdir(hpc_fol + 'extractedRGsyn/') == False:
    os.mkdir(hpc_fol + 'extractedRGsyn/')
os.chdir(hpc_fol + 'extractedRGsyn/')
for i in range(len(savedpts)):
    for tr in st:
        if tr.stats.xcoord == int(round(savedpts[i][0][0],-1)) and tr.stats.ycoord == int(round(savedpts[i][0][0],-1)):
            tr.stats.station = savedpts[i][1][1]
            tr.stats.channel = tr.stats.cmp
            tr.stats.network = savedpts[i][1][0]
            [dist,az,baz] = gps2dist_azimuth(savedpts[i][1][2], savedpts[i][1][3], tr.stats.lat, tr.stats.lon)
            tr.stats.location = dist
            print(dist)
            tr.write('Extracted_RGsyn_' + event + "_" + model + "_" + tr.stats.network + '_' + tr.stats.station + '_' + tr.stats.cmp + '.SAC',format='SAC')


    
    
    
    
    
    