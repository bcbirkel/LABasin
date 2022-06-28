#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 11:18:31 2020

@author: bcbirkel
"""

# This code reads in data from CESMD and converts it to SAC with appropriate headers. All CESMD files must be
# processed by CESMD (could work on raw, haven't checked), with all V2 files in a single directory. Choose 
# a different directory (at the bottom) for the SAC files to be written into.

from obspy import read, read_inventory
from obspy.signal.rotate import rotate2zne
import matplotlib.pyplot as plt
import numpy as np
import os
from obspy.core import Trace, Stream, AttribDict
import math
from obspy.core import UTCDateTime
from obspy.io.sac import SACTrace
from zipfile import ZipFile as unzip

#Ignore warnings due to python 2 and 3 conflict
import warnings
warnings.filterwarnings("ignore")


# %% Set parameters
#set event
events = ['lahabra_2014', 'beverlyhills_2001', 'chatsworth_2007', 'chinohills_2008', 'inglewood_2009', 'elmonte_2020']

event = events[1]

angled = Stream()

#set corners for map
llcrnrlon=-119
llcrnrlat=33
urcrnrlon=-117
urcrnrlat=34.5

#set bandpass frequencies
high_freq = 1/1
low_freq = 1/10

# %% Set paths
#set directories

path_dir_CESMD = '/Users/bcbirkel/Documents/GitHub/LABasin/CompiledEvents/' + event + '/CESMD/'
os.chdir(path_dir_CESMD)

v2fol = path_dir_CESMD + "V2/"
if os.path.isdir(v2fol) == True:
    os.rmdir(v2fol)
if os.path.isdir(v2fol) == False:
    os.mkdir(v2fol)

for root, dirs, files in os.walk((path_dir_CESMD), topdown=True):
   for name in files:
       if name.endswith(".zip") == True or name.endswith(".ZIP") == True:
           zipfol = os.path.join(root, name)
           with unzip(zipfol,'r') as zfol:
               zfol.extractall()
               for file in os.listdir("."):
                    if file.endswith(".V2") == True:
                        print(file)
                        os.rename(file, v2fol + file)



    # for fol in os.listdir(path_dir_CESMD):
    #     if fol.startswith(".") == False:
    #         print(fol)
    #         os.chdir(path_dir_CESMD + fol)
    #         for fol2 in os.listdir(path_dir_CESMD + fol):
    #             os.chdir(path_dir_CESMD + fol + fol2)
    #             for zipfol in os.listdir(path_dir_CESMD + fol):
    #                 print(zipfol)
    #                 if zipfol.endswith(".zip") == True or zipfol.endswith(".ZIP") == True:
    #                     with unzip(zipfol,'r') as zfol:
    #                         zfol.extractall()
    #                         print(os.listdir(path_dir_CESMD + fol))
    #                         for file in os.listdir(path_dir_CESMD + fol):
    #                             if file.endswith(".V2") == True:
    #                                 print(file)
    #                                 os.rename(file, v2fol + file)
            
# %%

if os.path.isdir(path_dir_CESMD + 'SAC/') == False:
    os.mkdir(path_dir_CESMD + 'SAC/')
    os.mkdir(path_dir_CESMD + 'SAC/acc/')
    os.mkdir(path_dir_CESMD + 'SAC/vel/')
    os.mkdir(path_dir_CESMD + 'SAC/disp/')

#read in files
file_list = os.listdir(v2fol)

#CESMD_stream = Stream()
for file in file_list:
    if file.startswith("CHAN") == False:
        os.chdir(v2fol)
        fname = open(file, "r")
        data = fname.read()
        lines = data.split('\n')
        
        deg = lines[0].split()[5]
        
        station = lines[5].split()[2]
        stlat = float(lines[5].split()[3][:-2])
        stlon = float(lines[5].split()[4][:-1])
        date = lines[4].split()[3].split('/')
        year = '20' + date[2].split(',')[0]
        month = date[0]
        day = date[1]
        time = lines[4].split()[4]
        datetime = year + "-" + month + "-" + day + "T" +  time
        starttime = UTCDateTime(datetime)
        
        tr = Trace()
        tr.stats.sac = AttribDict()
        tr.stats.npts = float(lines[45].split()[0])
        tr.stats.dt = float(lines[45].split()[8])
        
        
        nxt = 45 #beginning of first set of data points
        
        for ii in range(9):
            tr = Trace()
            npts = float(lines[nxt].split()[0])
            dt = float(lines[nxt].split()[8])
            datapts = []
            for line in lines[nxt+1:]:
                if len(datapts) < npts:
                    pts = line.replace('-',' -').split()
                    for i in range(len(pts)):
                        datapts.append(float(pts[i]))
                tr.data = np.asarray(datapts)
            
            tr.stats.sac = AttribDict()
            tr.stats.station = station
            tr.stats.network = "CE"
            tr.stats.delta = dt
            tr.stats.npts = npts
            
            tr.stats.sac.kstnm = station
            tr.stats.sac.kntwk = "CE"
            tr.stats.sac.stla = stlat
            tr.stats.sac.stlo = stlon
            tr.stats.sac.delta = dt
            tr.stats.sac.npts = npts
            tr.stats.sac.leven = True
            tr.stats.starttime = starttime 
            
            if ii%3 == 0:
                tr.stats.sac.kuser0 = "acc"
                tr.stats.sac.kuser1 = "cm/s^2"
                if ii == 3 or ii == 6:
                    deg = lines[nxt-45].split()[5]
                    print("current angle: " + deg)
            if ii%3 == 1:
                tr.stats.sac.kuser0 = "vel"
                tr.stats.sac.kuser1 = "cm/s"
            if ii%3 == 2:
                nxt = nxt + 46
                tr.stats.sac.kuser0 = "disp"
                tr.stats.sac.kuser1 = "cm"
                
            os.chdir('/Users/bcbirkel/Documents/GitHub/LABasin/CompiledEvents/' + event + '/CESMD/SAC/')
           
            if deg == "Up":
                    tr.stats.channel = "BHZ"
                    tr.stats.sac.cmp = "BHZ"
                    print("writing vertical comp")
                    tr.write(tr.stats.network + "_" + tr.stats.station + "_" + tr.stats.sac.cmp + "_" + str(tr.stats.sac.kuser0) + ".SAC", format="SAC")
            elif deg == "90":
                    tr.stats.channel = "BHE"
                    tr.stats.sac.cmp = "BHE"
                    print("writing East comp")
                    tr.write(tr.stats.network + "_" + tr.stats.station + "_" + tr.stats.sac.cmp + "_" + str(tr.stats.sac.kuser0) + ".SAC", format="SAC")
            elif deg == "360" or deg == "0":
                    tr.stats.channel = "BHN"
                    tr.stats.sac.cmp = "BHN"
                    print("writing North comp")
                    tr.write(tr.stats.network + "_" + tr.stats.station + "_" + tr.stats.sac.cmp + "_" + str(tr.stats.sac.kuser0) + ".SAC", format="SAC")
            elif deg == "H1" or deg == "H2":
                print("az unknown")
            else:
                tr.stats.cmpaz = deg
                print("odd angle: " + deg)
                angled.append(tr)   
            
            print(lines[nxt])            
            nxt = nxt + math.ceil(npts/8) + 1
    
for file in os.listdir(path_dir_CESMD + 'SAC/'):
    if file.endswith("acc.SAC"):
        os.rename(file, path_dir_CESMD + 'SAC/acc/' + file)
    if file.endswith("vel.SAC"):
        os.rename(file, path_dir_CESMD + 'SAC/vel/' + file)
    if file.endswith("disp.SAC"):
        os.rename(file, path_dir_CESMD + 'SAC/disp/' + file)

# os.system('say "your program is done"')

# os.mkdir(path_dir_CESMD + 'SAC/allSACfiles/')
# for file in os.listdir(path_dir_CESMD + 'SAC/'):
#     if file.endswith(".SAC"):
#         os.rename(file, path_dir_CESMD + 'SAC/allSACfiles/' + file)
        

# stns = []
# for tr in angled:
#     stns += tr.stats.station

# for tr in angled:
#     for st in stns:
#         if st == tr.stats.station:
            
    

# path_dir_CESMD_odd = '/Users/bcbirkel/Documents/GitHub/LABasin/CompiledEvents/' + event + '/CESMD/SAC/'
# stream_odd = read(path_dir_CESMD_odd + "*NOAZ*.SAC")

# for tr in stream_odd:
#     sta = tr.stats.station
#     angled.append(sta)

# angled_sta = ['13079', '14840', '23773', '24363', '24983']
# for file in file_list:
#     os.chdir(path_dir_CESMD)
#     fname = open(file, "r")
#     data = fname.read()
#     lines = data.split('\n')
    
#     deg = lines[0].split()[5]
#     station = lines[5].split()[2]
#     for st in angled_sta:
#         if station = sta:
            
            
            
#     for jj in range(9):
#         if tr.stats.station == sta:
#             if c == ch1:
#                 d1 = ch1.vel.data
#                 a1 = ch1.stats.az
#                 dip1 = 0
#                 if a1 == 'Up':
#                     a1 = 0
#                     dip1 = -90
#             if c == ch2:
#                 d2 = ch2.vel.data
#                 a2 = ch2.stats.az
#                 dip2 = 0
#                 if a2 == 'Up':
#                     a2 = 0
#                     dip2 = -90
#             if c == ch3:
#                 d3 = ch3.vel.data
#                 a3 = ch3.stats.az
#                 dip3 = 0
#                 if a3 == 'Up':
#                     a3 = 0
#                     dip3 = -90
#     for c in st:
#         rotated = rotate2zne(d1,a1,dip1,d2,a2,dip2,d3,a3,dip3)
#     print("rotated:" + rotated)
        
    
    # tr2 = Trace()
    # datapts = []
    # nxt = 46 + math.ceil(tr.stats.num_pts/8)
    # tr2.stats.num_pts = float(lines[nxt].split()[0])
    # tr2.stats.dt = float(lines[nxt].split()[8])
    # for line in lines[nxt+1:]:
    #     if len(datapts) < tr.stats.num_pts:
    #         pts = line.replace('-',' -').split()
    #         for i in range(len(pts)):
    #             datapts.append(float(pts[i]))
    #     tr2.data = np.asarray(datapts)
    # ch1.vel = tr2
    
    # tr3 = Trace()
    # datapts = []
    # nxt = nxt + math.ceil(tr2.stats.num_pts/8) + 1
    # tr3.stats.num_pts = float(lines[nxt].split()[0])
    # tr3.stats.dt = float(lines[nxt].split()[8])
    # for line in lines[nxt+1:]:
    #     if len(datapts) < tr.stats.num_pts:
    #         pts = line.replace('-',' -').split()
    #         for i in range(len(pts)):
    #             datapts.append(float(pts[i]))
    #     tr3.data = np.asarray(datapts)
    # ch1.disp = tr3
    

    # s2 = nxt + math.ceil(tr3.stats.num_pts/8) + 2
    
    # tr4= Trace()
    # ch2.stats.az = lines[s2+0].split()[5]
    # ch2.stats.station = lines[s2+5].split()[2]
    # ch2.stats.stlat = float(lines[s2+5].split()[3][:-2])
    # ch2.stats.stlon = float(lines[s2+5].split()[4][:-1])
    # ch2.stats.starttime = UTCDateTime(datetime)
    # tr4.stats.num_pts = float(lines[s2+45].split()[0])
    # tr4.stats.dt = float(lines[s2+45].split()[8])
    
    # datapts = []
    # for line in lines[s2+46:]:
    #     if len(datapts) < tr.stats.num_pts:
    #         pts = line.replace('-',' -').split()
    #         for i in range(len(pts)):
    #             datapts.append(float(pts[i]))
    #     tr4.data = np.asarray(datapts)
    # ch2.acc = tr4
    
    # tr5 = Trace()
    # datapts = []
    # nxt = s2 + 46 + math.ceil(tr4.stats.num_pts/8)
    # tr5.stats.num_pts = float(lines[nxt].split()[0])
    # tr5.stats.dt = float(lines[nxt].split()[8])
    # for line in lines[nxt+1:]:
    #     if len(datapts) < tr.stats.num_pts:
    #         pts = line.replace('-',' -').split()
    #         for i in range(len(pts)):
    #             datapts.append(float(pts[i]))
    #     tr5.data = np.asarray(datapts)
    # ch2.vel = tr5
    
    # tr6 = Trace()
    # datapts = []
    # nxt = nxt + math.ceil(tr5.stats.num_pts/8) + 1
    # tr6.stats.num_pts = float(lines[nxt].split()[0])
    # tr6.stats.dt = float(lines[nxt].split()[8])
    # for line in lines[nxt+1:]:
    #     if len(datapts) < tr.stats.num_pts:
    #         pts = line.replace('-',' -').split()
    #         for i in range(len(pts)):
    #             datapts.append(float(pts[i]))
    #     tr6.data = np.asarray(datapts)
    # ch2.disp = tr6

    # s3 = nxt + math.ceil(tr3.stats.num_pts/8) + 2
    
    # tr7= Trace()
    # ch3.stats.az = lines[s3+0].split()[5]
    # ch3.stats.station = lines[s3+5].split()[2]
    # ch3.stats.stlat = float(lines[s3+5].split()[3][:-2])
    # ch3.stats.stlon = float(lines[s3+5].split()[4][:-1])
    # ch3.stats.starttime = UTCDateTime(datetime)
    # tr7.stats.num_pts = float(lines[s3+45].split()[0])
    # tr7.stats.dt = float(lines[s3+45].split()[8])
    
    # datapts = []
    # for line in lines[s3+46:]:
    #     if len(datapts) < tr7.stats.num_pts:
    #         pts = line.replace('-',' -').split()
    #         for i in range(len(pts)):
    #             datapts.append(float(pts[i]))
    #     tr7.data = np.asarray(datapts)
    # ch3.acc = tr7
    
    # tr8 = Trace()
    # datapts = []
    # nxt = s3 + 46 + math.ceil(tr7.stats.num_pts/8)
    # tr8.stats.num_pts = float(lines[nxt].split()[0])
    # tr8.stats.dt = float(lines[nxt].split()[8])
    # for line in lines[nxt+1:]:
    #     if len(datapts) < tr8.stats.num_pts:
    #         pts = line.replace('-',' -').split()
    #         for i in range(len(pts)):
    #             datapts.append(float(pts[i]))
    #     tr8.data = np.asarray(datapts)
    # ch3.vel = tr8
    
    # tr9 = Trace()
    # datapts = []
    # nxt = nxt + math.ceil(tr8.stats.num_pts/8) + 1
    # tr9.stats.num_pts = float(lines[nxt].split()[0])
    # tr9.stats.dt = float(lines[nxt].split()[8])
    # for line in lines[nxt+1:]:
    #     if len(datapts) < tr9.stats.num_pts:
    #         pts = line.replace('-',' -').split()
    #         for i in range(len(pts)):
    #             datapts.append(float(pts[i]))
    #     tr9.data = np.asarray(datapts)
    # ch3.disp = tr9
    
    # st = [ch1, ch2, ch3]
    # for ch in st:
    #     if ch.stats.az == "Up" or ch.stats.az == "BHZ":
    #         ch.stats.cmp = "BHZ"  
    #         ch.stats.az = 0
    #     elif ch.stats.az == "360" or ch.stats.az == "BHN":
    #         ch.stats.cmp = "BHN"
    #         ch.stats.az = 0
    #     elif ch.stats.az == "90" or ch.stats.az == "BHE":
    #         ch.stats.cmp = "BHE"
    #         ch.stats.az = 90
    #     elif ch.stats.az == 0 or ch.stats.az == 90:
    #         print("no change to az or cmp")
    #     elif ch.stats.az == "H1" or ch.stats.az == "H2":
    #         print("az unknown")
    #     else:
    #         sta = ch.stats.station
    #         for c in st:
    #             if c.stats.station == sta:
    #                 if c == ch1:
    #                     d1 = ch1.vel.data
    #                     a1 = ch1.stats.az
    #                     dip1 = 0
    #                     if a1 == 'Up':
    #                         a1 = 0
    #                         dip1 = -90
    #                 if c == ch2:
    #                     d2 = ch2.vel.data
    #                     a2 = ch2.stats.az
    #                     dip2 = 0
    #                     if a2 == 'Up':
    #                         a2 = 0
    #                         dip2 = -90
    #                 if c == ch3:
    #                     d3 = ch3.vel.data
    #                     a3 = ch3.stats.az
    #                     dip3 = 0
    #                     if a3 == 'Up':
    #                         a3 = 0
    #                         dip3 = -90
    #         for c in st:
    #             rotated = rotate2zne(d1,a1,dip1,d2,a2,dip2,d3,a3,dip3)
    #         print("rotated:" + rotated)
        
    
    #     tr.stats.sac = AttribDict()
        
    #     #     ch.rotate2zne("RT->NE", back_azimuth=ch.stats.baz)
    
    # os.chdir('/Users/bcbirkel/Documents/GitHub/LABasin/CompiledEvents/' + event + '/CESMD/SAC/')
    
    # # ch1.write(ch1.stats.station + '_' + ch1.stats.cmp + "_all.SAC", format="SAC") 
    # # ch2.write(ch2.stats.station + '_' + ch2.stats.cmp + "_all.SAC", format="SAC") 
    # # ch3.write(ch3.stats.station + '_' + ch3.stats.cmp + "_all.SAC", format="SAC") 
    
    # # ch1.vel.write(ch1.stats.station + '_' + ch1.stats.cmp + "_vel.SAC", format="SAC") 
    # # ch2.vel.write(ch2.stats.station + '_' + ch2.stats.cmp + "_vel.SAC", format="SAC") 
    # # ch3.vel.write(ch3.stats.station + '_' + ch3.stats.cmp + "_vel.SAC", format="SAC")     
