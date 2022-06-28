#!/usr/bin/env python
# coding: utf-8

# Brianna Birkel - August 2020

from obspy import read
from obspy.signal.filter import bandpass
import pandas as pd
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from obspy.signal.invsim import simulate_seismometer, corn_freq_2_paz

#Ignore warnings due to python 2 and 3 conflict
import warnings
warnings.filterwarnings("ignore")


path_dir = '/Users/bcbirkel/Documents/Research/LABasin/rawIRISdata/20200422/'
stream = read(path_dir + "C*")
stream_CE = read(path_dir + "CE*")
stream_CI = read(path_dir + "CI*")


plt.figure(figsize=(10,5))
# setup mercator map projection.
m = Basemap(projection='merc',llcrnrlon=-120, llcrnrlat=33, urcrnrlon=-117, urcrnrlat=35,epsg=4269)
#http://server.arcgisonline.com/arcgis/rest/services
#EPSG Number of America is 4269
m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 2000, verbose= True)

evlat = stream_CI[0].stats.sac.evla; evlon = stream_CI[0].stats.sac.evlo

# plot stations
for tr in stream_CE:
    stlat = tr.stats.sac.stla; stlon = tr.stats.sac.stlo 
    #m.drawgreatcircle(stlon,stlat,evlon,evlat,linewidth=1,color='b')
    xx,yy = m(stlon,stlat)
    m.scatter(xx, yy, marker = "^" ,s=100, c="g" , edgecolors = "k", alpha = 1, label = "CE")
    
for tr in stream_CI:
    stlat = tr.stats.sac.stla; stlon = tr.stats.sac.stlo 
    #m.drawgreatcircle(stlon,stlat,evlon,evlat,linewidth=1,color='b')
    xx,yy = m(stlon,stlat)
    m.scatter(xx, yy, marker = "^" ,s=100, c="y" , edgecolors = "k", alpha = 1, label = "CI")

#Plot the event
xx,yy = m(evlon,evlat)
m.scatter(xx, yy, marker = "*" ,s=200, c="r" , edgecolors = "k", alpha = 1)    
    
#m.drawparallels(np.arange(-90,90,20),labels=[1,1,0,1])
plt.title("Event-station map")
#plt.legend(loc='best', bbox_to_anchor=(0.5, 0., 0.5, 0.5))
plt.show()


plt.figure(figsize=(10,5))

# =============================================================================
# #need to find a good way to sort by distance so plots aren't clumped together 
#i = 0
# for tr in stream_CI:
#     dist(i) = tr.stats.sac.dist
#     tr_num
# 
# stream_plot = sorted(stream_CI, key=lambda dist: tr.stats.sac.dist)
# =============================================================================

stream_plot = stream_CI[::10] #plot every 10th station

# inst2hz = corn_freq_2_paz(2.0)
# sts2 = {'gain': 60077000.0,
#         'poles': [(-0.037004+0.037016j),
#                   (-0.037004-0.037016j),
#                   (-251.33+0j),
#                   (-131.04-467.29j),
#                   (-131.04+467.29j)],
#         'sensitivity': 2516778400.0,
#         'zeros': [0j, 0j]}
stream_plot.normalize(global_max=True)

for tr in stream_plot:
    df = tr.stats.sampling_rate
    #tr.data = simulate_seismometer(tr.data, df)
    tr.filter("bandpass", freqmin=1/20, freqmax=1/2)
    dist = tr.stats.sac.dist
    plt.plot(tr.times(),tr.data*300+dist,c="k",linewidth=0.5)
#    plt.scatter(tr.stats.sac.t3,dist*0.01,marker=" | ",color="r")
plt.ylabel("distance from source (km)")    
plt.xlabel("time (s)")
#plt.ylim(84,77)
#plt.xlim(650,800)
plt.show()


# =============================================================================
# plt.figure(figsize=(10,5))
# for tr in stream:
#     tr.normalize()
#     dist = tr.stats.sac.dist*0.01
#     x = tr.times()
#     y = tr.data+dist
#     plt.fill_between(x,y, dist, y > dist, color='r', alpha = 0.8)
#     plt.fill_between(x,y, dist, y < dist, color='b', alpha = 0.8)
# plt.ylabel("x100 km")    
# #plt.ylim(84,77)
# #plt.xlim(650,800)
# plt.show()
# 
# 
# =============================================================================
# In[ ]:




