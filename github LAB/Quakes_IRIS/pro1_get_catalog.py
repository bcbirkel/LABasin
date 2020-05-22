#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 17:13:10 2019

@author: vidale
"""
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt
import os

client = Client('SCEDC')

os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB')

min_mag = 3.5
min_lat = 33.75
max_lat = 34.2
min_lon = -118.5
max_lon = -117.75
t1 = UTCDateTime("1996-01-01T00:00:00")
t2 = UTCDateTime("2020-01-01T00:00:00")

fname_cat = 'LAB.quakeml3'

catalog = client.get_events(starttime = t1, endtime = t2, minmagnitude = min_mag,
					minlatitude  = min_lat, maxlatitude  = max_lat,
					minlongitude = min_lon, maxlongitude = max_lon)
#ev_lon   = catalog[0].origins[0].longitude
#ev_lat   = catalog[0].origins[0].latitude
#ev_depth = catalog[0].origins[0].depth
#ev_t      = catalog[0].origins[0].time
#print('event:',catalog)
print(catalog.__str__(print_all=True))
plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = 12, 8
#catalog.plot(projection = 'local', resolution = 'h')
catalog.write(fname_cat, format='QUAKEML')
