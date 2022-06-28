#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 11:12:26 2020

@author: bcbirkel
"""

def run_RT_compare():
    
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    import sys
    from obspy import read_events, UTCDateTime, read, Trace, Stream
    from obspy.geodetics import gps2dist_azimuth
    
    	#%% Load data and synthetic waveforms
	st_dat = Stream()
	st_synE = Stream()
	st_synN = Stream()
	st_synZ = Stream()
	if vmodel == 'H':
		fname_dat     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_data/vel/ve_' + event_no + '.mseed'
		fname_synE     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/vel/cvms400-100/ve_'  + event_no + '_cvms400-100.mseed'
		fname_synN     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/vel/cvms400-100/vn_'  + event_no + '_cvms400-100.mseed'
		fname_synZ     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/vel/cvms400-100/vz_'  + event_no + '_cvms400-100.mseed'
	elif vmodel == 'S':
		fname_dat     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_data/vel/ve_' + event_no + '.mseed'
		fname_synE     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/vel/cvms426-223/ve_'  + event_no + '_cvms426-223.mseed'
		fname_synN     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/vel/cvms426-223/vn_'  + event_no + '_cvms426-223.mseed'
		fname_synZ     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/vel/cvms426-223/vz_'  + event_no + '_cvms426-223.mseed'
	st_dat=read(fname_dat)
	st_synE=read(fname_synE)
	st_synN=read(fname_synN)
	st_synZ=read(fname_synZ)
	print('In sgram file ' + str(st_synZ[609].data[0]) + '  ' + st_synZ[609].stats.station + '  ' + str(len(st_synZ)))
	print('In arrays st_name[609] ' + st_name[609] + ' st_name[610] ' + st_name[610])
	print('In arrays st_name[0] ' + st_name[0] + ' st_name[1] ' + st_name[1])

	print('1st data trace has : ' + str(len(st_synE[0].data)) + ' time pts ')
	print('synE has ' + str(len(st_synE)) + ' traces')
	print('synN has ' + str(len(st_synN)) + ' traces')
	print('synZ has ' + str(len(st_synZ)) + ' traces')
    
    plt.plot(data, figure='7')