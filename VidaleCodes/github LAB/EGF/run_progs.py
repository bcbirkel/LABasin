#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 16:55:50 2019

@author: vidale
"""

import os
os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/vidale/Documents/GitHub/LAB/EGF')

from plot_EGFs        import plot_EGFs

os.chdir('/Users/vidale/Documents/PyCode/LAB/Mseed')

#%% define common parameters, but have none so far

#%% Treatment to compare events
#%% Cull seismic section for common stations

all_indices = [0,1,2,4,5,6,7,8,9,10,11,12,14,15]
get_indices = [17]
plot_indices = [6]
freqmin = 0.1
freqmax = 0.5
filt = 1
norm = 1
plot_scale_fac = 3
max_plot_dist = 60
min_plot_dist = 0
start_buff = 5
end_buff = 100
#	component = 'Z'

#for i in get_indices:
#	print('Getting event ' + str(i))
#	get_Z(select_event=i)
#	get_N(select_event=i)
#	get_E(select_event=i)

station = 'PDR'
#freqmin = 0.1
#freqmax = 1.0
plot_EGFs(comp='ZZ',sta=station,freqmin=freqmin,freqmax=freqmax)
plot_EGFs(comp='TT',sta=station,freqmin=freqmin,freqmax=freqmax)

#for i in get_indices:
#	print('Getting event ' + str(i))
#	re_get_Z(select_event=i)
#	get_N(select_event=i)
#	get_E(select_event=i)
