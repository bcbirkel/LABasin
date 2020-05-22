#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 16:55:50 2019

@author: vidale
"""

import os
os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/bcbirkel/Documents/Research/LABasin/gitHub LAB/Quakes_IRIS')

from pro2_get_Z       import get_Z
from pro2_get_N       import get_N
from pro2_get_E       import get_E
from pro3_plot_sgrams import plot_sgrams
from re_pro2_get_Z    import re_get_Z

os.chdir('/Users/bcbirkel/Documents/Research/LABasin/PyCodeLAB/Mseed')

#%% define common parameters, but have none so far

#%% Treatment to compare events
#%% Cull seismic section for common stations

all_indices = [0,1,2,4,5,6,7,8,9,10,11,12,14,15]
get_indices = [17]
plot_indices = [14]
freqmin = 0.08
freqmax = 0.3
filt = 1
norm = 1
plot_scale_fac = 5
max_plot_dist = 80
min_plot_dist = 0
start_buff = 25
end_buff = 250
#	component = 'Z'

#for i in get_indices:
#	print('Getting event ' + str(i))
#	get_Z(select_event=i)
#	get_N(select_event=i)
#	get_E(select_event=i)

for i in plot_indices:
	print('Plotting event ' + str(i))
	plot_sgrams(select_event = i, component = 'Z', freqmin=freqmin, freqmax=freqmax, filt = filt,
			 norm = norm, plot_scale_fac = plot_scale_fac, start_buff = start_buff, end_buff = end_buff,
			 min_plot_dist = min_plot_dist, max_plot_dist = max_plot_dist)
#	plot_sgrams(select_event = i, component = 'N', freqmin=freqmin, freqmax=freqmax, filt = filt,
#			 norm = norm, plot_scale_fac = plot_scale_fac, start_buff = start_buff, end_buff = end_buff,
#			 min_plot_dist = min_plot_dist, max_plot_dist = max_plot_dist)
#	plot_sgrams(select_event = i, component = 'E', freqmin=freqmin, freqmax=freqmax, filt = filt,
#			 norm = norm, plot_scale_fac = plot_scale_fac, start_buff = start_buff, end_buff = end_buff,
#			 min_plot_dist = min_plot_dist, max_plot_dist = max_plot_dist)

#for i in get_indices:
#	print('Getting event ' + str(i))
#	re_get_Z(select_event=i)
#	get_N(select_event=i)
#	get_E(select_event=i)
