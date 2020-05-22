#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 15:48:32 2020

@author: bcbirkel
"""

#!/usr/bin/env python
# this program plots Ricardo's synthetic and data files
# John Vidale 6/2019

import os

#os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/bcbirkel/Documents/Research/LABasin/github LAB/Quakes')

from plot_data_vs_syn         import p_data_vs_syn
from plot_3Cdata              import p_3Cdata
from plot_3Csyn               import p_3Csyn
from plot_3Csyn               import p_3Csyn
from rot_data_RT              import rot_data_RT
from rot_syn_RT               import rot_syn_RT

ev_no = 0
min_d = 0
max_d = 60
plot_scale_fac = 0.05
start_buff = 10
end_buff = 110
filt = True
freq_min = 0.1
freq_max = 0.5
taper = True
taper_frac = 0.05
show_H = True
show_S = False
show_H_simp = False
show_S_simp = False
show_data = True
norm_each = False
basin = True
basin_width = 15
motion = 'vel'
com_time  = True
rotated = False
CI_only = True

#%% define common parameters
# 29 events

event_no = ['15481673', '14383980', '10410337', '14312160', '970387']
#event_no = ['9064568', '9093975', '9096972', '9140050', '9644101', '9703873', '9716853','9753489',
#	'9818433', '10216101', '10275733', '10370141', '10399889', '10403777', '10410337', '10530013',
#	'10541957', '10972299', '13692644', '14000376', '14116972', '14155260', '14239184', '14312160',
#	'14383980', '14494128', '14601172', '14738436', '15237281', '15481673']
vmodel = ['cvmhn', 'cvmhy', 'cvms400-100', 'cvms426-223']

#%% Rotate all data and look at it
#for evno in event_no:
#	rot_data_RT(event_no = evno, start_buff = start_buff, end_buff = end_buff, taper = taper, taper_frac = taper_frac,
#				  plot_scale_fac = plot_scale_fac, filt = filt, freq_min = freq_min, freq_max = freq_max,
#				  min_dist = min_d, max_dist = max_d, norm_each = norm_each, motion = motion, basin_width = basin_width, basin = basin)

# rotate one
#rot_data_RT(event_no = '9703873', start_buff = start_buff, end_buff = end_buff, taper = taper, taper_frac = taper_frac,
#			  plot_scale_fac = plot_scale_fac, filt = filt, freq_min = freq_min, freq_max = freq_max,
#			  min_dist = min_d, max_dist = max_d, norm_each = norm_each, motion = motion, basin_width = basin_width, basin = basin)

#%% Rotate synthetics and look at it
#for evno in event_no:
#	for vmod in vmodel:
#		rot_syn_RT(event_no = evno, start_buff = start_buff, end_buff = end_buff, taper = taper, taper_frac = taper_frac,
#					plot_scale_fac = plot_scale_fac, filt = filt, freq_min = freq_min, freq_max = freq_max, min_dist = min_d, max_dist = max_d,
#					vmodel = vmod,  basin_width = basin_width, basin = basin)

#%% SYN
# Harvard model
p_3Csyn(event_no = event_no[ev_no], start_buff = start_buff, end_buff = end_buff, taper = taper, taper_frac = taper_frac,
			  plot_scale_fac = plot_scale_fac, filt = filt, freq_min = freq_min, freq_max = freq_max, min_dist = min_d, max_dist = max_d,
			  vmodel = 'H', basin_width = basin_width, basin = basin, CI_only = CI_only)

# SCEC model
p_3Csyn(event_no = event_no[ev_no], start_buff = start_buff, end_buff = end_buff, taper = taper, taper_frac = taper_frac,
			  plot_scale_fac = plot_scale_fac, filt = filt, freq_min = freq_min, freq_max = freq_max, min_dist = min_d, max_dist = max_d,
			  vmodel = 'S', basin_width = basin_width, basin = basin, CI_only = CI_only)

"""
#%% DATA
p_3Cdata(event_no = event_no[ev_no], start_buff = start_buff, end_buff = end_buff, taper = taper, taper_frac = taper_frac,
			plot_scale_fac = plot_scale_fac, filt = filt, freq_min = freq_min, freq_max = freq_max,
			min_dist = min_d, max_dist = max_d, norm_each = norm_each, motion = motion,
			basin_width = basin_width, basin = basin, rotated = rotated, CI_only = CI_only)

# %%SYN vs DATA
p_data_vs_syn(event_no = event_no[ev_no], start_buff = start_buff, end_buff = end_buff, taper = taper, taper_frac = taper_frac,
			  plot_scale_fac = plot_scale_fac, filt = filt, freq_min = freq_min, freq_max = freq_max, min_dist = min_d, max_dist = max_d,
			  show_data = show_data, show_H = show_H, show_S = show_S, show_H_simp = show_H_simp, show_S_simp = show_S_simp,
			  compo = 'E', norm_each = norm_each, motion = motion, basin_width = basin_width, basin = basin,
			  com_time = com_time, rotated = rotated, CI_only = CI_only)

p_data_vs_syn(event_no = event_no[ev_no], start_buff = start_buff, end_buff = end_buff, taper = taper, taper_frac = taper_frac,
			  plot_scale_fac = plot_scale_fac, filt = filt, freq_min = freq_min, freq_max = freq_max, min_dist = min_d, max_dist = max_d,
			  show_data = show_data, show_H = show_H, show_S = show_S, show_H_simp = show_H_simp, show_S_simp = show_S_simp,
			  compo = 'N', norm_each = norm_each, motion = motion, basin_width = basin_width, basin = basin,
			  com_time = com_time, rotated = rotated, CI_only = CI_only)

p_data_vs_syn(event_no = event_no[ev_no], start_buff = start_buff, end_buff = end_buff, taper = taper, taper_frac = taper_frac,
			  plot_scale_fac = plot_scale_fac, filt = filt, freq_min = freq_min, freq_max = freq_max, min_dist = min_d, max_dist = max_d,
			  show_data = show_data, show_H = show_H, show_S = show_S, show_H_simp = show_H_simp, show_S_simp = show_S_simp,
			  compo = 'Z', norm_each = norm_each, motion = motion, basin_width = basin_width, basin = basin,
			  com_time = com_time, rotated = rotated, CI_only = CI_only)
"""

# returns to a convenient directory when done
os.chdir('/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo')
