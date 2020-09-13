#!/usr/bin/env python
# this program plots two selected file types
# John Vidale 6/2019

import os

#from plot_stations    import p_stations
#from plot_comm_stations    import p_comm_stations
from plot_compare_two_comps  import p_compare2

ev_no = 1
traces1 = 'CVM_S4'
traces2 = 'CVM_H'
#traces1 = 'Data' #data v CVM-H
traces2 = 'Data' #data v CVM-S4
# component1 = 'E'
# component2 = component1
# component2 = 'N'
# Traces:  Data, CVM_H, CVM_S4
min_d = 0
max_d = 55
plot_scale_fac = 0.04
start_buff = -10
end_buff = 80
filt = True
freq_min = 0.05
freq_max = 0.5
taper = True
taper_frac = 0.02
norm_each = False
dist_norm = True
basin = True
basin_width = 15
motion = 'vel'
com_time  = True
rotated = True
CI_only = False

#%% define common parameters
events = ['lahabra_2014', 'beverlyhills_2001', 'chatsworth_2007', 'chinohills_2008', 'inglewood_2009']
event_no = ['15481673','9703873','14312160','14383980','10410337']
vmodel = ['cvmhn', 'cvmhy', 'cvms400-100', 'cvms426-223']

#%% Map
# p_stations(event_no = event_no[ev_no], traces1 = traces1, traces2 = traces2, min_dist = min_d, max_dist = max_d, CI_only = CI_only, basin_width = basin_width, basin = basin)

# p_comm_stations(event_no = event_no[ev_no], traces1 = traces1, traces2 = traces2, min_dist = min_d, max_dist = max_d, CI_only = CI_only, basin_width = basin_width, basin = basin)

#%% New plotting software
p_compare2(eventname = events[ev_no], event_no = event_no[ev_no], traces1 = traces1, traces2 = traces2, start_buff = start_buff,
           component1 = 'Z', component2 = 'Z', end_buff = end_buff, taper = taper,
           taper_frac = taper_frac, plot_scale_fac = plot_scale_fac, filt = filt, freq_min = freq_min,
           freq_max = freq_max, min_dist = min_d, max_dist = max_d, basin_width = basin_width,
           basin = basin, norm_each = norm_each, dist_norm = dist_norm, CI_only = CI_only)

p_compare2(eventname = events[ev_no], event_no = event_no[ev_no], traces1 = traces1, traces2 = traces2, start_buff = start_buff,
            component1 = 'N', component2 = 'N', end_buff = end_buff, taper = taper,
            taper_frac = taper_frac, plot_scale_fac = plot_scale_fac, filt = filt, freq_min = freq_min,
            freq_max = freq_max, min_dist = min_d, max_dist = max_d, basin_width = basin_width,
            basin = basin, norm_each = norm_each, dist_norm = dist_norm, CI_only = CI_only)

p_compare2(eventname = events[ev_no], event_no = event_no[ev_no], traces1 = traces1, traces2 = traces2, start_buff = start_buff,
            component1 = 'E', component2 = 'E', end_buff = end_buff, taper = taper,
            taper_frac = taper_frac, plot_scale_fac = plot_scale_fac, filt = filt, freq_min = freq_min,
            freq_max = freq_max, min_dist = min_d, max_dist = max_d, basin_width = basin_width,
            basin = basin, norm_each = norm_each, dist_norm = dist_norm, CI_only = CI_only)

p_compare2(eventname = events[ev_no], event_no = event_no[ev_no], traces1 = traces1, traces2 = traces2, start_buff = start_buff,
            component1 = 'R', component2 = 'R', end_buff = end_buff, taper = taper,
            taper_frac = taper_frac, plot_scale_fac = plot_scale_fac, filt = filt, freq_min = freq_min,
            freq_max = freq_max, min_dist = min_d, max_dist = max_d, basin_width = basin_width,
            basin = basin, norm_each = norm_each, dist_norm = dist_norm, CI_only = CI_only)

p_compare2(eventname = events[ev_no], event_no = event_no[ev_no], traces1 = traces1, traces2 = traces2, start_buff = start_buff,
            component1 = 'T', component2 = 'T', end_buff = end_buff, taper = taper,
            taper_frac = taper_frac, plot_scale_fac = plot_scale_fac, filt = filt, freq_min = freq_min,
            freq_max = freq_max, min_dist = min_d, max_dist = max_d, basin_width = basin_width,
            basin = basin, norm_each = norm_each, dist_norm = dist_norm, CI_only = CI_only)

