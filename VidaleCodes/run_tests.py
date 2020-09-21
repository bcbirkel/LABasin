#!/usr/bin/env python
# this program plots two selected file types
# John Vidale 6/2019

import os

from plot_stations    import p_stations
from plot_comm_stations    import p_comm_stations
from plot_compare_two_comps  import p_compare2

ev_no = 24  # LaHabra 29, Inglewood 14, ChinoHills 24, Chatsworth 23, BeverlyHills 5
traces1 = 'Data'
traces2 = 'CVM_S4'
# component1 = 'E'
# component2 = component1
# component2 = 'N'
# Traces:  Data, CVM-H, CVM-S4
min_d = 10
max_d = 70
plot_scale_fac = 0.04
start_buff = -10
end_buff = 80
filt = True
freq_min = 0.05
freq_max = 0.333
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
event_no = ['9064568', '9093975', '9096972', '9140050', '9644101', '9703873', '9716853','9753489',
    '9818433', '10216101', '10275733', '10370141', '10399889', '10403777', '10410337', '10530013',
    '10541957', '10972299', '13692644', '14000376', '14116972', '14155260', '14239184', '14312160',
    '14383980', '14494128', '14601172', '14738436', '15237281', '15481673', '38824959', '39020663',
    '39322287', '39322767']

#%% Map
# p_stations(event_no = event_no[ev_no], traces1 = traces1, traces2 = traces2, min_dist = min_d, max_dist = max_d, CI_only = CI_only, basin_width = basin_width, basin = basin)

# p_comm_stations(event_no = event_no[ev_no], traces1 = traces1,
#                 traces2 = traces2, min_dist = min_d, max_dist = max_d,
#                 CI_only = CI_only, basin_width = basin_width, basin = basin)

#%% New plotting software
fig_inc = 0
p_compare2(event_no = event_no[ev_no], traces1 = traces1, traces2 = traces2, start_buff = start_buff,
            component1 = 'Z', component2 = 'Z', end_buff = end_buff, taper = taper,
            taper_frac = taper_frac, plot_scale_fac = plot_scale_fac, filt = filt, freq_min = freq_min,
            freq_max = freq_max, min_dist = min_d, max_dist = max_d, basin_width = basin_width,
            basin = basin, norm_each = norm_each, dist_norm = dist_norm, CI_only = CI_only, fig_inc = fig_inc)

# p_compare2(event_no = event_no[ev_no], traces1 = traces1, traces2 = traces2, start_buff = start_buff,
#             component1 = 'N', component2 = 'E', end_buff = end_buff, taper = taper,
#             taper_frac = taper_frac, plot_scale_fac = plot_scale_fac, filt = filt, freq_min = freq_min,
#             freq_max = freq_max, min_dist = min_d, max_dist = max_d, basin_width = basin_width,
#             basin = basin, norm_each = norm_each, dist_norm = dist_norm, CI_only = CI_only, fig_inc = fig_inc)

# p_compare2(event_no = event_no[ev_no], traces1 = traces1, traces2 = traces2, start_buff = start_buff,
#             component1 = 'E', component2 = 'N', end_buff = end_buff, taper = taper,
#             taper_frac = taper_frac, plot_scale_fac = plot_scale_fac, filt = filt, freq_min = freq_min,
#             freq_max = freq_max, min_dist = min_d, max_dist = max_d, basin_width = basin_width,
#             basin = basin, norm_each = norm_each, dist_norm = dist_norm, CI_only = CI_only, fig_inc = fig_inc)

p_compare2(event_no = event_no[ev_no], traces1 = traces1, traces2 = traces2, start_buff = start_buff,
            component1 = 'R', component2 = 'R', end_buff = end_buff, taper = taper,
            taper_frac = taper_frac, plot_scale_fac = plot_scale_fac, filt = filt, freq_min = freq_min,
            freq_max = freq_max, min_dist = min_d, max_dist = max_d, basin_width = basin_width,
            basin = basin, norm_each = norm_each, dist_norm = dist_norm, CI_only = CI_only, fig_inc = fig_inc)

p_compare2(event_no = event_no[ev_no], traces1 = traces1, traces2 = traces2, start_buff = start_buff,
            component1 = 'T', component2 = 'T', end_buff = end_buff, taper = taper,
            taper_frac = taper_frac, plot_scale_fac = plot_scale_fac, filt = filt, freq_min = freq_min,
            freq_max = freq_max, min_dist = min_d, max_dist = max_d, basin_width = basin_width,
            basin = basin, norm_each = norm_each, dist_norm = dist_norm, CI_only = CI_only, fig_inc = fig_inc)

