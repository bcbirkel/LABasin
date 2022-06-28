#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:04:13 2020

@author: bcbirkel
"""

from run_JS_syn import run_JS_syn
from compare_test1 import syn_compare
#from get_stations import get_stations

model = 'H'
ev_no = 0
min_d = 0
max_d = 55
plot_scale_fac = 0.05
plot_fac_js = .05
start_buff = 10
end_buff = 110
filt = True
freq_min = 0.1
freq_max = 0.2
taper = True
taper_frac = 0.05
show_H = True
show_S = False
show_H_simp = False
show_S_simp = False
show_data = True
norm_each = False
basin = False
basin_width = 15
motion = 'vel'
com_time  = True
rotated = False
CI_only = False
overlay = True
flip_z = False


if ev_no == 0: # la habra
    print ("setting vars for event 0")
    max_d = 55
#    plot_scale_fac = .05 
#    plot_fac_js = .03 
    flip_z = True
    CI_only = True

if ev_no == 1: # chino
    print ("setting vars for event 1")
    max_d = 60
#    plot_scale_fac = .05
#    plot_fac_js = .03
    CI_only = True
    
if ev_no == 2: # inglewood
    print ("setting vars for event 2")
    min_d = 5
    max_d = 20
#    plot_scale_fac = .05
#    plot_fac_js = .05
    
if ev_no == 3: #chatsworth
    print ("setting vars for event 3")
    max_d = 55
#    plot_scale_fac = .08
#    plot_fac_js = .03
    CI_only = False

if ev_no == 4: # bev hills
    print ("setting vars for event 4")
    max_d = 55
#    plot_scale_fac = .08
#    plot_fac_js = .05
    CI_only = False

event_no = ['15481673', '14383980', '10410337', '14312160', '9703873']
#event_no = ['9064568', '9093975', '9096972', '9140050', '9644101', '9703873', '9716853','9753489',
#	'9818433', '10216101', '10275733', '10370141', '10399889', '10403777', '10410337', '10530013',
#	'10541957', '10972299', '13692644', '14000376', '14116972', '14155260', '14239184', '14312160',
#	'14383980', '14494128', '14601172', '14738436', '15237281', '15481673']
vmodel = ['cvmhn', 'cvmhy', 'cvms400-100', 'cvms426-223']
if model == 'H':
    vmodel = vmodel[1]
elif model== 'S':
    vmodel = vmodel[3]
else:
    print("fix model param to be 'H' or 'S'")


run_JS_syn(event_num=ev_no,vmodel=vmodel)

#test getstations -- also is inside run_JS_syn
#update_sta_code_js, update_sta_seq_no_js, update_sta_code_rt = get_stations()
#print(type(update_sta_code_js))
#print(update_sta_code_rt)

#syn_compare(event_no = event_no[ev_no], start_buff = start_buff, end_buff = end_buff, taper = taper, taper_frac = taper_frac,
#			  plot_scale_fac = plot_scale_fac, plot_fac_js = plot_fac_js, filt = filt, freq_min = freq_min, freq_max = freq_max, min_dist = min_d, max_dist = max_d,
#			  vmodel = 'S', basin_width = basin_width, basin = basin, CI_only = CI_only, overlay=overlay, flip_z=flip_z)

#%%

# NEED TO MAP JS STATIONS OVER FOR PLOT_SYN_JS TO WORK
