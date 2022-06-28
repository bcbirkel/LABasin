#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 19:15:59 2020

@author: bcbirkel

%% Run script to plot Jon Stewart's synthetics 
"""
def run_JS_syn(event_num, vmodel):
    import os
    os.chdir('/Users/bcbirkel/Documents/Research/LABasin/stewart_synthetics/')
    from plot_syn_JS import plot_syn_JS
    from get_stations import get_stations
    
    ## Set up events/stations
    # events overlapping with Taborda - la habra, chino hills, inglewood, chatsworth, beverly hills
    
    event_no = ['pw_102', '1002', '1011', '1019', '1036']
    
    #stations overlapping with Taborda:
    get_stations()
    sta_code_js, sta_seq_js, sta_code_rt = get_stations()
    print(sta_seq_js)
    
    ## VARIABLES 
    eventfile = event_no[event_num]
    vmodel = vmodel
    print(vmodel)
    sta_num = 0
    
    if eventfile == 'pw_102':
        event_title = 'lahabra'
        event_lat = 33.9325
        event_lon = -117.9158
    elif eventfile == '1002':
        event_title = 'chinohills'
        event_lat = 33.9465
        event_lon = -117.7667
    elif eventfile == '1011':
        event_title = 'inglewood'
        event_lat = 33.9377
        event_lon = -118.3357
    elif eventfile == '1019':
        event_title = 'chatsworth'
        event_lat = 34.2983
        event_lon = -118.6255
    elif eventfile == '1036':
        event_title = 'beverlyhills'
        event_lat = 34.0541
        event_lon = -118.3929
    else: 
        print('unknown event file')
    
    time, stN, stE, stZ = plot_syn_JS(eventfile = eventfile, event_title = event_title, vmodel = vmodel, stations = sta_seq_js, sta_names = sta_code_js, event_lat = event_lat, event_lon = event_lon)
    
    return sta_code_rt, sta_num, time, stN, stE, stZ

