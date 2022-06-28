#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 15:48:14 2020

@author: bcbirkel
"""

def syn_compare(event_no, start_buff, end_buff, taper, taper_frac,
            plot_scale_fac, plot_fac_js, filt, freq_min, freq_max, min_dist, max_dist,
            vmodel, basin_width, basin, CI_only, overlay, flip_z):

    from obspy import UTCDateTime
    from obspy import Stream
    from obspy import read
    from obspy import Trace
    from obspy.core.trace import Stats
    from obspy.geodetics import gps2dist_azimuth
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    import time

    import sys # don't show any warnings
    import warnings

    if not sys.warnoptions:
        warnings.simplefilter("ignore")

    show_data = 0
    start_time_wc = time.time()
    bad_sta = True        
    
    if event_no == '15481673':
        eventname = 'pw_102'
    if event_no == '14383980':
        eventname = '1002'
    if event_no == '10410337':
        eventname = '1011'
        bad_sta = True
    if event_no == '14312160':
        eventname = '1019'
        bad_sta = False
    if event_no == '9703873':
        eventname = '1036'
        

    #%% find event details for origin time, lat, lon
    ev_file = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/event_list.txt'
    file_ev = open(ev_file, 'r')
    for line in file_ev:           # pull numbers off the rest of the lines
        split_line = line.split()
        event = split_line[0]
        if event == event_no:
            ev_lat       = float(split_line[3])
            ev_lon       = float(split_line[2])
            t1           = UTCDateTime(split_line[5])
            date_label  = split_line[5][0:10]
            year1        = split_line[5][0:4]
    print(event_no + str(t1) + ' ' + date_label + ' ' + year1 + '  ' + str(ev_lat) + ' ' + str(ev_lon))

    #%% 
    badt_file = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_data/bad_trace.txt'
    file_badt = open(badt_file, 'r')
    badt_lines = file_badt.readlines()
    badt_event   = []
    badt_station = []
    badt_compo   = []

    for line in badt_lines:           # pull numbers off all the lines
        split_line = line.split()
        badt_event.append(  split_line[0])
        badt_station.append(split_line[1])
        badt_compo.append(  split_line[2])
    print(str(len(badt_event)) + ' bad traces in list')

    #%% Open station location file
    sta_file = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/ricardo_stations.txt'
    file_st = open(sta_file, 'r')
    line = file_st.readline()      # read first line to skip header information
    lines = file_st.readlines()
    print(str(len(lines)) + ' stations read from ' + sta_file)

    # Load station coords into arrays, many more stations than used
    station_index = range(len(lines))
    st_num   = []
    st_netw  = []
    st_name  = []
    st_dist  = []
    st_az    = []
    st_baz   = []
    st_lat   = []
    st_lon   = []
    for line in lines:
        split_line = line.split()
        st_num.append( split_line[0])
        st_netw.append(split_line[2])
        st_name.append(split_line[3])
        st_lat.append( split_line[4])
        st_lon.append( split_line[5])
        distance = gps2dist_azimuth( ev_lat, ev_lon, float(split_line[4]), float(split_line[5])) # Get traveltime and azimuth
        st_dist.append(distance[0]/1000.) # azimuth
        st_az.append(distance[1]) # azimuth
        st_baz.append(distance[2]) # back-azimuth
    print('number of stations in list is ' + str(len(st_num)) + ' or ' + str(station_index))

    #%% Load data and synthetic waveforms
    st_dat = Stream()
    st_synE = Stream()
    st_synN = Stream()
    st_synZ = Stream()
    js_synE = Stream()
    js_synN = Stream()
    js_synZ = Stream()
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
        fname_jsT      = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/compare_syn/time-'+ eventname +'.mseed'
        fname_jsE      = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/compare_syn/E_bp-'+ eventname +'.mseed'
        fname_jsN      = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/compare_syn/N_bp-'+ eventname +'.mseed'
        fname_jsZ      = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/compare_syn/Z_bp-'+ eventname +'.mseed'
    st_dat=read(fname_dat)
    st_synE=read(fname_synE)
    st_synN=read(fname_synN)
    st_synZ=read(fname_synZ)
    st_jsT=read(fname_jsT)
    st_jsE=read(fname_jsE)
    st_jsN=read(fname_jsN)
    st_jsZ=read(fname_jsZ)
    print('In sgram file ' + str(st_synZ[609].data[0]) + '  ' + st_synZ[609].stats.station + '  ' + str(len(st_synZ)))
    print('In arrays st_name[609] ' + st_name[609] + ' st_name[610] ' + st_name[610])
    print('In arrays st_name[0] ' + st_name[0] + ' st_name[1] ' + st_name[1])

    print('1st data trace has : ' + str(len(st_synE[0].data)) + ' time pts ')
    print('synE has ' + str(len(st_synE)) + ' traces')
    print('synN has ' + str(len(st_synN)) + ' traces')
    print('synZ has ' + str(len(st_synZ)) + ' traces')
    print('jsZ has ' + str(len(st_jsZ)) + ' traces')
    print(st_jsZ[0])
    #%%
    # select data by distance (and azimuth?), and cull synthetics to match data
    st_dat_select = Stream()
    st_synE_select = Stream()
    st_synN_select = Stream()
    st_synZ_select = Stream()
    js_synE_select = Stream()
    js_synN_select = Stream()
    js_synZ_select = Stream()
    for tr in st_dat: # examine traces one by one
        if tr.stats.network == 'CI' or CI_only == False:
            for ii in range(len(st_name)):  # find matching entry in station roster
                if (tr.stats.network == st_netw[ii] and tr.stats.station == st_name[ii]): # find station in inventory
                    if (st_dist[ii] < max_dist) and (st_dist[ii] > min_dist): # exclude stations too close or too far
                        # basin is roughly within 15 km of these 3 stations
                        distance = gps2dist_azimuth( 33.99053, -118.36171, float(st_lat[ii]), float(st_lon[ii])) # station BHP
                        dist1 = distance[0]/1000  # convert m to km
                        distance = gps2dist_azimuth( 33.88110, -118.17568, float(st_lat[ii]), float(st_lon[ii])) # station LTP
                        dist2 = distance[0]/1000
                        distance = gps2dist_azimuth( 33.80776, -117.98116, float(st_lat[ii]), float(st_lon[ii])) # station BRE
                        dist3 = distance[0]/1000
    #                    print(tr.stats.station + ' ' + str(dist1) + ' ' + str(dist2) + ' ' + str(dist3))
                        if basin == False or (dist1 < basin_width) or (dist2 < basin_width) or (dist3 < basin_width):  # keep stations only within X km of basin axis
    #                        print('selected: ' + tr.stats.station)
    #                        tr.stats.distance = st_dist[ii] # add distance to trace metadata
    #                        st_datE_select += tr
                            tr.stats.distance = st_dist[ii] # add distance to trace metadata
                            st_synE[ii].stats.distance = st_dist[ii]
                            st_synN[ii].stats.distance = st_dist[ii]
                            st_synZ[ii].stats.distance = st_dist[ii]
                            st_dat_select += tr
                            st_synE_select += st_synE[ii]
                            st_synN_select += st_synN[ii]
                            st_synZ_select += st_synZ[ii]
#                            for jstr in st_jsE:
#                                print(jstr.stats.get('station'))
#                                if jstr.stats.get('station') == tr.stats.station:
#                                    js_synE_select += jstr
#                            for jstr in st_jsN:
#                                if jstr.stats.get('station') == tr.stats.station:
#                                    js_synN_select += jstr
#                            for jstr in st_jsZ:
#                                if jstr.stats.get('station') == tr.stats.station:
#                                    js_synZ_select += jstr
                            print(tr.stats.station + ' counter ' + str(ii) + ' ' + 'num ' + str(int(st_num[ii])))
                        else:
                            print('Not in basin: '  + tr.stats.station)
                    else:
                        print('Too far: '  + tr.stats.station)       
    print('now data has ' + str(len(st_dat_select)) + ' traces')
    print('now synE has '  + str(len(st_synE_select)) + ' traces')
    print('now synN has '  + str(len(st_synN_select)) + ' traces')
    print('now synZ has '  + str(len(st_synZ_select)) + ' traces')

    #%%  reject data and sythetics on bad trace list, either individual components or A for all components
    ## when culling for like stations, may be missing stations in which ricardo's sta are denoted by letters.
    
    if bad_sta == True:
        st_dat_good = Stream()
        for tr in st_dat_select: # examine traces one by one
            do_write = 1
            for ii in range(len(badt_event)):
                if event_no == badt_event[ii] and tr.stats.station == badt_station[ii]: # find station in inventory
                    if badt_compo[ii] == 'A' or badt_compo[ii] == 'E' or badt_compo[ii] == 'N' or badt_compo[ii] == 'Z':
                        do_write = 0
            if do_write == 1:
                st_dat_good += tr
        print('After rejecting labeled bad traces ones, dat has '       + str(len(st_dat_good))       + ' traces')
    
        st_synE_good = Stream()
        js_synE_good = Stream()
        for tr in st_synE_select: # examine traces one by one
            do_write = 1
            for ii in range(len(badt_event)):
                if event_no == badt_event[ii] and tr.stats.station == badt_station[ii]: # find station in inventory
                    if badt_compo[ii] == 'E' or badt_compo[ii] == 'A':
                        do_write = 0
            if do_write == 1:
                for jstr in st_jsE:
                    if jstr.stats.get('station') == tr.stats.station:
                        js_synE_good += jstr
                        st_synE_good += tr
        print('After rejecting labeled bad traces ones, synE has '       + str(len(st_synE_good))       + ' traces')
        print('After rejecting labeled bad traces ones, js_synE has '       + str(len(js_synE_good))       + ' traces')
    
        st_synN_good = Stream()
        js_synN_good = Stream()
        for tr in st_synN_select: # examine traces one by one
            do_write = 1
            for ii in range(len(badt_event)):
                if event_no == badt_event[ii] and tr.stats.station == badt_station[ii]: # find station in inventory
                    if badt_compo[ii] == 'N' or badt_compo[ii] == 'A':
                        do_write = 0
            if do_write == 1:
                for jstr in st_jsN:
                    #print(tr.stats.station)
                    #print(jstr.stats.get('station'))
                    if jstr.stats.get('station') == tr.stats.station:
                        js_synN_good += jstr
                        st_synN_good += tr
        print('After rejecting labeled bad traces ones, synN has '       + str(len(st_synN_good))       + ' traces')
        print('After rejecting labeled bad traces ones, js_synN has '       + str(len(js_synN_good))       + ' traces')
        
        st_synZ_good = Stream()
        js_synZ_good = Stream()
        for tr in st_synZ_select: # examine traces one by one
            do_write = 1
            for ii in range(len(badt_event)):
                if event_no == badt_event[ii] and tr.stats.station == badt_station[ii]: # find station in inventory
                    if badt_compo[ii] == 'Z' or badt_compo[ii] == 'A':
                        do_write = 0
            if do_write == 1:
                for jstr in st_jsZ:
                    if jstr.stats.get('station') == tr.stats.station:
                        js_synZ_good += jstr
                        st_synZ_good += tr
        print('After rejecting labeled bad traces ones, synZ has '       + str(len(st_synZ_good))       + ' traces')
        print('After rejecting labeled bad traces ones, js_synZ has '       + str(len(js_synZ_good))       + ' traces')
      
    else:
        st_dat_good = st_dat_select
        st_synE_good = st_synE_select
        st_synN_good = st_synN_select
        st_synZ_good = st_synZ_select
        js_synE_good = js_synE_select
        js_synN_good = js_synN_select
        js_synZ_good = js_synZ_select
        
        print('now js_synE has ' + str(len(js_synE_good)) + ' traces')
        print('now synE has '  + str(len(st_synE_good)) + ' traces')

    
    #%%  detrend, taper, filter
    if taper:
        st_dat_good.detrend( type='simple')
        st_synE_good.detrend(type='simple')
        st_synN_good.detrend(type='simple')
        st_synZ_good.detrend(type='simple')
#        js_synE_good.detrend(type='simple')
#        js_synN_good.detrend(type='simple')
#        js_synZ_good.detrend(type='simple')
    if filt:
        st_dat_good.filter( 'bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
        st_synE_good.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
        st_synN_good.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
        st_synZ_good.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
        js_synE_good.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
        js_synN_good.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
        js_synZ_good.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
    if taper:
        st_dat_good.detrend( type='simple')
        st_synE_good.detrend(type='simple')
        st_synN_good.detrend(type='simple')
        st_synZ_good.detrend(type='simple')
        js_synE_good.detrend(type='simple')
        js_synN_good.detrend(type='simple')
        js_synZ_good.detrend(type='simple')

    # JS synthetics have 32000 data points; ricardos have 1001 
    ##  INSTEAD OF DECIMATING DATA, INTERPOLATE TIME.  y = y1 + ((y2 - y1) / (x2 - x1)) * (x - x1)
    #jsT_good = st_synZ_good[0].interpolate(method="linear", npts=32000, sampling_rate=640)
#    for tr in st_jsT:
#        for ii in tr:
#            jsT_good[ii] = st_jsT[ii]*0.05
#    
#    print(jsT_good)
        
#    js_synE_good.decimate(16)
#    js_synN_good.decimate(16)
#    js_synZ_good.decimate(16)
#    js_synE_good.decimate(2)
#    js_synN_good.decimate(2)
#    js_synZ_good.decimate(2)
    #%%
    # plot traces
    if vmodel == 'S':
        fig_index = 8
    elif vmodel == 'H':
        fig_index = 9
    plt.close(fig_index)
    plt.figure(fig_index,figsize=(10,8))
    plt.xlim(-start_buff,end_buff)
    plt.ylim(min_dist,max_dist)

    # find max
    maxE = 0
    js_maxE = 0
    for tr in st_synE_good:
        tr_max = max(abs(tr.data))*tr.stats.distance
        if tr_max > maxE:
            maxE = tr_max
         
    count = 0
    for tr in js_synE_good:
        tr_max = max(abs(tr.data))*st_synE_good[count].stats.distance
        if tr_max > maxE:
            js_maxE = tr_max
        count += 1
    
    
    maxN = 0
    js_maxN = 0
    count = 0    
    for tr in st_synN_good:
        tr_max = max(abs(tr.data))*tr.stats.distance
        if tr_max > maxN:
            maxN = tr_max
            
    for tr in js_synN_good:
        tr_max = max(abs(tr.data))*st_synE_good[count].stats.distance
        if tr_max > maxE:
            js_maxN = tr_max
        count += 1

    maxZ = 0
    js_maxZ = 0
    count = 0
    for tr in st_synZ_good:
        tr_max = max(abs(tr.data))*tr.stats.distance
        if tr_max > maxZ:
            maxZ = tr_max
            
    for tr in js_synZ_good:
        tr_max = max(abs(tr.data))*st_synE_good[count].stats.distance
        if tr_max > maxE:
            js_maxZ = tr_max
        count += 1
            
    print('Max E, N, and Z synthetic are ' + str(maxE) + '  ' + str(maxN) + '  ' + str(maxZ))
    print('Max E, N, and Z synthetic for JS are ' + str(js_maxE) + '  ' + str(js_maxN) + '  ' + str(js_maxZ))

    max_all = max(maxZ, maxN, maxE)
    plot_fac = plot_scale_fac * (max_dist - min_dist) / max_all
    
    js_max_all = max(js_maxZ, js_maxN, js_maxE)
    js_plot_fac = plot_fac_js * (max_dist - min_dist) / js_max_all

    if show_data:
        for tr in st_dat_good:
            dist_offset = tr.stats.distance # km
            ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)
    #        if red_plot == 1:
    #            shift = red_time + (dist_offset - red_dist) * red_slow
    #            ttt = ttt - shift
    #        plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
    #            - tr.data.min()) + dist_offset, color = 'green')
            plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
                - tr.data.min()) + dist_offset, color = 'black')
            print(str(tr.stats.distance) + ' distance ' + tr.stats.station + ' station')    #plt.title(fname1)

    print(st_jsT)
    #   print(jsT_good)
    
#    time = Stats()
#    for tr in st_synE_good:
#        time.distance += tr.stats.distance
#    js_synE_good.stats = Stats()
#    js_synE_good.stats.distance = st_synE_good.stats.distance
    
    #print labels whether or not data is shown
    for tr in st_dat_good:
        dist_offset = tr.stats.distance # km

    if overlay == True:
        #Plot East comp
        plt.close(1)
        plt.figure(1,figsize=(10,8))
        plt.xlim(-start_buff,end_buff)
        plt.ylim(min_dist,max_dist)
        count = 0
        
        for tr in st_synE_good:
            dist_offset = tr.stats.distance
            ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)#    These lines used to cause a crash in Spyder
            plt.plot(ttt, (tr.data * plot_fac * tr.stats.distance) + dist_offset, color = 'red')
            plt.text(s = tr.stats.network + ' ' + tr.stats.station ,x = end_buff*0.95,y = dist_offset
            + max_dist*0.015, color = 'black')  #label traces

        for jstr in js_synE_good:   
            dist_offset = st_synE_good[count].stats.distance
            ttt = np.arange(0, 400, 0.05)
            #ttt = np.arange(len(jstr)) * (0.05) + (st_synE_good[count].stats.starttime - t1 - 50)#    These lines used to cause a crash in Spyder
            plt.plot(ttt, (jstr.data[0:8000] * js_plot_fac * st_synE_good[count].stats.distance) + dist_offset, color = 'blue')
            plt.text(s = st_synE_good[count].stats.network + ' ' + st_synE_good[count].stats.station ,x = end_buff*0.95,y = dist_offset
            + max_dist*0.015, color = 'black')  #label traces
    
            print(count)
            print(dist_offset)
            print(jstr.data[0:8000])
            print(jstr.stats.get('station'))
            print(st_synE_good[count].stats.station)
            print(st_synE_good[count].stats.distance)
            
            count += 1
            
        print(st_synE_good[0].stats.delta)
        print(st_synE_good[1].stats.delta)
    
        plt.xlabel('Time (s)')
        plt.ylabel('Epicentral distance from event (km)')
        plt.title('Jon Stewart (blue) vs Ricardo Taborda (red) Syn ' + date_label + ' ' + event_no + ' comp: East ' + vmodel)
        plt.show()
        
        #plot North comp
        plt.close(2)
        plt.figure(2,figsize=(10,8))
        plt.xlim(-start_buff,end_buff)
        plt.ylim(min_dist,max_dist)
        count = 0
        
        for tr in st_synN_good:
            dist_offset = tr.stats.distance
            ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)#    These lines used to cause a crash in Spyder
            plt.plot(ttt, (tr.data * plot_fac * tr.stats.distance) + dist_offset, color = 'red')
            plt.text(s = tr.stats.network + ' ' + tr.stats.station ,x = end_buff*0.95,y = dist_offset
            + max_dist*0.015, color = 'black')  #label traces

        for jstr in js_synN_good:   
            dist_offset = st_synE_good[count].stats.distance
            ttt = np.arange(0, 400, 0.05)
            #ttt = np.arange(len(jstr)) * (0.05) + (st_synE_good[count].stats.starttime - t1 - 50)#    These lines used to cause a crash in Spyder
            plt.plot(ttt, (jstr.data[0:8000] * js_plot_fac * st_synE_good[count].stats.distance) + dist_offset, color = 'blue')
            plt.text(s = st_synE_good[count].stats.network + ' ' + st_synE_good[count].stats.station ,x = end_buff*0.95,y = dist_offset
            + max_dist*0.015, color = 'black')  #label traces
    
            print(count)
            print(dist_offset)
            print(jstr.data[0:8000])
            print(jstr.stats.get('station'))
            print(st_synE_good[count].stats.station)
            print(st_synE_good[count].stats.distance)
            
            count += 1
            
        print(st_synE_good[0].stats.delta)
        print(st_synE_good[1].stats.delta)
    
        plt.xlabel('Time (s)')
        plt.ylabel('Epicentral distance from event (km)')
        plt.title('Jon Stewart (blue) vs Ricardo Taborda (red) Syn ' + date_label + ' ' + event_no + ' comp: North ' + vmodel)
        plt.show()
        
        plt.close(3)
        plt.figure(3,figsize=(10,8))
        plt.xlim(-start_buff,end_buff)
        plt.ylim(min_dist,max_dist)
        count = 0
        
        for tr in st_synZ_good:
            dist_offset = tr.stats.distance
            ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)#    These lines used to cause a crash in Spyder
            plt.plot(ttt, (tr.data * plot_fac * tr.stats.distance) + dist_offset, color = 'red')
            plt.text(s = tr.stats.network + ' ' + tr.stats.station ,x = end_buff*0.95,y = dist_offset
            + max_dist*0.015, color = 'black')  #label traces

        for jstr in js_synZ_good:   
            if flip_z == False:
                dist_offset = st_synE_good[count].stats.distance
                ttt = np.arange(0, 400, 0.05)
                #ttt = np.arange(len(jstr)) * (0.05) + (st_synE_good[count].stats.starttime - t1 - 50)#    These lines used to cause a crash in Spyder
                plt.plot(ttt, (jstr.data[0:8000] * js_plot_fac * st_synE_good[count].stats.distance) + dist_offset, color = 'blue')
                plt.text(s = st_synE_good[count].stats.network + ' ' + st_synE_good[count].stats.station ,x = end_buff*0.95,y = dist_offset
                + max_dist*0.015, color = 'black')  #label traces
        
                print(count)
                print(dist_offset)
                print(jstr.data[0:8000])
                print(jstr.stats.get('station'))
                print(st_synE_good[count].stats.station)
                print(st_synE_good[count].stats.distance)
                
                count += 1
            else:
                dist_offset = st_synE_good[count].stats.distance
                ttt = np.arange(0, 400, 0.05)
                #ttt = np.arange(len(jstr)) * (0.05) + (st_synE_good[count].stats.starttime - t1 - 50)#    These lines used to cause a crash in Spyder
                plt.plot(ttt, (jstr.data[0:8000] * -1 * js_plot_fac * st_synE_good[count].stats.distance) + dist_offset, color = 'blue')
                plt.text(s = st_synE_good[count].stats.network + ' ' + st_synE_good[count].stats.station ,x = end_buff*0.95,y = dist_offset
                + max_dist*0.015, color = 'black')  #label traces
        
                print(count)
                print(dist_offset)
                print(jstr.data[0:8000])
                print(jstr.stats.get('station'))
                print(st_synE_good[count].stats.station)
                print(st_synE_good[count].stats.distance)
                
                count += 1
            
        print(st_synE_good[0].stats.delta)
        print(st_synE_good[1].stats.delta)
    
        plt.xlabel('Time (s)')
        plt.ylabel('Epicentral distance from event (km)')
        plt.title('Jon Stewart (blue) vs Ricardo Taborda (red) Syn ' + date_label + ' ' + event_no + ' comp: Z ' + vmodel)
        plt.show()
    else:
        for tr in st_synE_good:
            dist_offset = tr.stats.distance
            ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)#    These lines used to cause a crash in Spyder
            plt.plot(ttt, (tr.data * plot_fac * tr.stats.distance) + dist_offset, color = 'red')
            plt.text(s = tr.stats.network + ' ' + tr.stats.station ,x = end_buff*0.95,y = dist_offset
            + max_dist*0.015, color = 'black')  #label traces
     
        plt.xlabel('Time (s)')
        plt.ylabel('Epicentral distance from event (km)')
    #    plt.title(fname1[8:18] + ' vs ' + fname2[8:18])
        plt.title('Ricardo Syn ' + date_label + ' ' + event_no + ' E-black N-green Z-red of vmodel ' + vmodel)
        
        plt.close(10)
        plt.figure(10,figsize=(10,8))
        plt.xlim(-start_buff,end_buff)
        plt.ylim(min_dist,max_dist)
        count = 0
        for jstr in js_synE_good:   
            dist_offset = st_synE_good[count].stats.distance
            ttt = np.arange(0, 400, 0.05)
            #ttt = np.arange(len(jstr)) * (0.05) + (st_synE_good[count].stats.starttime - t1 - 50)#    These lines used to cause a crash in Spyder
            plt.plot(ttt, (jstr.data[0:8000] * js_plot_fac * st_synE_good[count].stats.distance) + dist_offset, color = 'blue')
            plt.text(s = st_synE_good[count].stats.network + ' ' + st_synE_good[count].stats.station ,x = end_buff*0.95,y = dist_offset
            + max_dist*0.015, color = 'black')  #label traces
    
            print(count)
            print(dist_offset)
            print(jstr.data[0:8000])
            print(jstr.stats.get('station'))
            print(st_synE_good[count].stats.station)
            print(st_synE_good[count].stats.distance)
            
            count += 1
            
        print(st_synE_good[0].stats.delta)
        print(st_synE_good[1].stats.delta)
    #    for tr in st_dat_good:
    #        dist_offset = tr.stats.distance # km
    #        plt.text(s = tr.stats.network + ' ' + tr.stats.station ,x = end_buff*0.95,y = dist_offset
    #                + max_dist*0.015, color = 'black')  #label traces
    #
    #    for tr in st_synN_good:
    #        dist_offset = tr.stats.distance
    #        ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)#    These lines used to cause a crash in Spyder
    #        plt.plot(ttt, (tr.data * plot_fac * tr.stats.distance) + dist_offset, color = 'green')
    #
    #    for tr in st_synZ_good:
    #        dist_offset = tr.stats.distance
    #        ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)#    These lines used to cause a crash in Spyder
    #        plt.plot(ttt, (tr.data * plot_fac * tr.stats.distance) + dist_offset, color = 'red')
    
    
        plt.xlabel('Time (s)')
        plt.ylabel('Epicentral distance from event (km)')
    #    plt.title(fname1[8:18] + ' vs ' + fname2[8:18])
        plt.title('Jon Stewart Syn ' + date_label + ' ' + event_no + ' E-black N-green Z-red of vmodel ' + vmodel)
        plt.show()
    
        #  Save processed files
    #    fname1 = 'Pro_Files/HD' + date_label1 + 'sel.mseed'
    #    fname2 = 'Pro_Files/HD' + date_label2 + 'sel.mseed'
    #    st1good.write(fname1,format = 'MSEED')
    #    st2good.write(fname2,format = 'MSEED')
    
        elapsed_time_wc = time.time() - start_time_wc
        print('This job took ' + str(elapsed_time_wc) + ' seconds')
        #os.system('say "Done"')