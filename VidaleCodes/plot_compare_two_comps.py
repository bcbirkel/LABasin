#!/usr/bin/env python
# program to compare two sets of seismograms for LA basin motions
# John Vidale 9/2020

def p_compare2(event_no, traces1, traces2, start_buff, end_buff, taper, taper_frac,
            plot_scale_fac, filt, freq_min, freq_max, min_dist, max_dist,
            norm_each, dist_norm, basin_width, basin, CI_only, component1, component2,
            fig_inc):
    # component 1 is E, 2 is N, 3 is Z
    # H is cvmhy, H_simp is cvmhn
    # S is cvms426-223, S_simp is cvms400-100
    # norm_each == 1 scales the peak of each trace to the same size
    #    otherwise, the amplitude is normalized to 10 km, with the scale set with plot_scale_fac

    from obspy import UTCDateTime
    from obspy import Stream
    from obspy import read
    from obspy.geodetics import gps2dist_azimuth
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    import time

    import sys # don't show any warnings
    import warnings

    if not sys.warnoptions:
        warnings.simplefilter("ignore")

    start_time_wc = time.time()

    # Box limits
    Lon_min=-118.8
    Lat_min=33.5
    Lon_max=-117.5
    Lat_max=34.6

    if event_no == '15481673':  # no 29
        event_name = 'LaHabra'
    elif event_no == '10410337':  # no 14
        event_name = 'Inglewood'
    elif event_no == '14383980':  # no 29
        event_name = 'ChinoHills'
    elif event_no == '14312160':  # no 29
        event_name = 'Chatsworth'
    elif event_no == '9703873':  # no 29
        event_name = 'BeverlyHills'
    else:
        print('Number does not correspond to a valid event')
        sys.exit()

    verbose = False # print every distance of a station

    if component2 == 'same':
        component2 = component1

    #%% find event details for origin time, lat, lon
    ev_file = '/Users/vidale/Documents/PyCode/LAB/Compare/event_list.txt'
    file_ev = open(ev_file, 'r')
    for line in file_ev:           # pull numbers off all the lines
        split_line = line.split()
        event = split_line[0]
        if event == event_no:
            ev_lat       = float(split_line[3])
            ev_lon       = float(split_line[2])
            t1           = UTCDateTime(split_line[5])
            date_label  = split_line[5][0:10]
            year1        = split_line[5][0:4]
    print(event_name + ': ev_no ' + event_no + ' ' + str(t1) + ' Lat-Lon ' + str(ev_lat) + ' ' + str(ev_lon))

    #%% find event details for origin time, lat, lon
    badt_file = '/Users/vidale/Documents/PyCode/LAB/Compare/bad_trace.txt'
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
    sta_file = '/Users/vidale/Documents/PyCode/LAB/Compare/stations.txt'
    file_st = open(sta_file, 'r')
    line = file_st.readline()      # read first line to skip header information
    lines = file_st.readlines()
    # print(str(len(lines)) + ' stations read from ' + sta_file)

    # Load station coords into arrays, many more stations than used
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
#        print('Event ' + str(ev_lat) + ' ' + str(ev_lon) + ' station ' + split_line[4] + ' ' + split_line[5] + ' distance ' + str(distance[0]))
        st_dist.append(distance[0]/1000.) # distance
        st_az.append(distance[1]) # azimuth
        st_baz.append(distance[2]) # back-azimuth
    print('number of stations in list is ' + str(len(st_num)))

    #%% Load data and synthetic waveforms
    sgrams1 = Stream()
    if traces1 == 'Data':
        fname_datZ = '/Users/vidale/Documents/PyCode/LAB/Compare/' + event_name + '/Data/' + component1 + '.mseed'
    elif traces1 == 'CVM_S4':
        fname_datZ = '/Users/vidale/Documents/PyCode/LAB/Compare/' + event_name + '/CVM_S4/' + component1 + '.mseed'
    elif traces1 == 'CVM_H':
        fname_datZ = '/Users/vidale/Documents/PyCode/LAB/Compare/' + event_name + '/CVM_H/' + component1 + '.mseed'
    sgrams1 = read(fname_datZ)

    # if traces1 == 'Ricardo_syn':
    #     lenH = len(sgrams1)     # correct Hercules units from from meters to cm
    #     for ii in range(lenH):  # assume trace count of all synthetic files is same
    #         sgrams1[ii].data      = 100 * sgrams1[ii].data

    len1_in = len(sgrams1)
    print(traces1 + ' has ' + str(len1_in) + ' traces')
    if len(sgrams1) == 0:
        print('No stations found in traces1, quit now!')
        sys.exit()
    print('1st gather, 1st data trace has : ' + str(len(sgrams1[0].data)) + ' time pts ')
    print('1st gather, 1st trace starts at ' + str(sgrams1[0].stats.starttime) + ', event at ' + str(t1))

    if traces2 != 'no':
        sgrams2 = Stream()
        if traces2 == 'Data':
            fname_datZ = '/Users/vidale/Documents/PyCode/LAB/Compare/' + event_name + '/Data/' + component2 + '.mseed'
        elif traces2 == 'CVM_S4':
            fname_datZ = '/Users/vidale/Documents/PyCode/LAB/Compare/' + event_name + '/CVM_S4/' + component2 + '.mseed'
        elif traces2 == 'CVM_H':
            fname_datZ = '/Users/vidale/Documents/PyCode/LAB/Compare/' + event_name + '/CVM_H/' + component2 + '.mseed'
        sgrams2 = read(fname_datZ)

    # if traces2 == 'Ricardo_syn':
    #     lenH = len(sgrams2)     # correct Hercules units from from meters to cm
    #     for ii in range(lenH):  # assume trace count of all synthetic files is same
    #         sgrams2[ii].data      = 100 * sgrams2[ii].data

    if traces2 != 'no':
        len2_in = len(sgrams2)
        print(traces2 + ' has ' + str(len2_in) + ' traces')
        if len(sgrams2) == 0:
            print('No stations found in traces1, quit now!')
            sys.exit()
        print('2nd gather, 1st trace has : ' + str(len(sgrams2[0].data)) + ' time pts ')
        print('2nd gather, 1st trace starts at ' + str(sgrams2[0].stats.starttime) + ', event at ' + str(t1))

    #%%  Keep only data that is not on list of bad traces
    # either individual components or A for all components
    sgrams1_good = Stream()
    for tr in sgrams1: # examine traces one by one
        do_write = 1
        for ii in range(len(badt_event)):
            if event_no == badt_event[ii] and tr.stats.station == badt_station[ii]: # find station in inventory
                if badt_compo[ii] == 'Z' or badt_compo[ii] == 'A':
                    do_write = 0
        if do_write == 1:
            sgrams1_good += tr
    print('After rejecting listed bad traces, ' + str(len(sgrams1_good)) + ' out of ' + str(len1_in) + ' traces remain in 1st gather.')

    if traces2 != 'no' or '':
        sgrams2_good = Stream()
        for tr in sgrams2: # examine traces one by one
            do_write = 1
            for ii in range(len(badt_event)):
                if event_no == badt_event[ii] and tr.stats.station == badt_station[ii]: # find station in inventory
                    if badt_compo[ii] == 'Z' or badt_compo[ii] == 'A':
                        do_write = 0
            if do_write == 1:
                sgrams2_good += tr
        print('After rejecting listed bad traces, ' + str(len(sgrams2_good)) + ' out of ' + str(len2_in) + ' traces remain in 2nd gather.')

    if len(sgrams1_good) == 0:
        print('No common stations found, quit now!')
        sys.exit()
#%%  Keep only common traces
    if traces2 != 'no':
        sgrams1_comm = Stream()
        sgrams2_comm = Stream()
        for tr1 in sgrams1_good: # examine traces one by one
            found_it = 0
            if tr1.stats.station in st_name:  # find station in station list
                found_it = 1
                ii = st_name.index(tr1.stats.station)
                for tr2 in sgrams2_good: # examine traces one by one
                    if tr2.stats.station == tr1.stats.station:  # find station in station list
                        in_CI = (tr1.stats.network == 'CI' or tr2.stats.network == 'CI')
                        if in_CI == True or CI_only == False:
                            sgrams1_comm += tr1
                            sgrams2_comm += tr2
            if found_it == 0:
                print('Not on station list! station ' + tr1.stats.station + ' network ' + tr1.stats.network)
        print('1st, 2nd sets of seismograms have ' + str(len(sgrams1_comm)) + ' and ' + str(len(sgrams2_comm)) + ' common stations')
    else:
        sgrams1_comm = Stream()
        sgrams1_comm = sgrams1_good.copy()

#%% Select data and syn, only need one, don't need 2nd gather for map
    # select data by distance (and azimuth?), and cull synthetics to match data

    sgrams1_select = Stream()
    for tr in sgrams1_comm: # examine traces one by one
        if (tr.stats.network == 'CI' or CI_only == False or traces1 == 'CVM_H' or traces1 == 'CVM_S4'):
            found_it = 0
            for ii in range(len(st_name)):  # find matching entry in station roster
                if tr.stats.station == st_name[ii]: # find station in inventory
                    found_it = 1
                    flat = float(st_lat[ii])
                    flon = float(st_lon[ii])
                    if st_dist[ii] < max_dist and (flat < Lat_max and flat > Lat_min and
                                                   flon < Lon_max and flon > Lon_min):
                        # exclude stations too far away or outside the map box
                        # basin is roughly within 15 km of these 3 stations
                        distance = gps2dist_azimuth( 33.99053, -118.36171, float(flat), float(flon)) # station BHP
                        dist1 = distance[0]/1000  # convert m to km
                        distance = gps2dist_azimuth( 33.88110, -118.17568, float(flat), float(flon)) # station LTP
                        dist2 = distance[0]/1000
                        distance = gps2dist_azimuth( 33.80776, -117.98116, float(flat), float(flon)) # station BRE
                        dist3 = distance[0]/1000
    #                    print(tr.stats.station + ' ' + str(dist1) + ' ' + str(dist2) + ' ' + str(dist3))
    #                    keep stations only within X km of basin axis
                        if basin == False or (dist1 < basin_width) or (dist2 < basin_width) or (dist3 < basin_width):
    #                        print('selected: ' + tr.stats.station)
                            tr.stats.distance = st_dist[ii]
                            tr.stats.stla = st_lat[ii]
                            tr.stats.stlo = st_lon[ii]
                            sgrams1_select += tr
            # if found_it == 0:
            #     print('Not on station list! station ' + tr.stats.station + ' network ' + tr.stats.network)
    print('Sgrams 1: Within ' + str(max_dist) + ' km distance and optional basin culling, ' + str(len(sgrams1_select)) + ' traces remain for plotting.')

    if traces2 != 'no':
        sgrams2_select = Stream()
        for tr in sgrams2_comm: # examine traces one by one
            if (tr.stats.network == 'CI' or CI_only == False or traces2 == 'CVM_H' or traces2 == 'CVM_S4'):
                found_it = 0
                for ii in range(len(st_name)):  # find matching entry in station roster
                    if tr.stats.station == st_name[ii]: # find station in inventory
                        found_it = 1
                        flat = float(st_lat[ii])
                        flon = float(st_lon[ii])
                        if st_dist[ii] < max_dist and (flat < Lat_max and flat > Lat_min and
                                                       flon < Lon_max and flon > Lon_min):
                            # exclude stations too far away or outside the map box
                            # basin is roughly within 15 km of these 3 stations
                            distance = gps2dist_azimuth( 33.99053, -118.36171, float(flat), float(flon)) # station BHP
                            dist1 = distance[0]/1000  # convert m to km
                            distance = gps2dist_azimuth( 33.88110, -118.17568, float(flat), float(flon)) # station LTP
                            dist2 = distance[0]/1000
                            distance = gps2dist_azimuth( 33.80776, -117.98116, float(flat), float(flon)) # station BRE
                            dist3 = distance[0]/1000
        #                    print(tr.stats.station + ' ' + str(dist1) + ' ' + str(dist2) + ' ' + str(dist3))
                            if basin == False or (dist1 < basin_width) or (dist2 < basin_width) or (dist3 < basin_width):  # keep stations only within X km of basin axis
        #                        print('selected: ' + tr.stats.station)
                                tr.stats.distance = st_dist[ii]
                                tr.stats.stla = st_lat[ii]
                                tr.stats.stlo = st_lon[ii]
                                sgrams2_select += tr
                # if found_it == 0:
                #     print('Not on station list! station ' + tr.stats.station + ' network ' + tr.stats.network)
        print('Sgrams 2: Within ' + str(max_dist) + ' km distance and optional basin culling, ' + str(len(sgrams2_select)) + ' traces remain for plotting.')
    if len(sgrams1_select) == 0:
        print('No common stations found, quit now!')
        sys.exit()

    #%%  detrend, taper, filter
    if taper:
        sgrams1_select.taper( taper_frac)
        if traces2 != 'no':
            sgrams2_select.taper( taper_frac)
    if filt:
        sgrams1_select.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
        if traces2 != 'no':
            sgrams2_select.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
    if taper:
        sgrams1_select.taper( taper_frac)
        if traces2 != 'no':
            sgrams2_select.taper( taper_frac)

    #%% plot traces, different windows in case multiple components are to be compared
    if component1 == 'Z':
        fig_index = 13 + fig_inc
    elif component1 == 'N':
        fig_index = 14 + fig_inc
    elif component1 == 'E':
        fig_index = 15 + fig_inc
    elif component1 == 'R':
        fig_index = 16 + fig_inc
    elif component1 == 'T':
        fig_index = 17 + fig_inc

    plt.close(fig_index)
    plt.figure(fig_index,figsize=(10,8))
    plt.xlim(start_buff,end_buff)
    plt.ylim(min_dist,max_dist)

    # find max of absolute amplitude
    maxD1 = 0
    for tr in sgrams1_select:
        tr_max = max(abs(tr.data))
        if tr_max > maxD1:
            maxD1 = tr_max

    if traces2 != 'no':
        maxD2 = 0
        for tr in sgrams2_select:
            tr_max = max(abs(tr.data))
            if tr_max > maxD2:
                maxD2 = tr_max
        maxD = max(maxD1, maxD2)
#    print('Max data and H, S, H_simp, and S_simp synthetics are ' + str(maxD) + '  ' + str(maxSH) + '  ' + str(maxSS)+ '  ' + str(maxSH_simp)+ '  ' + str(maxSS_simp))

    # find max normalized to 10 km distance (i.e., amp divided by distance, assumes 1/R amp fall-off)
    maxD_N = 0
    for tr in sgrams1_select:
        if verbose == True:
            print(f'Distance is {tr.stats.distance:6.4f}' + '  ' + tr.stats.station)
        tr_max = max(abs(tr.data))*tr.stats.distance
        if tr_max > maxD_N:
            maxD_N = tr_max

    print(f'Peak amp of trace1 is {maxD1:6.4f}')
    if traces2 != 'no':
        print(f'Peak amp of trace2 is {maxD2:6.4f}')
        print(f'Peak amp of both traces is {maxD:6.4f}')
    print(f'Normed amp of data is {maxD_N:6.4f}')
    plot_fac = plot_scale_fac * (max_dist - min_dist) / maxD_N

    for tr in sgrams1_select:
        dist_offset = tr.stats.distance # km
        if dist_offset > min_dist and dist_offset < max_dist:
            pnum = max(abs(tr.data))
            printee2 = f'amp {pnum:8.3f}  {tr.stats.network}  {tr.stats.station}'
            plt.text(s = printee2,x = end_buff*0.8,y = dist_offset
                    + (max_dist - min_dist)*0.02, color = 'black')  #label traces and note amplitude (cm, cm/s, cm/s^^2)
            if traces1 == 'CVM_H' or traces1 == 'CVM_S4':
                tr.stats.starttime = t1
            ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)
            color_p = 'red'
            if traces2 == 'no':
                color_p = 'black'
            if norm_each:
                plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
                    - tr.data.min()) + dist_offset, color = color_p)
            elif dist_norm == True:
                plt.plot(ttt, (tr.data * plot_fac * tr.stats.distance) + dist_offset, color = color_p)
            else:
                plt.plot(ttt, (tr.data * plot_fac) + dist_offset, color = color_p)

    if traces2 != 'no':
        for tr in sgrams2_select:
            dist_offset = tr.stats.distance
            if dist_offset > min_dist and dist_offset < max_dist:
                if traces2 == 'CVM_H' or traces2 == 'CVM_S4':
                    tr.stats.starttime = t1
                ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)
                if norm_each:
                    plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
                        - tr.data.min()) + dist_offset, color = 'green')
                elif dist_norm == True:
                    plt.plot(ttt, (tr.data * plot_fac * tr.stats.distance) + dist_offset, color = 'green')
                else:
                    plt.plot(ttt, (tr.data * plot_fac) + dist_offset, color = 'green')

    plt.xlabel('Time (s)')
    plt.ylabel('Epicentral distance from event (km)')
#    plt.title(fname1[8:18] + ' vs ' + fname2[8:18])
    if traces2 != 'no':
        plt.title(date_label + ' red ' + traces1 + ' ' + component1 + ' green ' + traces2 + ' ' + component2)
    else:
        plt.title(date_label + ' ' + traces1 + ' ' + component1)
    plt.show()

    elapsed_time_wc = time.time() - start_time_wc
    print('This job took ' + str(elapsed_time_wc) + ' seconds')
    os.system('say "Done"')