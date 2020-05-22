#!/usr/bin/env python
# input is set of LAB strong motion traces
# this program tapers and filters
# This programs deals with a single event.
# John Vidale 5/2019

def plot_sgrams(select_event=14, filt = 1, freqmin = 0.05, freqmax = 1.0, component = 'Z',
                norm = 1, plot_scale_fac = .03, start_buff = 5, end_buff = 20,
                min_plot_dist = 6.45, max_plot_dist = 6.55):
    import os
    from obspy import Stream
    from obspy import read_events, read
    import matplotlib.pyplot as plt
    from obspy.geodetics import gps2dist_azimuth
    import numpy as np

    os.environ['PATH'] += os.pathsep + '/usr/local/bin'
    os.chdir('/Users/bcbirkel/Documents/Research/LABasin/PyCodeLAB/Mseed')

#    component = 'Z'
#    select_event = 14
#    norm = 1
#    plot_scale_fac = .03
#    max_plot_dist = 6.55
#    min_plot_dist = 6.45
#    start_buff = 5
#    end_buff = 20

#%% stations to exclude
    basin_sta = ['SMF2','LAF','USC','LCG','WTT','COO','WTT2','LGB','LTP','STS','DLA','BRE','OGC','FUL','LLS','SAN','LBW1','LBW2','5499','CRF']
    bad_sta = ['BVH','LAX','KIK','CAC','24851','14405','24400','WTT2','PDR','HLL','5402','GVR','MIKB','RHC2','PEM','BRE','GSA']
#    bad_sta   = ['PASC','DJJ','RHC','OGC','SMF','CAC','DSN5','14405','24400','KIK','MIKB','PEM']
    bad_loc   = ['01','02','10','30','40','41','B0']
#   takenout_sta = ['BHP','PDR']


    # Load station coords into arrays
    sta_file = '/Users/bcbirkel/Documents/Research/LABasin/PyCodeLAB/sta_list.txt'
    with open(sta_file, 'r') as file:
        lines = file.readlines()
    print(str(len(lines)) + ' stations read from ' + sta_file)

    station_index = range(len(lines))
    st_names = []
    st_lats  = []
    st_lons  = []
    for ii in station_index:
        line = lines[ii]
        split_line = line.split()
        st_names.append(split_line[0])
        st_lats.append( split_line[1])
        st_lons.append( split_line[2])

    # Load seismograms into array
    st = Stream()
    fname = '/Users/bcbirkel/Documents/Research/LABasin/PyCodeLAB/Mseed/IRISdata/event' + str(select_event) + '/event' + str(select_event) + component + '_chosen.mseed'
    st = read(fname)


    if select_event > 15:
        fname_inv = 'IRISdata/LAB.QUAKEML2'
        LAB = read_events(fname_inv, format='QUAKEML')
        if select_event == 16:
            t = LAB[1].origins[0].time
            ref_lat = LAB[1].origins[0].latitude
            ref_lon = LAB[1].origins[0].longitude
        elif select_event == 17:
            t = LAB[5].origins[0].time
            ref_lat = LAB[5].origins[0].latitude
            ref_lon = LAB[5].origins[0].longitude
        elif select_event == 18:
            t = LAB[14].origins[0].time
            ref_lat = LAB[14].origins[0].latitude
            ref_lon = LAB[14].origins[0].longitude
        elif select_event == 19:
            t = LAB[13].origins[0].time
            ref_lat = LAB[13].origins[0].latitude
            ref_lon = LAB[13].origins[0].longitude
    else:
        fname_inv = 'IRISdata/LAB.QUAKEML'
        LAB = read_events(fname_inv, format='QUAKEML')
        t = LAB[select_event].origins[0].time
        ref_lat = LAB[select_event].origins[0].latitude
        ref_lon = LAB[select_event].origins[0].longitude

    print(f'source lat is {ref_lat:.2f}  source lon is {ref_lon:.2f}')

    st.detrend()
    st.taper(0.1)
    if filt == 1:
        st.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=False)

    # make figure
#    fig = plt.figure(4)
#    plt.title('Trace viewer')
#    plt.xlabel('seconds')
#    plt.ylabel('stations')
#    st.plot(size=(1500,1400), equal_scale=False)

    print('Made it to here! Timedate:  ' + str(t) + 'trace 0 datetime:  ' + str(st[0].stats.starttime))

    if component == 'Z':
        plt.close(3)
        fig_index = 3
    elif component == 'N':
        plt.close(4)
        fig_index = 4
    elif component == 'E':
        plt.close(5)
        fig_index = 5
    else:
        print('Invalid component')
    plt.figure(fig_index,figsize=(15,10))
    plt.xlim(-start_buff,end_buff)
    plt.ylim(min_plot_dist,max_plot_dist)

    max_amp = 0
    for tr in st:
        if (tr.stats.station not in bad_sta):
            if (tr.stats.location not in bad_loc):
                if tr.stats.station in st_names:
                    sta_index = st_names.index(tr.stats.station)
                    stalat = float(st_lats[sta_index])
                    stalon = float(st_lons[sta_index]) # look up lat & lon again to find distance
                    distance = gps2dist_azimuth(stalat,stalon,ref_lat,ref_lon) # Get traveltimes again, hard to store
                max_local = np.max(np.abs(tr.data))*distance[0]
                if max_local > max_amp:
                    max_amp = max_local
    print('Also made it to here!')
    for tr in st:
        if (tr.stats.station not in bad_sta):
            if (tr.stats.location not in bad_loc):
                if tr.stats.station in st_names:
                    sta_index = st_names.index(tr.stats.station)
                    stalat = float(st_lats[sta_index])
                    stalon = float(st_lons[sta_index]) # look up lat & lon again to find distance
                    distance = gps2dist_azimuth(stalat,stalon,ref_lat,ref_lon) # Get traveltimes again, hard to store
                else:
                    print(tr.stats.station + ' is not in station list')

                dist_offset = distance[0]/1000 # correct to km
#                print(f'{tr.stats.station:s} {tr.stats.channel:s} distance is {dist_offset:.1f} slat is {stalat:.3f}  slon is {stalon:.3f}')
#                print(str(tr.stats.starttime) + '  ' + str(t))
                ttt = tr.stats.starttime - t + np.arange(len(tr.data)) * tr.stats.delta

                if tr.stats.station in basin_sta: # color basin stations red
                    if norm == 1:  # trace normalize amplitude
                        plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
                            - tr.data.min()) + dist_offset, color = 'red')
                    else:  # use common amplitude scale, multiply by distance
                        plt.plot(ttt, (tr.data - np.median(tr.data))*(2*distance[0]*plot_scale_fac /max_amp)
                         + dist_offset, color = 'red')
                    plt.text(s = tr.stats.station,x = end_buff*0.95,y = dist_offset
                            + (max_plot_dist-min_plot_dist)*0.015, color = 'red')  #label traces
                else:  # color stations outside basin black
                    if norm == 1:  # trace normalize amplitude
                        plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
                            - tr.data.min()) + dist_offset, color = 'black')
                    else:
                        plt.plot(ttt, (tr.data - np.median(tr.data))*(2*distance[0]*plot_scale_fac /max_amp)
                         + dist_offset, color = 'black')
                    plt.text(s = tr.stats.station,x = end_buff*0.95,y = dist_offset
                        + (max_plot_dist-min_plot_dist)*0.015, color = 'black')  #label traces

    date_title = str(t)
    plt.title('Beverly Hills, Mag 4.2 - ' + date_title[0:19] + ' ' + '- ' + component + '-component acceleration')
#    plt.title('Event ' + str(select_event) + ' ' + date_title[0:19] + ' ' + ' ' + component + ' acc.')
    plt.xlabel('Time (s)')
    plt.ylabel('Distance (km)')
