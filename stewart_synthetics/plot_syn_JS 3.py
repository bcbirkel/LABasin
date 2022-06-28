#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 15:30:49 2019

@author: bcbirkel
"""

# Plot ASCII files to see Jon Stewart's synthetics
# units are cm/s
# Brianna Birkel - 11/2019
#     from __future__ import print_function

def plot_syn_JS(eventfile, event_title, vmodel, stations, sta_names):
    
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    import sys
    from obspy import read_events, UTCDateTime, read, Trace, Stream
    from obspy.geodetics import gps2dist_azimuth
    from obspy.core.util.attribdict import AttribDict
    from obspy.core.trace import Stats
    from obspy.io.sac import SACTrace

    
# %% VARIABLES
#    freqmin = 0.1
#    freqmax = 0.5
#    start_buff = 0
#    end_buff = 150
#    min_plot_dist = 0
#    max_plot_dist = 100
#    scale_fact = 1
    
    count = 0
    # crap for file naming to work right
    if vmodel == 'S':
        mod = 'si'
        mod_title = 'S4'
    elif vmodel == 'H':
        mod = 'h'
        mod_title = 'H'
    
# %% SETUP 
    dir = '/Users/bcbirkel/Documents/Research/LABasin/stewart_synthetics/CVM-' + \
    mod_title + '/' + eventfile + '_' + event_title + '/Vel/'
    
    print(dir)
    os.chdir(dir)
    
    event_files = os.listdir(dir)
    num_files = len(event_files)
    print(num_files) 
    print(event_files)
    #stastr = str(stations)
    #print(stastr)
    
    stT_all = Stream()
    stE_all = Stream()
    stN_all = Stream()
    stZ_all = Stream()
    
# %% initialize file/variables    

    for ii in range(num_files):
        fname = event_files[ii]
        for sta in range(len(stations)):
            os.chdir(dir)
            stastr=str(stations[sta])
            sta_name = str(sta_names[sta])
            if ('.bbp' and stastr) in fname:
                #print(fname)
                #print(os.getcwd())
                with open(fname, 'r') as file:
                    lines = file.readlines()
                    #print('check1')
                     #print(sta_file)
                     #print(str(len(lines)) + ' stations read from ' + file)
                     #INTIALIZE
                    time = []
                    comp_N = []
                    comp_E = []
                    comp_Z = []
                    all_text = '' 
                    # Remove comments, put txt file into lists for time/comps
                    num_lines = range(len(lines))
                    #print('check2')
                
# %% put data into list
                    for line in open(fname):
                        if not line.startswith("#"):
                            all_text = all_text + line
                            #print(all_text)
                            line = line.split('\t')
                            #print(line)
#                            for num in range(0,4):
#                                line[num] = float(line[num])
                            time.append(line[0])
                            comp_N.append(line[1])
                            comp_E.append(line[2])
                            comp_Z.append(line[3])
                        if "longitude" in line:
                            line = line.split()
                            sta_lon = line[2]
                        if "latitude" in line:
                            line = line.split()
                            sta_lat = line[2]
                            
    
                    #scale up synthetics - unused currently
#                    comp_N = [scale_fact*x for x in comp_N] 
#                    comp_E = [scale_fact*x for x in comp_E] 
#                    comp_Z = [scale_fact*x for x in comp_Z] 
                    

    # %% fill header            
     
                    # Fill header attributes
                    print('station code going into header is ' + sta_name)
                    stats = Stats()
                    header={'delta': (float(time[1])-float(time[0])), 'npts': len(time), 'sampling_rate': 1/(float(time[1])-float(time[0])),'kstnm': sta_name, \
                            'stla': sta_lat, 'stlo': sta_lon}
                    #stats.network= 'CI'
                    #stats.channel = 'ZNE'
                    stats.npts = len(time)
                    stats.sampling_rate = 1/(float(time[1])-float(time[0]))
                    stats.station = sta_name
                    #stats.starttime = time[0]
                    #stats.endtime = time[range(time)]
                    stats.stla = sta_lat
                    stats.stlo = sta_lon
                    print(stats.station)

                    # set current time
                    #stats['starttime'] = UTCDateTime()
                    #***GET TIME FROM STA FILE
                    
    # %% set time/component lists as NumPy arrays
                    time_ms = np.asarray(time)
                    N_ms = np.asarray(comp_N)
                    E_ms = np.asarray(comp_E)
                    Z_ms = np.asarray(comp_Z)
                    
                    #type(N_ms)
                    print(E_ms)
                    #print('check')
                    
    # %% Stream objects for time + each comp
                    
                    os.chdir('/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/compare_syn')
                    
                    #stT = Stream([Trace(data=time_ms, header=stats)])
                    #stT_all += stT
                    stT = SACTrace(data=time_ms, **header)
                    
                    #stN = Stream([Trace(data=N_ms, header=stats)])
                    stN.detrend()
#                    stN.taper(0.1)
#                    stN.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=False)
                    #stN_all += stN
                    stN = SACTrace(data=N_ms, **header)
                    
                    #stE = Stream([Trace(data=E_ms, header=stats)])
                    #stE.detrend()
                    #stE.taper(0.1)
                    #stE.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=False)
                    #stE_all += stE
                    stE = SACTrace(data=E_ms, **header)
                    
                    
                    #stZ = Stream([Trace(data=Z_ms, header=stats)])
                    stZ.detrend()
                    #stZ.taper(0.1)
                    #stZ.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=False)
                    #stZ_all += stZ
                    #stE.write("Z_bp-" + stastr + "-" + event_title + ".mseed", format='MSEED', reclen=256)
                    stZ = SACTrace(data=Z_ms, **header)
                    
                    # Show that it worked, convert NumPy character array back to string
                    #st1 = read("data.mseed")
                    #print(st1[0].data.tostring())
                    #print(stE)
                    count += 1
    
        # %% convert to SAC and save files
                    os.chdir('/Users/bcbirkel/Documents/Research/LABasin/Graves_SAC/working')
                    stN.write("Syn_" + event_title + "_" + sta_name + "_N.SAC", format='SAC')
                    stE.write("Syn_" + event_title + "_" + sta_name + "_E.SAC", format='SAC')
                    stZ.write("Syn_" + event_title + "_" + sta_name + "_Z.SAC", format='SAC')
    
    
    print(os.getcwd())
    print("for loop ran " + str(count) + " times")
    print(stE_all)
    print(stE)
    print("check here")
    print(eventfile)
    #print(stN.stats.code)
    os.chdir('/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/compare_syn')
    #stT_all.write("time-" + eventfile + ".mseed", format='MSEED', reclen=256, ) 
    stN_all.write("N_bp-" + eventfile +  ".mseed", format='MSEED', reclen=256)
    stE_all.write("E_bp-" + eventfile +  ".mseed", format='MSEED', reclen=256)
    stZ_all.write("Z_bp-" + eventfile +  ".mseed", format='MSEED', reclen=256)
            
    return time, stN, stE, stZ

            
    
    """                
    # %% Plotting
                    
                    if numTr == 'one':
                        if component == 'Z':
                            plt.close(3)
                            fig_index = 3
                            comp = stZ[0]
                        elif component == 'N':
                            plt.close(4)
                            fig_index = 4
                            comp = stN[0]
                        elif component == 'E':
                            plt.close(5)
                            fig_index = 5
                            comp = stE[0]
                        elif component == 'all':
                            plt.close(6)
                            fig_index = 6
                        else:
                            print('Invalid component')
                        plt.figure(fig_index,figsize=(15,10))
                        plt.xlim(-start_buff,end_buff)
                        #plt.ylim(min_plot_dist,max_plot_dist)
                        if component == 'all':
                            plt.plot(time,stN[0],stE[0],stZ[0])
                        else:
                            plt.plot(time,comp)
                    elif numTr == 'all':
                        for ii in range(7,10):
                            if ii == 7:
                                comp = stN[0]
                                name = 'N'
                            elif ii == 8:
                                comp = stE[0]
                                name = 'E'
                            elif ii == 9:
                                comp = stZ[0]
                                name = 'Z'
                            plt.close(ii)
                            plt.figure(ii,figsize=(15,10))
                            plt.xlim(-start_buff,end_buff)
                            plt.plot(time,comp)
                            plt.title(event_title + ' - ' + name + '-component velocity, model = ' 
                                      + vmodel + ', sta = ' + station)
                            plt.xlabel('Time (s)')
                            plt.ylabel('Velocity (cm/s)')
                    else:
                        print('numTr must be one or all')
                        
        return time, stN[0], stE[0], stZ[0]
            


# list to str
time_str = ' '.join(map(str, time))
comp_N_str = ' '.join(map(str, comp_N))
comp_E_str = ' '.join(map(str, comp_E))
comp_Z_str = ' '.join(map(str, comp_Z))
# =============================================================================
# for row in sta_file:
#     line = row.split('\n')
#     line = ' '.join(line)
#     if not line.startswith("#"):
#         line = line.split('\t')
#         for num in range(0,4):
#             line[num] = float(line[num])
#         time.append(line[0])
#         comp_N.append(line[1])
#         comp_E.append(line[2])
#         comp_Z.append(line[3])
# =============================================================================

#for row in num_lines:
#    line = row.split('\n')
#    line = ' '.join(line)
#    if not line.startswith("#"):
#        line = line.split('\t')
#        all_text == all_text + line

#time_ms = np.fromstring(time, dtype='float', sep=' ')
#N_ms = np.fromstring(comp_N, dtype='float', sep=' ')
#E_ms = np.fromstring(comp_E, dtype='float', sep=' ')
#Z_ms = np.fromstring(comp_Z, dtype='float', sep=' ')

st = Stream([Trace(data=mseed_data, header=stats)])
# write as ASCII file (encoding=0)
#st.write("data.mseed", format='MSEED', encoding=0, reclen=256) 
st.write("data.mseed", format='MSEED', reclen=256) 



   



#plt.title('Beverly Hills, Mag 4.2 - ' + date_title[0:19] + ' ' + '- ' + component + '-component acceleration')

#    plt.title('Event ' + str(select_event) + ' ' + date_title[0:19] + ' ' + ' ' + component + ' acc.')

#max_amp = 0
#for tr in st:
#    if (tr.stats.station not in bad_sta):
#        if (tr.stats.location not in bad_loc):
#            if tr.stats.station in st_names:
#                sta_index = st_names.index(tr.stats.station)
#                stalat = float(st_lats[sta_index])
#                stalon = float(st_lons[sta_index]) # look up lat & lon again to find distance
#                distance = gps2dist_azimuth(stalat,stalon,ref_lat,ref_lon) # Get traveltimes again, hard to store
#            max_local = np.max(np.abs(tr.data))*distance[0]
#            if max_local > max_amp:
#                max_amp = max_local
##print('Also made it to here!')
#for tr in st:
#    if (tr.stats.station not in bad_sta):
#        if (tr.stats.location not in bad_loc):
#            if tr.stats.station in st_names:
#                sta_index = st_names.index(tr.stats.station)
#                stalat = float(st_lats[sta_index])
#                stalon = float(st_lons[sta_index]) # look up lat & lon again to find distance
#                distance = gps2dist_azimuth(stalat,stalon,ref_lat,ref_lon) # Get traveltimes again, hard to store
#            else:
#                print(tr.stats.station + ' is not in station list')
#
#            dist_offset = distance[0]/1000 # correct to km
##                print(f'{tr.stats.station:s} {tr.stats.channel:s} distance is {dist_offset:.1f} slat is {stalat:.3f}  slon is {stalon:.3f}')
##                print(str(tr.stats.starttime) + '  ' + str(t))
#            ttt = tr.stats.starttime - t + np.arange(len(tr.data)) * tr.stats.delta
#
#            if tr.stats.station in basin_sta: # color basin stations red
#                if norm == 1:  # trace normalize amplitude
#                    plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
#                        - tr.data.min()) + dist_offset, color = 'red')
#                else:  # use common amplitude scale, multiply by distance
#                    plt.plot(ttt, (tr.data - np.median(tr.data))*(2*distance[0]*plot_scale_fac /max_amp)
#                     + dist_offset, color = 'red')
#                plt.text(s = tr.stats.station,x = end_buff*0.95,y = dist_offset
#                        + (max_plot_dist-min_plot_dist)*0.015, color = 'red')  #label traces
#            else:  # color stations outside basin black
#                if norm == 1:  # trace normalize amplitude
#                    plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
#                        - tr.data.min()) + dist_offset, color = 'black')
#                else:
#                    plt.plot(ttt, (tr.data - np.median(tr.data))*(2*distance[0]*plot_scale_fac /max_amp)
#                     + dist_offset, color = 'black')
#                plt.text(s = tr.stats.station,x = end_buff*0.95,y = dist_offset
#                    + (max_plot_dist-min_plot_dist)*0.015, color = 'black')  #label traces
#
#date_title = str(t)
#plt.title('Beverly Hills, Mag 4.2 - ' + date_title[0:19] + ' ' + '- ' + component + '-component acceleration')
##    plt.title('Event ' + str(select_event) + ' ' + date_title[0:19] + ' ' + ' ' + component + ' acc.')
#plt.xlabel('Time (s)')
#plt.ylabel('Distance (km)')
#
# 
## 
## fig = plt.figure()
## 
## ax1 = fig.add_subplot()
## ax2 = fig.add_subplot()
## ax3 = fig.add_subplot()
## 
## ax1.set_title("Velocity - North-South")    
## ax1.set_xlabel('Time (s)')
## ax1.set_ylabel('Velocity (m/s)')
## 
## ax1.plot(time[0:1000],comp_N[0:1000])
## ax2.plot(time[0:1000],comp_E[0:1000])
## ax3.plot(time[0:1000],comp_Z[0:1000])
## 
## leg = ax1.legend()
## 
## plt.show()
## =============================================================================
"""