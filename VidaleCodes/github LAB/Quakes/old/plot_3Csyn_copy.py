#!/usr/bin/env python
# program to plot synthetics for LA basin motions
#  seems to have option to also plot E component of data, perhaps needed to know which stations have data
# John Vidale 6/2019

def p_3Csyn_copy(event_no, start_buff, end_buff, taper, taper_frac,
			plot_scale_fac, filt, freq_min, freq_max, min_dist, max_dist,
			vmodel):

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

	show_data = 0
	start_time_wc = time.time()

	#%% find event details for origin time, lat, lon
	ev_file = '/Users/vidale/Documents/PyCode/LAB/ricardo/event_list.txt'
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

	#%% Open station location file
	sta_file = '/Users/vidale/Documents/PyCode/LAB/ricardo/ricardo_stations.txt'
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
	if vmodel == 'H':
		fname_dat     = '/Users/vidale/Documents/PyCode/LAB/ricardo/Mseed_data/vel/ve_' + event_no + '.mseed'
		fname_synE     = '/Users/vidale/Documents/PyCode/LAB/ricardo/Mseed_syn/vel/cvms400-100/ve_'  + event_no + '_cvms400-100.mseed'
		fname_synN     = '/Users/vidale/Documents/PyCode/LAB/ricardo/Mseed_syn/vel/cvms400-100/vn_'  + event_no + '_cvms400-100.mseed'
		fname_synZ     = '/Users/vidale/Documents/PyCode/LAB/ricardo/Mseed_syn/vel/cvms400-100/vz_'  + event_no + '_cvms400-100.mseed'
	elif vmodel == 'S':
		fname_dat     = '/Users/vidale/Documents/PyCode/LAB/ricardo/Mseed_data/vel/ve_' + event_no + '.mseed'
		fname_synE     = '/Users/vidale/Documents/PyCode/LAB/ricardo/Mseed_syn/vel/cvms426-223/ve_'  + event_no + '_cvms426-223.mseed'
		fname_synN     = '/Users/vidale/Documents/PyCode/LAB/ricardo/Mseed_syn/vel/cvms426-223/vn_'  + event_no + '_cvms426-223.mseed'
		fname_synZ     = '/Users/vidale/Documents/PyCode/LAB/ricardo/Mseed_syn/vel/cvms426-223/vz_'  + event_no + '_cvms426-223.mseed'
	st_dat=read(fname_dat)
	st_synE=read(fname_synE)
	st_synN=read(fname_synN)
	st_synZ=read(fname_synZ)
	print('In sgram file ' + str(st_synZ[609].data[0]) + '  ' + st_synZ[609].stats.station + '  ' + str(len(st_synZ)))
	print('In arrays st_name[609] ' + st_name[609] + ' st_name[610] ' + st_name[610])
	print('In arrays st_name[0] ' + st_name[0] + ' st_name[1] ' + st_name[1])

	print('1st data trace has : ' + str(len(st_synE[0].data)) + ' time pts ')
	print('synE has ' + str(len(st_synE)) + ' traces')
	print('synN has ' + str(len(st_synN)) + ' traces')
	print('synZ has ' + str(len(st_synZ)) + ' traces')

	#%%
	# select data by distance (and azimuth?), and cull synthetics to match data
	st_dat_select = Stream()
	st_synE_select = Stream()
	st_synN_select = Stream()
	st_synZ_select = Stream()
	for tr in st_dat: # examine traces one by one
		for ii in range(len(st_name)):  # find matching entry in station roster
			if (tr.stats.network == st_netw[ii] and tr.stats.station == st_name[ii]): # find station in inventory
#				if tr.stats.station == 'LKL':
#					print('counter ' + str(ii) + ' name at counter ' + st_name[ii] +
#					   ' name at counter-1 ' + st_name[ii-1] + ' ' + 'num ' + str(int(st_num[ii])) + ' 609 ' +
#						 str(st_synZ[609].data[0]) + ' 610 ' + str(st_synZ[610].data[0]))
#					print('counter ' + str(ii) + ' dist at counter ' + str(st_dist[ii]) +
#					   ' dist at counter-1 ' + str(st_dist[ii-1]) + ' lat at counter ' +
#					   str(st_lat[ii]) + ' lat at counter-1 ' + str(st_lat[ii-1]))
				if (st_dist[ii] < max_dist) and (st_dist[ii] > min_dist): # exclude stations too close or too far
					tr.stats.distance = st_dist[ii] # add distance to trace metadata
					# I don't understand the -1 in the following 7 lines, but it works
					st_synE[ii-1].stats.distance = st_dist[ii]
					st_synN[ii-1].stats.distance = st_dist[ii]
					st_synZ[ii-1].stats.distance = st_dist[ii]
					st_dat_select += tr
					st_synE_select += st_synE[ii-1]
					st_synN_select += st_synN[ii-1]
					st_synZ_select += st_synZ[ii-1]
					print(tr.stats.station + ' counter ' + str(ii) + ' ' + 'num ' + str(int(st_num[ii])))
				else:
					print('Too far: '  + tr.stats.station)
	print('now data has ' + str(len(st_dat_select)) + ' traces')
	print('now synH has '  + str(len(st_synE_select)) + ' traces')
	print('now synS has '  + str(len(st_synN_select)) + ' traces')
	print('now synS has '  + str(len(st_synZ_select)) + ' traces')

	#%%  detrend, taper, filter
	if taper != 0:
		st_dat_select.detrend(type='simple')
		st_synE_select.detrend(type='simple')
		st_synN_select.detrend(type='simple')
		st_synZ_select.detrend(type='simple')
	if filt != 0:
		st_dat_select.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=2, zerophase=True)
		st_synE_select.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=2, zerophase=True)
		st_synN_select.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=2, zerophase=True)
		st_synZ_select.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=2, zerophase=True)
	if taper != 0:
		st_dat_select.detrend(type='simple')
		st_synE_select.detrend(type='simple')
		st_synN_select.detrend(type='simple')
		st_synZ_select.detrend(type='simple')

	#%%
	# plot traces
	if vmodel == 'S':
		fig_index = 8
	elif vmodel == 'H':
		fig_index = 10
	plt.close(fig_index)
	plt.figure(fig_index,figsize=(20,10))
	plt.xlim(-start_buff,end_buff)
	plt.ylim(min_dist,max_dist)

	# find max
	maxE = 0
	for tr in st_synE_select:
		tr_max = max(abs(tr.data))
		if tr_max > maxE:
			maxE = tr_max
	maxN = 0
	for tr in st_synN_select:
		tr_max = max(abs(tr.data))
		if tr_max > maxN:
			maxN = tr_max
	maxZ = 0
	for tr in st_synZ_select:
		tr_max = max(abs(tr.data))
		if tr_max > maxZ:
			maxZ = tr_max
	print('Max E, N, and Z synthetic are ' + str(maxE) + '  ' + str(maxN) + '  ' + str(maxZ))

	if show_data:
		for tr in st_dat_select:
			dist_offset = tr.stats.distance # km
			ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)
	#		if red_plot == 1:
	#			shift = red_time + (dist_offset - red_dist) * red_slow
	#			ttt = ttt - shift
	#		plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
	#			- tr.data.min()) + dist_offset, color = 'green')
			plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
				- tr.data.min()) + dist_offset, color = 'black')
			print(str(tr.stats.distance) + ' distance ' + tr.stats.station + ' station')	#plt.title(fname1)

	#print labels whether or not data is shown
	for tr in st_dat_select:
		dist_offset = tr.stats.distance # km
		plt.text(s = tr.stats.network + ' ' + tr.stats.station ,x = end_buff*0.95,y = dist_offset
			    + max_dist*0.015, color = 'black')  #label traces
#
	for tr in st_synE_select:
		dist_offset = tr.stats.distance
		ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)#	These lines used to cause a crash in Spyder
		plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
			- tr.data.min()) + dist_offset, color = 'black')

	for tr in st_dat_select:
		dist_offset = tr.stats.distance # km
		plt.text(s = tr.stats.network + ' ' + tr.stats.station ,x = end_buff*0.95,y = dist_offset
			    + max_dist*0.015, color = 'black')  #label traces

	for tr in st_synN_select:
		dist_offset = tr.stats.distance
		ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)#	These lines used to cause a crash in Spyder
		plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
			- tr.data.min()) + dist_offset, color = 'green')

	for tr in st_synZ_select:
		dist_offset = tr.stats.distance
		ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)#	These lines used to cause a crash in Spyder
		plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
			- tr.data.min()) + dist_offset, color = 'red')

	plt.xlabel('Time (s)')
	plt.ylabel('Epicentral distance from event (km)')
#	plt.title(fname1[8:18] + ' vs ' + fname2[8:18])
	plt.title(date_label + ' ' + event_no + ' 3 components of vmodel ' + vmodel)
	plt.show()

	#  Save processed files
#	fname1 = 'Pro_Files/HD' + date_label1 + 'sel.mseed'
#	fname2 = 'Pro_Files/HD' + date_label2 + 'sel.mseed'
#	st1good.write(fname1,format = 'MSEED')
#	st2good.write(fname2,format = 'MSEED')

	elapsed_time_wc = time.time() - start_time_wc
	print('This job took ' + str(elapsed_time_wc) + ' seconds')
	os.system('say "Done"')