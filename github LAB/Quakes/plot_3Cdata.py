#!/usr/bin/env python
# program to compare data and synthetics for LA basin motions
# John Vidale 6/2019

def p_3Cdata(event_no, start_buff, end_buff, taper, taper_frac, norm_each,
			plot_scale_fac, filt, freq_min, freq_max, min_dist, max_dist,
			 motion, basin_width, basin, rotated, CI_only):

	from obspy import UTCDateTime
	from obspy import Stream, Trace
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

	#%% find event details for origin time, lat, lon
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
#		print('Event ' + str(ev_lat) + ' ' + str(ev_lon) + ' station ' + split_line[4] + ' ' + split_line[5] + ' distance ' + str(distance[0]))
		st_dist.append(distance[0]/1000.) # azimuth
		st_az.append(distance[1]) # azimuth
		st_baz.append(distance[2]) # back-azimuth
	print('number of stations in list is ' + str(len(st_num)) + ' or ' + str(station_index))
#	print('lat is  ' + str(st_lat[807]) + ',  lon is ' + str(st_lon[807]) + ' distance is ' + str(st_dist[807]))

	#%% Load data and synthetic waveforms
	st_datE = Stream()
	st_datN = Stream()
	st_datZ = Stream()
	nstart = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_data/'
	if rotated:
		nstart = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_data/rot/v'
		fname_datE = nstart + 'r_' + event_no + '.mseed'
		fname_datN = nstart + 't_' + event_no + '.mseed'
		fname_datZ = nstart + 'z_' + event_no + '.mseed'
	elif motion == 'acc':
		fname_datE = nstart + 'acc/ve_' + event_no + '.mseed'
		fname_datN = nstart + 'acc/vn_' + event_no + '.mseed'
		fname_datZ = nstart + 'acc/vz_' + event_no + '.mseed'
	elif motion == 'vel':
		fname_datE = nstart + 'vel/ve_' + event_no + '.mseed'
		fname_datN = nstart + 'vel/vn_' + event_no + '.mseed'
		fname_datZ = nstart + 'vel/vz_' + event_no + '.mseed'
	elif motion == 'disp':
		fname_datE = nstart + 'disp/ve_' + event_no + '.mseed'
		fname_datN = nstart + 'disp/vn_' + event_no + '.mseed'
		fname_datZ = nstart + 'disp/vz_' + event_no + '.mseed'
	else:
		print('error: trying to write unknown component of motion')
	st_datE=read(fname_datE)
	st_datN=read(fname_datN)
	st_datZ=read(fname_datZ)

	print('1st data trace has : ' + str(len(st_datE[0].data)) + ' time pts ')
	print('E has ' + str(len(st_datE)) + ' traces')
	print('N has ' + str(len(st_datN)) + ' traces')
	print('Z has ' + str(len(st_datZ)) + ' traces')
	print('E trace starts at ' + str(st_datE[0].stats.starttime) + ', event at ' + str(t1))
	print('N trace starts at ' + str(st_datN[0].stats.starttime) + ', event at ' + str(t1))
	print('Z trace starts at ' + str(st_datZ[0].stats.starttime) + ', event at ' + str(t1))

		#%%  reject data on bad trace list, either individual components or A for all components
	st_datE_good = Stream()
	for tr in st_datE: # examine traces one by one
		do_write = 1
		for ii in range(len(badt_event)):
			if event_no == badt_event[ii] and tr.stats.station == badt_station[ii]: # find station in inventory
				if badt_compo[ii] == 'E' or badt_compo[ii] == 'A':
					do_write = 0
		if do_write == 1:
			st_datE_good += tr
	print('After rejecting labeled bad traces ones, E data has '       + str(len(st_datE_good))       + ' traces')

	st_datN_good = Stream()
	for tr in st_datN: # examine traces one by one
		do_write = 1
		for ii in range(len(badt_event)):
			if event_no == badt_event[ii] and tr.stats.station == badt_station[ii]: # find station in inventory
				if badt_compo[ii] == 'N' or badt_compo[ii] == 'A':
					do_write = 0
		if do_write == 1:
			st_datN_good += tr
	print('After rejecting labeled bad traces ones, N data has '       + str(len(st_datN_good))       + ' traces')

	st_datZ_good = Stream()
	for tr in st_datZ: # examine traces one by one
		do_write = 1
		for ii in range(len(badt_event)):
			if event_no == badt_event[ii] and tr.stats.station == badt_station[ii]: # find station in inventory
				if badt_compo[ii] == 'Z' or badt_compo[ii] == 'A':
					do_write = 0
		if do_write == 1:
			st_datZ_good += tr
	print('After rejecting labeled bad traces ones, Z data has '       + str(len(st_datZ_good))       + ' traces')

#%%
	# select data by distance (and azimuth?), and cull synthetics to match data
	st_datE_select = Stream()
	st_datN_select = Stream()
	st_datZ_select = Stream()
	for tr in st_datE_good: # examine traces one by one
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
	#					print(tr.stats.station + ' ' + str(dist1) + ' ' + str(dist2) + ' ' + str(dist3))
						if basin == False or (dist1 < basin_width) or (dist2 < basin_width) or (dist3 < basin_width):  # keep stations only within X km of basin axis
	#						print('selected: ' + tr.stats.station)
							tr.stats.distance = st_dist[ii] # add distance to trace metadata
							st_datE_select += tr
	for tr in st_datN_good: # examine traces one by one
		if tr.stats.network == 'CI' or CI_only == False:
			for ii in range(len(st_name)):  # find matching entry in station roster
				if (tr.stats.network == st_netw[ii] and tr.stats.station == st_name[ii]): # find station in inventory
					if st_dist[ii] < max_dist: # exclude stations too far away
						# basin is roughly within 15 km of these 3 stations
						distance = gps2dist_azimuth( 33.99053, -118.36171, float(st_lat[ii]), float(st_lon[ii])) # station BHP
						dist1 = distance[0]/1000  # convert m to km
						distance = gps2dist_azimuth( 33.88110, -118.17568, float(st_lat[ii]), float(st_lon[ii])) # station LTP
						dist2 = distance[0]/1000
						distance = gps2dist_azimuth( 33.80776, -117.98116, float(st_lat[ii]), float(st_lon[ii])) # station BRE
						dist3 = distance[0]/1000
	#					print(tr.stats.station + ' ' + str(dist1) + ' ' + str(dist2) + ' ' + str(dist3))
						if basin == False or (dist1 < basin_width) or (dist2 < basin_width) or (dist3 < basin_width):  # keep stations only within X km of basin axis
	#						print('selected: ' + tr.stats.station)
							tr.stats.distance = st_dist[ii]
							st_datN_select += tr
	for tr in st_datZ_good: # examine traces one by one
		if tr.stats.network == 'CI' or CI_only == False:
			for ii in range(len(st_name)):  # find matching entry in station roster
				if (tr.stats.network == st_netw[ii] and tr.stats.station == st_name[ii]): # find station in inventory
					if st_dist[ii] < max_dist: # exclude stations too far away
						# basin is roughly within 15 km of these 3 stations
						distance = gps2dist_azimuth( 33.99053, -118.36171, float(st_lat[ii]), float(st_lon[ii])) # station BHP
						dist1 = distance[0]/1000  # convert m to km
						distance = gps2dist_azimuth( 33.88110, -118.17568, float(st_lat[ii]), float(st_lon[ii])) # station LTP
						dist2 = distance[0]/1000
						distance = gps2dist_azimuth( 33.80776, -117.98116, float(st_lat[ii]), float(st_lon[ii])) # station BRE
						dist3 = distance[0]/1000
	#					print(tr.stats.station + ' ' + str(dist1) + ' ' + str(dist2) + ' ' + str(dist3))
						if basin == False or (dist1 < basin_width) or (dist2 < basin_width) or (dist3 < basin_width):  # keep stations only within X km of basin axis
	#						print('selected: ' + tr.stats.station)
							tr.stats.distance = st_dist[ii]
							st_datZ_select += tr
	print('now E has ' + str(len(st_datE_select)) + ' traces')
	print('now N has ' + str(len(st_datN_select)) + ' traces')
	print('now Z has ' + str(len(st_datZ_select)) + ' traces')

	#%%  detrend, taper, filter
	if taper:
#		st_datE_select.detrend(type='simple')
#		st_datN_select.detrend(type='simple')
#		st_datZ_select.detrend(type='simple')
		st_datE_select.taper(taper_frac)
		st_datN_select.taper(taper_frac)
		st_datZ_select.taper(taper_frac)
	if filt:
		st_datE_select.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
		st_datN_select.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
		st_datZ_select.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
	if taper:
		st_datE_select.taper(taper_frac)
		st_datN_select.taper(taper_frac)
		st_datZ_select.taper(taper_frac)
#		st_datE_select.detrend(type='simple')
#		st_datN_select.detrend(type='simple')
#		st_datZ_select.detrend(type='simple')

	#%%
	# plot traces
	fig_index = 7
	plt.close(fig_index)
	plt.figure(fig_index,figsize=(10,8))
	plt.xlim(-start_buff,end_buff)
	plt.ylim(min_dist,max_dist)

	# find maxs
	maxE = 0
	for tr in st_datE_select:
		tr_max = max(abs(tr.data))*tr.stats.distance
		if tr_max > maxE:
			maxE = tr_max
	maxN = 0
	for tr in st_datN_select:
		tr_max = max(abs(tr.data))*tr.stats.distance
		if tr_max > maxN:
			maxN = tr_max
	maxZ = 0
	for tr in st_datZ_select:
		tr_max = max(abs(tr.data))*tr.stats.distance
		if tr_max > maxZ:
			maxZ = tr_max
	max_all = max(maxZ, maxN, maxE)
#	print(f'Normed amp of data and H, S, H_simp, and S_simp synthetics are {maxD_N:6.4f} {maxSH_N:6.4f} {maxSS_N:6.4f} {maxSH_N_simp:6.4f} + {maxSS_N_simp:6.4f}')
	print('Max E N Z all are ' + str(maxE) + '  ' + str(maxN) + '  ' + str(maxZ) + '  ' + str(max_all))

	plot_fac = plot_scale_fac * (max_dist - min_dist) / max_all

	for tr in st_datE_select:
		dist_offset = tr.stats.distance # km
		ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)
		if norm_each:
			plt.plot(ttt, (tr.data - np.median(tr.data))*plot_fac /(tr.data.max()
			- tr.data.min()) + dist_offset, color = 'black')
		else:
			plt.plot(ttt, (tr.data * plot_fac * tr.stats.distance) + dist_offset, color = 'black')

	for tr in st_datN_select:
		dist_offset = tr.stats.distance # km
		ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)
		if norm_each:
			plt.plot(ttt, (tr.data - np.median(tr.data))*plot_fac /(tr.data.max()
				- tr.data.min()) + dist_offset, color = 'green')
		else:
			plt.plot(ttt, (tr.data * plot_fac * tr.stats.distance) + dist_offset, color = 'green')

	for tr in st_datZ_select:
		dist_offset = tr.stats.distance # km
		ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)
		if norm_each:
			plt.plot(ttt, (tr.data - np.median(tr.data))*plot_fac /(tr.data.max()
				- tr.data.min()) + dist_offset, color = 'red')
		else:
			plt.plot(ttt, (tr.data * plot_fac * tr.stats.distance) + dist_offset, color = 'red')

	for tr in st_datE_select:
		dist_offset = tr.stats.distance # km
		plt.text(s = tr.stats.network + ' ' + tr.stats.channel + ' ' + tr.stats.station,x = end_buff*0.95,y = dist_offset
			    + max_dist*0.015, color = 'black')  #label traces

	plt.xlabel('Time (s)')
	plt.ylabel('Epicentral distance from event (km)')
#	plt.title(fname1[8:18] + ' vs ' + fname2[8:18])
	if rotated:
		if motion == 'acc':
			plt.title('Data ' + date_label + '  ' + event_no + ' 3C Acc  Radial-black Trans-green Z-red')
		elif motion == 'vel':
			plt.title('Data ' + date_label + '  ' + event_no + ' 3C Vel  Radial-black Trans-green Z-red')
		elif motion == 'disp':
			plt.title('Data ' + date_label + '  ' + event_no + ' 3C Disp  Radial-black Trans-green Z-red')
		else:
			print('error: trying to write unknown component of motion')
	else:
		if motion == 'acc':
			plt.title('Data ' + date_label + '  ' + event_no + ' 3C Acc  E-black N-green Z-red')
		elif motion == 'vel':
			plt.title('Data ' + date_label + '  ' + event_no + ' 3C Vel  E-black N-green Z-red')
		elif motion == 'disp':
			plt.title('Data ' + date_label + '  ' + event_no + ' 3C Disp  E-black N-green Z-red')
		else:
			print('error: trying to write unknown component of motion')
	plt.show()

	#  Save processed files
#	fname1 = 'Pro_Files/HD' + date_label1 + 'sel.mseed'
#	fname2 = 'Pro_Files/HD' + date_label2 + 'sel.mseed'
#	st1good.write(fname1,format = 'MSEED')
#	st2good.write(fname2,format = 'MSEED')

	elapsed_time_wc = time.time() - start_time_wc
	print('This job took ' + str(elapsed_time_wc) + ' seconds')
	#os.system('say "Done"')