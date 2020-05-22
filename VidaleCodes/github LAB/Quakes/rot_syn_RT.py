#!/usr/bin/env python
# program to rotate synthetics to radial and transverse for LA basin motions
# John Vidale 7/2019

def rot_syn_RT(event_no, start_buff, end_buff, taper, taper_frac,
			plot_scale_fac, filt, freq_min, freq_max, min_dist, max_dist,
			vmodel, basin_width, basin):

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

	show_data = 0
	start_time_wc = time.time()
	print('This is event number: ' + str(event_no) + ' vmod is ' + vmodel)

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
		st_dist.append(distance[0]/1000.) # azimuth
		st_az.append(distance[1]) # azimuth
		st_baz.append(distance[2]) # back-azimuth
	print('number of stations in list is ' + str(len(st_num)) + ' or ' + str(station_index))

	#%% Load data and synthetic waveforms
	st_dat = Stream()
	st_synE = Stream()
	st_synN = Stream()
	st_synZ = Stream()
	if vmodel == 'cvmhy':
		fname_dat     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_data/vel/ve_' + event_no + '.mseed'
		fname_synE     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/vel/cvmhy/ve_'  + event_no + '_cvmhy.mseed'
		fname_synN     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/vel/cvmhy/vn_'  + event_no + '_cvmhy.mseed'
		fname_synZ     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/vel/cvmhy/vz_'  + event_no + '_cvmhy.mseed'
	elif vmodel == 'cvms426-223':
		fname_dat     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_data/vel/ve_' + event_no + '.mseed'
		fname_synE     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/vel/cvms426-223/ve_'  + event_no + '_cvms426-223.mseed'
		fname_synN     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/vel/cvms426-223/vn_'  + event_no + '_cvms426-223.mseed'
		fname_synZ     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/vel/cvms426-223/vz_'  + event_no + '_cvms426-223.mseed'
	elif vmodel == 'cvmhn':
		fname_dat     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_data/vel/ve_' + event_no + '.mseed'
		fname_synE     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/vel/cvmhn/ve_'  + event_no + '_cvmhn.mseed'
		fname_synN     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/vel/cvmhn/vn_'  + event_no + '_cvmhn.mseed'
		fname_synZ     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/vel/cvmhn/vz_'  + event_no + '_cvmhn.mseed'
	elif vmodel == 'cvms400-100':
		fname_dat     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_data/vel/ve_' + event_no + '.mseed'
		fname_synE     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/vel/cvms400-100/ve_'  + event_no + '_cvms400-100.mseed'
		fname_synN     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/vel/cvms400-100/vn_'  + event_no + '_cvms400-100.mseed'
		fname_synZ     = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/vel/cvms400-100/vz_'  + event_no + '_cvms400-100.mseed'
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

	#%% select data by distance (and azimuth?), and cull synthetics to match data
	st_dat_select = Stream()
	st_synE_select = Stream()
	st_synN_select = Stream()
	st_synZ_select = Stream()
	for tr in st_dat: # examine traces one by one
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
#						tr.stats.distance = st_dist[ii] # add distance to trace metadata
#						st_datE_select += tr
						tr.stats.distance = st_dist[ii] # add distance to trace metadata
						st_synE[ii].stats.distance = st_dist[ii]
						st_synN[ii].stats.distance = st_dist[ii]
						st_synZ[ii].stats.distance = st_dist[ii]
						st_dat_select += tr
						st_synE_select += st_synE[ii]
						st_synN_select += st_synN[ii]
						st_synZ_select += st_synZ[ii]
#						print(tr.stats.station + ' counter ' + str(ii) + ' ' + 'num ' + str(int(st_num[ii])))
					else:
						print('Not in basin: '  + tr.stats.station)
				else:
					print('Too far: '  + tr.stats.station)
	print('now data has ' + str(len(st_dat_select)) + ' traces')
	print('now synH has '  + str(len(st_synE_select)) + ' traces')
	print('now synS has '  + str(len(st_synN_select)) + ' traces')
	print('now synS has '  + str(len(st_synZ_select)) + ' traces')

	#%%  reject data and sythetics on bad trace list, either individual components or A for all components
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
	for tr in st_synE_select: # examine traces one by one
		do_write = 1
		for ii in range(len(badt_event)):
			if event_no == badt_event[ii] and tr.stats.station == badt_station[ii]: # find station in inventory
				if badt_compo[ii] == 'E' or badt_compo[ii] == 'N' or badt_compo[ii] == 'A':
					do_write = 0
		if do_write == 1:
			st_synE_good += tr
	print('After rejecting labeled bad traces ones, synE has '       + str(len(st_synE_good))       + ' traces')

	st_synN_good = Stream()
	for tr in st_synN_select: # examine traces one by one
		do_write = 1
		for ii in range(len(badt_event)):
			if event_no == badt_event[ii] and tr.stats.station == badt_station[ii]: # find station in inventory
				if badt_compo[ii] == 'E' or badt_compo[ii] == 'N' or badt_compo[ii] == 'A':
					do_write = 0
		if do_write == 1:
			st_synN_good += tr
	print('After rejecting labeled bad traces ones, synN has '       + str(len(st_synN_good))       + ' traces')

	st_synZ_good = Stream()
	for tr in st_synZ_select: # examine traces one by one
		do_write = 1
		for ii in range(len(badt_event)):
			if event_no == badt_event[ii] and tr.stats.station == badt_station[ii]: # find station in inventory
				if badt_compo[ii] == 'Z' or badt_compo[ii] == 'A':
					do_write = 0
		if do_write == 1:
			st_synZ_good += tr
	print('After rejecting labeled bad traces ones, synZ has '       + str(len(st_synZ_good))       + ' traces')

#%%  rotate horizontals, so N becomes radial, E becomes transverse
	trN     = Trace()
	trE     = Trace()
	ntrace = len(st_synE_good)
	for iii in range(ntrace):  # loop over all good stations
#		print('In rotate loop, ntrace is ' + str(ntrace))
		for ii in range(len(st_name)):  # find matching entry in station roster
			if (st_synE_good[iii].stats.network == st_netw[ii] and st_synE_good[iii].stats.station == st_name[ii]): # find station in inventory
					distance = gps2dist_azimuth( ev_lat, ev_lon, float(st_lat[ii]), float(st_lon[ii])) # station 14825
					baz = distance[2]  # back-azimuth
		ba = np.radians(baz)
		trN = st_synN_good[iii].copy()
		trE = st_synE_good[iii].copy()
#		print(st_synE_good[iii].stats.station + ' ' + str(baz))
		st_synN_good[iii].data = - (trE.data * np.sin(ba)) - (trN.data * np.cos(ba))
		st_synE_good[iii].data = - (trE.data * np.cos(ba)) + (trN.data * np.sin(ba))

	pathRT = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo /Mseed_syn/rot/' + vmodel + '/v'
	st_synE_good.write(pathRT + 't_' + event_no + '_' + vmodel + '.mseed', format='MSEED')
	st_synN_good.write(pathRT + 'r_' + event_no + '_' + vmodel + '.mseed', format='MSEED')
	st_synZ_good.write(pathRT + 'z_' + event_no + '_' + vmodel + '.mseed', format='MSEED')

	#%%  detrend, taper, filter
	if taper:
		st_dat_good.detrend( type='simple')
		st_synE_good.detrend(type='simple')
		st_synN_good.detrend(type='simple')
		st_synZ_good.detrend(type='simple')
	if filt:
		st_dat_good.filter( 'bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
		st_synE_good.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
		st_synN_good.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
		st_synZ_good.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
	if taper:
		st_dat_good.detrend( type='simple')
		st_synE_good.detrend(type='simple')
		st_synN_good.detrend(type='simple')
		st_synZ_good.detrend(type='simple')

	#%%
	# plot traces
	if vmodel == 'cvms426-223':
		fig_index = 10
	elif vmodel == 'cvmhy':
		fig_index = 11
	elif vmodel == 'cvms400-100':
		fig_index = 12
	elif vmodel == 'cvmhn':
		fig_index = 13
	plt.close(fig_index)
	plt.figure(fig_index,figsize=(10,8))
	plt.xlim(-start_buff,end_buff)
	plt.ylim(min_dist,max_dist)

	# find max
	maxE = 0
	for tr in st_synE_good:
		tr_max = max(abs(tr.data))*tr.stats.distance
		if tr_max > maxE:
			maxE = tr_max
	maxN = 0
	for tr in st_synN_good:
		tr_max = max(abs(tr.data))*tr.stats.distance
		if tr_max > maxN:
			maxN = tr_max
	maxZ = 0
	for tr in st_synZ_good:
		tr_max = max(abs(tr.data))*tr.stats.distance
		if tr_max > maxZ:
			maxZ = tr_max
	print('Max E, N, and Z synthetic are ' + str(maxE) + '  ' + str(maxN) + '  ' + str(maxZ))

	max_all = max(maxZ, maxN, maxE)
	plot_fac = plot_scale_fac * (max_dist - min_dist) / max_all

	if show_data:
		for tr in st_dat_good:
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
	for tr in st_dat_good:
		dist_offset = tr.stats.distance # km
		plt.text(s = tr.stats.network + ' ' + tr.stats.station ,x = end_buff*0.95,y = dist_offset
			    + max_dist*0.015, color = 'black')  #label traces
#
	for tr in st_synE_good:
		dist_offset = tr.stats.distance
		ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)#	These lines used to cause a crash in Spyder
		plt.plot(ttt, (tr.data * plot_fac * tr.stats.distance) + dist_offset, color = 'black')

	for tr in st_dat_good:
		dist_offset = tr.stats.distance # km
		plt.text(s = tr.stats.network + ' ' + tr.stats.station ,x = end_buff*0.95,y = dist_offset
			    + max_dist*0.015, color = 'black')  #label traces

	for tr in st_synN_good:
		dist_offset = tr.stats.distance
		ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)#	These lines used to cause a crash in Spyder
		plt.plot(ttt, (tr.data * plot_fac * tr.stats.distance) + dist_offset, color = 'green')

	for tr in st_synZ_good:
		dist_offset = tr.stats.distance
		ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)#	These lines used to cause a crash in Spyder
		plt.plot(ttt, (tr.data * plot_fac * tr.stats.distance) + dist_offset, color = 'red')

	plt.xlabel('Time (s)')
	plt.ylabel('Epicentral distance from event (km)')
#	plt.title(fname1[8:18] + ' vs ' + fname2[8:18])
	plt.title(date_label + ' ' + event_no + ' black E green N red Z of vmodel ' + vmodel)
	plt.show()

	#  Save processed files
#	fname1 = 'Pro_Files/HD' + date_label1 + 'sel.mseed'
#	fname2 = 'Pro_Files/HD' + date_label2 + 'sel.mseed'
#	st1good.write(fname1,format = 'MSEED')
#	st2good.write(fname2,format = 'MSEED')

	elapsed_time_wc = time.time() - start_time_wc
	print('This job took ' + str(elapsed_time_wc) + ' seconds')
	#os.system('say "Done"')