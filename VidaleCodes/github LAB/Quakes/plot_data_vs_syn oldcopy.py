#!/usr/bin/env python
# program to compare data and synthetics for LA basin motions
# John Vidale 6/2019

def p_data_vs_syn(event_no, start_buff, end_buff, taper, taper_frac,
			plot_scale_fac, filt, freq_min, freq_max, min_dist, max_dist,
			show_data, show_H, show_S, show_H_simp, show_S_simp,
			compo, norm_each, motion, basin_width, basin,
			com_time):
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

	#%% find event details for origin time, lat, lon
	ev_file = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/event_list.txt'
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
		st_dist.append(distance[0]/1000.) # distance
		st_az.append(distance[1]) # azimuth
		st_baz.append(distance[2]) # back-azimuth
	print('number of stations in list is ' + str(len(st_num)) + ' or ' + str(station_index))

	#%% Load data and all synthetic waveforms for this event
	st_dat       = Stream()
	st_synH      = Stream()
	st_synS      = Stream()
	st_synH_simp = Stream()
	st_synS_simp = Stream()

	# choose E, N, or Z component
	if compo == 'E':
		comp_pre = 've'
	elif compo == 'N':
		comp_pre = 'vn'
	elif compo == 'Z':
		comp_pre = 'vz'
	else:
		print(compo + ' is not a valid component')

	# choose acceleration, velocity, or displacement
	if motion == 'disp':
		motion_pre = 'disp'
	elif motion == 'vel':
		motion_pre = 'vel'
	elif motion == 'acc':
		motion_pre = 'acc'
	else:
		print(motion + ' is not a valid type of motion')

	fname_dat       = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_data/'            + motion_pre + '/' + comp_pre + '_' + event_no + '.mseed'
	fname_synH      = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/' + motion_pre + '/cvmhy/'       + comp_pre + '_' + event_no + '_cvmhy.mseed'
	fname_synS      = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/' + motion_pre + '/cvms426-223/' + comp_pre + '_' + event_no + '_cvms426-223.mseed'
	fname_synH_simp = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/' + motion_pre + '/cvmhn/'       + comp_pre + '_' + event_no + '_cvmhn.mseed'
	fname_synS_simp = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/' + motion_pre + '/cvms400-100/' + comp_pre + '_' + event_no + '_cvms400-100.mseed'

	st_dat       = read(fname_dat)
	st_synH      = read(fname_synH)
	st_synS      = read(fname_synS)
	st_synH_simp = read(fname_synH_simp)
	st_synS_simp = read(fname_synS_simp)

	lenH = len(st_synH)     # correct units from Hercules from meters to cm
	for ii in range(lenH):  # assume trace count of all synthetic files is same
		st_synH[ii].data      = 100 * st_synH[ii].data
		st_synS[ii].data      = 100 * st_synS[ii].data
		st_synH_simp[ii].data = 100 * st_synH_simp[ii].data
		st_synS_simp[ii].data = 100 * st_synS_simp[ii].data

	print('1st data trace has : ' + str(len(st_dat[0].data)) + ' time pts ')
	print('data has ' + str(len(st_dat)) + ' traces')
	print('synthetic gathers have ' + str(lenH) + ' traces')
	print('data trace starts at ' + str(st_dat[0].stats.starttime) + ', event at ' + str(t1))
	print('synthetic traces start at ' + str(st_synH[0].stats.starttime) + ', event at ' + str(t1))

	#%%  reject data on bad trace list, either individual components or A for all components
	st_dat_good = Stream()
	for tr in st_dat: # examine traces one by one
		do_write = 1
		for ii in range(len(badt_event)):
			if event_no == badt_event[ii] and tr.stats.station == badt_station[ii]: # find station in inventory
				if badt_compo[ii] == compo or badt_compo[ii] == 'A':
					do_write = 0
			if tr.stats.station == '24400':
				print(event_no + ' ' + badt_event[ii] + ' ' + tr.stats.station + ' ' + badt_station[ii] + ' ' + compo + ' ' + badt_compo[ii])
		if do_write == 1:
			st_dat_good += tr
	print('After rejecting labeled bad traces ones, data has '       + str(len(st_dat_good))       + ' traces')

	#%%  select data by distance (and azimuth?), and cull synthetics to match data
	st_dat_select       = Stream()
	st_synH_select      = Stream()
	st_synS_select      = Stream()
	st_synH_simp_select = Stream()
	st_synS_simp_select = Stream()
	for tr in st_dat_good: # examine traces one by one
		for ii in range(len(st_name)):  # find matching entry in station roster
			if tr.stats.network == st_netw[ii] and tr.stats.station == st_name[ii]: # find station in inventory
				if (st_dist[ii] < max_dist) and (st_dist[ii] > min_dist): # exclude stations too close or too far
					distance = gps2dist_azimuth( 33.976, -118.209, float(st_lat[ii]), float(st_lon[ii])) # station 14825
					dist1 = distance[0]/1000  # convert m to km
					distance = gps2dist_azimuth( 33.785, -117.890, float(st_lat[ii]), float(st_lon[ii])) # station 13893
					dist2 = distance[0]/1000
					distance = gps2dist_azimuth( 33.874, -118.039, float(st_lat[ii]), float(st_lon[ii])) # station CRF
					dist3 = distance[0]/1000
					if com_time: # trim data and syn traces to common time interval
						s_t = tr.stats.starttime
						e_t = tr.stats.endtime
						# all syns have common time interval
						if st_synH[ii].stats.starttime > s_t:
							s_t = st_synH[ii].stats.starttime
						if st_synH[ii].stats.endtime   < e_t:
							e_t = st_synH[ii].stats.endtime
						tr.trim(starttime=s_t,endtime = e_t)
						st_synH[ii].trim(starttime=s_t,endtime = e_t)
						st_synS[ii].trim(starttime=s_t,endtime = e_t)
						st_synH_simp[ii].trim(starttime=s_t,endtime = e_t)
						st_synS_simp[ii].trim(starttime=s_t,endtime = e_t)
#					print(tr.stats.station + ' ' + str(dist1) + ' ' + str(dist2) + ' ' + str(dist3))
					if basin == False or (dist1 < basin_width) or (dist2 < basin_width) or (dist3 < basin_width):  # keep stations only within X km of basin axis
#						print('selected: ' + tr.stats.station)
						tr.stats.distance = st_dist[ii] # add distance to trace metadata
						st_synH[ii].stats.distance      = st_dist[ii]
						st_synS[ii].stats.distance      = st_dist[ii]
						st_synH_simp[ii].stats.distance = st_dist[ii]
						st_synS_simp[ii].stats.distance = st_dist[ii]
						st_dat_select += tr
						st_synH_select      += st_synH[ii]
						st_synS_select      += st_synS[ii]
						st_synH_simp_select += st_synH_simp[ii]
						st_synS_simp_select += st_synS_simp[ii]
	print('now data has '       + str(len(st_dat_select))       + ' traces')
	print('now synH has '       + str(len(st_synH_select))      + ' traces')
	print('now synS has '       + str(len(st_synS_select))      + ' traces')
	print('now synH_simp has '  + str(len(st_synH_simp_select)) + ' traces')
	print('now synS_simp has '  + str(len(st_synS_simp_select)) + ' traces')


	#%%  detrend, taper, filter
	if taper:
		st_dat_select.taper(      taper_frac)
		st_synH_select.taper(     taper_frac)
		st_synS_select.taper(     taper_frac)
		st_synH_simp_select.taper(taper_frac)
		st_synS_simp_select.taper(taper_frac)
	if filt:
		st_dat_select.filter(      'bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
		st_synH_select.filter(     'bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
		st_synS_select.filter(     'bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
		st_synH_simp_select.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
		st_synS_simp_select.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
	if taper:
		st_dat_select.taper(      taper_frac)
		st_synH_select.taper(     taper_frac)
		st_synS_select.taper(     taper_frac)
		st_synH_simp_select.taper(taper_frac)
		st_synS_simp_select.taper(taper_frac)

	#%% plot traces
	if compo == 'E':
		fig_index = 13
	elif compo == 'N':
		fig_index = 14
	elif compo == 'Z':
		fig_index = 15
	plt.close(fig_index)
	plt.figure(fig_index,figsize=(10,8))
	plt.xlim(-start_buff,end_buff)
	plt.ylim(min_dist,max_dist)

	# find max of absolute amplitude
	maxD = 0
	for tr in st_dat_select:
		tr_max = max(abs(tr.data))
		if tr_max > maxD:
			maxD = tr_max
#	print('Max data and H, S, H_simp, and S_simp synthetics are ' + str(maxD) + '  ' + str(maxSH) + '  ' + str(maxSS)+ '  ' + str(maxSH_simp)+ '  ' + str(maxSS_simp))
#	print(f'Amp of data and H, S, H_simp, and S_simp synthetics are {maxD:6.4f} {maxSH:6.4f} {maxSS:6.4f} {maxSH_simp:6.4f} + {maxSS_simp:6.4f}')

	# find max normalized to 10 km distance (i.e., amp divided by distance, assumes 1/R amp fall-off)
	maxD_N = 0
	for tr in st_dat_select:
		tr_max = max(abs(tr.data))*tr.stats.distance
		if tr_max > maxD_N:
			maxD_N = tr_max
	print(f'Normed amp of data is {maxD_N:6.4f}')
	plot_fac = plot_scale_fac * (max_dist - min_dist) / maxD_N

	if show_data:
		for tr in st_dat_select:
			dist_offset = tr.stats.distance # km
			ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)
			if norm_each:
				plt.plot(ttt, (tr.data - np.median(tr.data))*plot_fac /(tr.data.max()
				- tr.data.min()) + dist_offset, color = 'black')
			else:
				plt.plot(ttt, (tr.data * plot_fac * tr.stats.distance) + dist_offset, color = 'black')
#			print(str(tr.stats.distance) + ' distance ' + tr.stats.station + ' station')	#plt.title(fname1)

	for tr in st_dat_select:
		dist_offset = tr.stats.distance # km
		pnum = max(abs(tr.data))
		printee2 = f'amp {pnum:6.4f}  {tr.stats.network}  {tr.stats.station}'
		plt.text(s = printee2,x = end_buff*0.8,y = dist_offset
			    + (max_dist - min_dist)*0.02, color = 'black')  #label traces

	if show_H:
		for tr in st_synH_select:
			dist_offset = tr.stats.distance
			ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)#	These lines used to cause a crash in Spyder
			if norm_each:
				plt.plot(ttt, (tr.data - np.median(tr.data))*plot_fac /(tr.data.max()
					- tr.data.min()) + dist_offset, color = 'red')
			else:
				plt.plot(ttt, (tr.data * plot_fac * tr.stats.distance) + dist_offset, color = 'red')

	if show_S:
		for tr in st_synS_select:
			dist_offset = tr.stats.distance
			ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)#	These lines used to cause a crash in Spyder
			if norm_each:
				plt.plot(ttt, (tr.data - np.median(tr.data))*plot_fac /(tr.data.max()
					- tr.data.min()) + dist_offset, color = 'green')
			else:
				plt.plot(ttt, (tr.data * plot_fac * tr.stats.distance) + dist_offset, color = 'green')

	if show_H_simp:
		for tr in st_synH_simp_select:
			dist_offset = tr.stats.distance
			ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)#	These lines used to cause a crash in Spyder
			if norm_each:
				plt.plot(ttt, (tr.data - np.median(tr.data))*plot_fac /(tr.data.max()
					- tr.data.min()) + dist_offset, color = 'orange')
			else:
				plt.plot(ttt, (tr.data * plot_fac * tr.stats.distance) + dist_offset, color = 'orange')

	if show_S_simp:
		for tr in st_synS_simp_select:
			dist_offset = tr.stats.distance
			ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)#	These lines used to cause a crash in Spyder
			if norm_each:
				plt.plot(ttt, (tr.data - np.median(tr.data))*plot_fac /(tr.data.max()
					- tr.data.min()) + dist_offset, color = 'purple')
			else:
				plt.plot(ttt, (tr.data * plot_fac * tr.stats.distance) + dist_offset, color = 'purple')

	plt.xlabel('Time (s)')
	plt.ylabel('Epicentral distance from event (km)')
#	plt.title(fname1[8:18] + ' vs ' + fname2[8:18])
	plt.title(date_label + ' ' + compo + ' ' + motion + ' component: data vs synthetics')
	plt.show()

	#  Save processed files
#	fname1 = 'Pro_Files/HD' + date_label1 + 'sel.mseed'
#	fname2 = 'Pro_Files/HD' + date_label2 + 'sel.mseed'
#	st1good.write(fname1,format = 'MSEED')
#	st2good.write(fname2,format = 'MSEED')

	elapsed_time_wc = time.time() - start_time_wc
	print('This job took ' + str(elapsed_time_wc) + ' seconds')
	os.system('say "Done"')