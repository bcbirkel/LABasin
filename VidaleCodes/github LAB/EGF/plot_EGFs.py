#!/usr/bin/env python
# input is set of hinet traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# This programs deals with a single event.
# John Vidale 2/2019

def plot_EGFs(comp='ZZ',sta='USC',freqmin=0.2,freqmax=1.0):
	import os
	from obspy import Stream
	from obspy import read
	import matplotlib.pyplot as plt
	from obspy.geodetics import gps2dist_azimuth
	import numpy as np

	os.environ['PATH'] += os.pathsep + '/usr/local/bin'
	#os.chdir('/Users/vidale/Documents/PyCode/LAB/Spare/SAC')

	plot_scale_fac = 5
	flip_time = 0
	max_plot_dist = 70
	max_plot_time = 200

	#freqmin = 0.5
	#freqmax = 2.0
	#comp = 'ZZ'
	#sta  = 'USC'
	summer  = 1
	basin_only = 1

	#basin_sta = ['USC', 'GR2', 'BHP', 'WTT2', 'LGB','LTP', 'STS', 'DLA', 'BRE', 'OGC', 'FUL', 'LLS', 'SAN']
	basin_sta = ['SMF2','PDR','LAF','USC','LCG','BHP','WTT','COO','WTT2','LGB','LTP','STS','DLA','BRE','OGC','FUL','LLS','SAN','LBW1','LBW2','5499','CRF']

	# Load station coords into arrays
	sta_file = '/Users/vidale/Documents/PyCode/LAB/sta_list.txt'
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
	if summer == 1:
		fname = '/Users/vidale/Documents/PyCode/LAB/Mseed/' + comp + '/' + sta + '_' + comp + '_sum.mseed'
	else:
		fname = '/Users/vidale/Documents/PyCode/LAB/Mseed/' + comp + '/' + sta + '_' + comp + '.mseed'
	st = read(fname)

	if sta in st_names:
		sta_index = st_names.index(sta)
		ref_lat = float(st_lats[sta_index])
		ref_lon = float(st_lons[sta_index])
	else:
		print(sta + ' is not in station list')

	# make figure
	if comp == 'ZZ':
		plt.close(7)
		fig_index = 7
	elif comp == 'TT':
		plt.close(8)
		fig_index = 8
	plt.figure(fig_index,figsize=(15,10))
	plt.xlim(-max_plot_time,max_plot_time)
	plt.ylim(0,max_plot_dist)
	for tr in st:
		if flip_time == 1:
			tr.data = tr.data[::-1] # flip time series
		tr.detrend()
		tr.taper(0.1)
		tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=False)

		if tr.stats.station in st_names:
			sta_index = st_names.index(tr.stats.station)
			stalat = float(st_lats[sta_index])
			stalon = float(st_lons[sta_index]) # look up lat & lon again to find distance
			distance = gps2dist_azimuth(stalat,stalon,ref_lat,ref_lon) # Get traveltimes again, hard to store
			print(str(distance[0]) + ' is distance')
		else:
			print(tr.stats.station + ' is not in station list')

		dist_offset = distance[0]/1000 # correct to km
		ttt = np.arange(len(tr.data)) * tr.stats.delta - 300

		if tr.stats.station in basin_sta:
			sta_index = basin_sta.index(tr.stats.station)  # just to test algorithm
	#		print(tr.stats.station + ' index is ' + str(sta_index))
			plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
				- tr.data.min()) + dist_offset, color = 'red')
			plt.text(s = tr.stats.station,x = max_plot_time*0.95,y = dist_offset +
				max_plot_dist*0.015, color = 'red')

		else:
			plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
				- tr.data.min()) + dist_offset, color = 'black')
			plt.text(s = tr.stats.station,x = max_plot_time*0.95,y = dist_offset +
				max_plot_dist*0.015, color = 'black')

	if comp == 'TT':
		plt.title('Love wave basin profile')
	else:
		plt.title('Rayleigh wave basin profile')
	plt.xlabel('Time (s)')
	plt.ylabel('Distance from ' + sta + ' (km)')
