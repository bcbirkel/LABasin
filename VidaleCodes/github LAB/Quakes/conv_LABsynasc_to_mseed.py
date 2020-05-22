#!/usr/bin/env python
# This programs reads in ASCII then writes out an mseed file.
#  units are meters
# John Vidale 6/2019

def con_synasc_mseed(event_no, vmodel, delta, motion):

	import numpy as np
	from obspy import UTCDateTime, read, Trace, Stream
	import os
	import matplotlib.pyplot as plt

	os.chdir('/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/spare_files/synthetics/' + event_no + '/outputfiles/stations/'+ vmodel)
	delta = 0.1

	# find event details for origin time
	ev_file = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/event_list.txt'
	file_ev = open(ev_file, 'r')
	for line in file_ev:           # pull numbers off the rest of the lines
		split_line = line.split()
		event = split_line[0]
		if event == event_no:
			t1           = UTCDateTime(split_line[5])
			date_label1  = split_line[5][0:10]
			year1        = split_line[5][0:4]

	print(event_no, vmodel, str(t1) + ' ' + date_label1 + ' ' + year1)

	#%% Open station location file
	sta_file = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/ricardo_stations.txt'
	file_st = open(sta_file, 'r')
	line = file_st.readline()      # read first line to skip header information
	lines = file_st.readlines()
	print(str(len(lines)) + ' stations read from ' + sta_file)

	# Load station coords into arrays, many more stations than used
#	line_cnt = 0; # keep track of the line number in station file
#	st_line  = []
	st_num   = [] # keep track of the station index
	st_netw  = []
	st_name  = []
	st_lat   = []
	st_lon   = []
	for line in lines:
#		st_line.append(line_cnt)
#		line_cnt = line_cnt+1
		split_line = line.split()
		st_num.append(  split_line[0])
		st_netw.append( split_line[2])
		st_name.append( split_line[3])
		lat = float(split_line[4])
		st_lat.append(lat)
		lon = float(split_line[5])
		st_lon.append(lon)
	print('number of stations in list is ' + str(len(st_num)))

	#%%  Read in seismograms
	st_ev_ve = Stream()
	st_ev_vn = Stream()
	st_ev_vz = Stream()
	cnt_no  = 0
	cnt_yes = 0
	for cnt in range(len(st_num)):  # cnt is line number in file, st_num is string with file name number
		labelcnt = st_num[cnt]
		eqsta_file = 'station.' + labelcnt
		file_exists = os.path.isfile(eqsta_file)
		if file_exists:
			file = open(eqsta_file, 'r')
			cnt_yes = cnt_yes + 1
			line = file.readline()      # skip header line, which has no useful information
			datave = []
			datavn = []
			datavz = []
			for line in file:           # pull numbers off the rest of the lines
				split_line = line.split()
				if motion == 'disp':  # get displacement values
					datavn.append(float(split_line[1]))
					datave.append(float(split_line[2]))
					datavz.append(float(split_line[3]))
				elif motion == 'vel':  # get velocity values
					datavn.append(float(split_line[4]))
					datave.append(float(split_line[5]))
					datavz.append(float(split_line[6]))
				elif motion == 'acc':  # get acceleration values
					datavn.append(float(split_line[7]))
					datave.append(float(split_line[8]))
					datavz.append(float(split_line[9]))
				else:
					print('error: asking for unknown component of motion')
#			print('station ' + st_name[cnt] + ' ' + ' st_num ' + labelcnt + ' cnt ' + str(cnt))
			# convert list to array
			ve_array = np.array(datave)
			vn_array = np.array(datavn)
			vz_array = np.array(datavz)
			# Fill header attributes
			stats_e = {'network': st_netw[cnt], 'station': st_name[cnt], 'location': '',
			         'channel': 'HNE', 'npts': len(ve_array), 'delta': delta,
			         'mseed': {'dataquality': 'D'}}
			stats_n = {'network': st_netw[cnt], 'station': st_name[cnt], 'location': '',
			         'channel': 'HNN', 'npts': len(vn_array), 'delta': delta,
			         'mseed': {'dataquality': 'D'}}
			stats_z = {'network': st_netw[cnt], 'station': st_name[cnt], 'location': '',
			         'channel': 'HNZ', 'npts': len(vz_array), 'delta': delta,
			         'mseed': {'dataquality': 'D'}}

			# set current time
			stats_e['starttime'] = t1
			stats_n['starttime'] = t1
			stats_z['starttime'] = t1

			st_ve = Stream([Trace(data=ve_array, header=stats_e)])
			st_vn = Stream([Trace(data=vn_array, header=stats_n)])
			st_vz = Stream([Trace(data=vz_array, header=stats_z)])

			st_ev_ve.append(st_ve[0])
			st_ev_vn.append(st_vn[0])
			st_ev_vz.append(st_vz[0])
		else:
			cnt_no = cnt_no + 1

	if motion == 'disp':
		pathn = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/disp'
	elif motion == 'vel':
		pathn = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/vel'
	elif motion == 'acc':
		pathn = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_syn/acc'
	else:
		print('error: trying to write unknown component of motion')

	st_ev_ve.write(pathn + '/' + vmodel + '/ve' + '_' + event_no + '_' + vmodel + '.mseed', format='MSEED')
	st_ev_vn.write(pathn + '/' + vmodel + '/vn' + '_' + event_no + '_' + vmodel + '.mseed', format='MSEED')
	st_ev_vz.write(pathn + '/' + vmodel + '/vz' + '_' + event_no + '_' + vmodel + '.mseed', format='MSEED')

	#st_ev_vz.plot()

	os.system('say "Done"')