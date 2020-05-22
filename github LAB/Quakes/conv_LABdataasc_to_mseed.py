
#!/usr/bin/env python
# This programs reads in ASCII then writes out an mseed file.  Units are cm.
# John Vidale 6/2019

def con_dataasc_mseed(event_no, delta, motion, taper, taper_fac):

	import numpy as np
	from obspy import UTCDateTime, read, Trace, Stream
	import os
	import matplotlib.pyplot as plt
	import datetime
	import time
	os.chdir('/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/processed-data/' + event_no)

	# find event details for origin time, only use in writing output filenames
	ev_file = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/event_list.txt'
	file_ev = open(ev_file, 'r')
	for line in file_ev:           # pull numbers off the rest of the lines
		split_line = line.split()
		event = split_line[0]
		if event == event_no:
			t1           = UTCDateTime(split_line[5])
			date_label1  = split_line[5][0:10]
			year1        = split_line[5][0:4]

	print(event_no, str(t1) + ' ' + date_label1 + ' ' + year1)

	st_ev_ve = Stream()
	st_ev_vn = Stream()
	st_ev_vz = Stream()
	cnt_no  = 0
	cnt_yes = 0

	sta_file = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/ricardo_stations.txt'
	file_sta = open(sta_file, 'r')
	line = file_sta.readline()      # read first line to skip header information
	for line in file_sta:           # pull numbers off the rest of the lines
		split_line = line.split()
		sta_num   = split_line[0]
		netw      = split_line[2]
		staname      = split_line[3]
		latitude  = split_line[4]
		longitude = split_line[5]

		good_file = 'bad'
		check_Z = "no"  # will Z need timing correction?
		if netw == 'CE' or netw == 'CI' or netw == 'NP' or netw == 'WR' or netw == 'ZY':
			fname = netw + '.' + staname + '.HN.her'  # prefer HN channels
			exists = os.path.isfile(fname)
			if exists:
				good_file = fname
				chan = 'HN'
				check_Z = 'HN'
			else:
				fname = netw + '.' + staname + '.BH.her'  # BH are next best
				exists = os.path.isfile(fname)
				if exists:
					good_file = fname
					chan = 'BH'
					check_Z = 'BH'
				else:
					fname = netw + '.' + staname + '.V1.her'  # then V1
					exists = os.path.isfile(fname)
					if exists:
						good_file = fname
						chan = 'V1'
					else:
						fname = netw + '.' + staname + '.V2.her'  # then V2
						exists = os.path.isfile(fname)
						if exists:
							good_file = fname
							chan = 'V2'
						else:
							fname = netw + '.' + staname + '.HNN.00.ascii'  # finally these
							exists = os.path.isfile(fname)
							if exists:
								good_file = fname
								chan = 'HN'
							else:
								fname = netw + '.' + staname + '.HNN.01.ascii'  # finally these
								exists = os.path.isfile(fname)
								if exists:
									good_file = fname
									chan = 'HN'
								else:
									fname = netw + '.' + staname + '.HNN.02.ascii'  # finally these
									exists = os.path.isfile(fname)
									if exists:
										good_file = fname
										chan = 'HN'
									else:
										fname = netw + '.' + staname + '.HNN.05.ascii'  # finally these
										exists = os.path.isfile(fname)
										if exists:
											good_file = fname
											chan = 'HN'
										else:
											fname = netw + '.' + staname + '.HNN.06.ascii'  # finally these
											exists = os.path.isfile(fname)
											if exists:
												good_file = fname
												chan = 'HN'
											else:
												fname = netw + '.' + staname + '.HNN.10.ascii'  # finally these
												exists = os.path.isfile(fname)
												if exists:
													good_file = fname
													chan = 'HN'
												else:
													fname = netw + '.' + staname + '.HNN.30.ascii'  # finally these
													exists = os.path.isfile(fname)
													if exists:
														good_file = fname
														chan = 'HN'
		if good_file == 'bad':
			cnt_no = cnt_no + 1
		else:
			cnt_yes = cnt_yes + 1

		if good_file != 'bad':  # if a trace exists for this station, get it

			if(check_Z != 'no'):  # if Z timing is likely in error, wrongly assigned the E start time
				fnameZ = 'ascii/' + event_no + '.' + netw + '.' + staname + '.' + check_Z + 'Z' + '.ascii'
				fnameE = 'ascii/' + event_no + '.' + netw + '.' + staname + '.' + check_Z + 'E' + '.ascii'
#				existsZ = os.path.isfile(fnameZ)
#				existsN = os.path.isfile(fnameN)

				fileZ = open(fnameZ, 'r')
				fileE = open(fnameE, 'r')
				lineZ = fileZ.readline()  # first line has start time of trace
				lineE = fileE.readline()  # first line has start time of trace
				split_lineZ = lineZ.split()
				split_lineE = lineE.split()
				tZ = UTCDateTime(split_lineZ[4])
				tE = UTCDateTime(split_lineE[4])
				Z_shift = tZ - tE  # shift to correct Z timing error
				print('file ' + fnameZ + ' ' + str(tZ) + '  ' + str(tE) + ' has shift ' + str(Z_shift))

			file = open(good_file, 'r')

			line = file.readline()  # first line has start time of trace
			split_line = line.split()

			if chan == 'V1' or chan == 'V2':  # Hellish date format in V1
				t1 = datetime.datetime.strptime(split_line[4],'%m/%d/%y,%H:%M:%S.%f')
#				print('Station: ' + staname + ' ' + split_line[4] +  ' datetime: ' + str(t1))
			else:
				t1           = UTCDateTime(split_line[4])
#				print('Station: ' + staname + ' ' + split_line[4] +  ' datetime: ' + str(t1))

			line = file.readline()      # skip header, second line
			timeline = []
			datave = []
			datavn = []
			datavz = []
			for line in file:           # pull velocity numbers off the rest of the lines
				split_line = line.split()
				timeline.append(float(split_line[0]))
				if motion == 'disp':  # get displacement values
					datave.append(float(split_line[1]))
					datavn.append(float(split_line[2]))
					datavz.append(float(split_line[3]))
				elif motion == 'vel':  # get velocity values
					datave.append(float(split_line[4]))
					datavn.append(float(split_line[5]))
					datavz.append(float(split_line[6]))
				elif motion == 'acc':  # get acceleration values
					datave.append(float(split_line[7]))
					datavn.append(float(split_line[8]))
					datavz.append(float(split_line[9]))
				else:
					print('error: asking for unknown component of motion')

			# convert lists to arrays
			ve_array = np.array(datave)
			vn_array = np.array(datavn)
			vz_array = np.array(datavz)

			# calculate conversion to 10 sps, adjacent samples have round-off in Ricardos files
			dt = (timeline[10] - timeline[0])/10.
			dec_fac = delta/dt
#			print(netw + staname + 'dec_factor is ' + str(dec_fac))

			# Fill header attributes
			stats_e = {'network': netw, 'station': staname, 'location': '',
			         'channel': chan, 'npts': len(ve_array), 'delta': dt,
			         'mseed': {'dataquality': 'D'}}
			stats_n = {'network': netw, 'station': staname, 'location': '',
			         'channel': chan, 'npts': len(vn_array), 'delta': dt,
			         'mseed': {'dataquality': 'D'}}
			stats_z = {'network': netw, 'station': staname, 'location': '',
			         'channel': chan, 'npts': len(vz_array), 'delta': dt,
			         'mseed': {'dataquality': 'D'}}

			# set current time
			stats_e['starttime'] = t1
			stats_n['starttime'] = t1
			if(check_Z == 'no'):
				stats_z['starttime'] = t1
			else:
				stats_z['starttime'] = t1 + Z_shift  # correct Ricardo's error

			st_ve = Stream([Trace(data=ve_array, header=stats_e)])
			st_vn = Stream([Trace(data=vn_array, header=stats_n)])
			st_vz = Stream([Trace(data=vz_array, header=stats_z)])

			st_ve.detrend()
			st_vn.detrend()
			st_vz.detrend()
			if taper == 1:
				st_ve.taper(taper_fac)
				st_vn.taper(taper_fac)
				st_vz.taper(taper_fac)

			int_dec_fac = round(dec_fac) # decimate to 10 sps
			if abs((int_dec_fac/dec_fac) - 1) > 0.001:
				print(netw + ' ' + staname + ': check decimation factor, computed is ' + str(dec_fac) + ', but rounded is ' + str(int_dec_fac))
			if int_dec_fac == 20: # decimate only deals with factors 16 or less
				st_ve.decimate(5)
				st_vn.decimate(5)
				st_vz.decimate(5)
				st_ve.decimate(4)
				st_vn.decimate(4)
				st_vz.decimate(4)
			elif int_dec_fac < 16:
				st_ve.decimate(int_dec_fac)
				st_vn.decimate(int_dec_fac)
				st_vz.decimate(int_dec_fac)
			else:
				print('re-code - decimation factor, computed is ' + str(int_dec_fac) + ', bigger than obspy limit of 16')

#			st_ve.taper(0.3)
#			st_vn.taper(0.3)
#			st_vz.taper(0.3)
			st_ev_ve.append(st_ve[0])
			st_ev_vn.append(st_vn[0])
			st_ev_vz.append(st_vz[0])

	#st_z.append(st_e[0])
	if motion == 'disp':
		pathn = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_data/disp/'
	elif motion == 'vel':
		pathn = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_data/vel/'
	elif motion == 'acc':
		pathn = '/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo/Mseed_data/acc/'
	else:
		print('error: trying to write unknown component of motion')
	st_ev_ve.write(pathn + 've' + '_' + event_no + '.mseed', format='MSEED')
	st_ev_vn.write(pathn + 'vn' + '_' + event_no + '.mseed', format='MSEED')
	st_ev_vz.write(pathn + 'vz' + '_' + event_no + '.mseed', format='MSEED')
	#st_ev_vz.plot()

# debugging lines
#	for tr in st_ev_ve:
#		if tr.stats.station == '24861' or tr.stats.station == 'VCS':
#			print('Station: ' + tr.stats.station + ' datetime: ' + str(tr.stats.starttime))

	print('stations: ' + str(cnt_yes) + ' present, ' + str(cnt_no) + ' absent, out of 808 attempts')
	os.system('say "Done"')

	## Show that it worked, convert NumPy character array back to string
	#st1 = read("weather.mseed")
	#print(st1[0].data.tostring())