#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 17:13:10 2019

@author: vidale
"""
def re_get_Z(select_event=0):
	from obspy.clients.fdsn import Client
	import matplotlib.pyplot as plt
	from obspy import read_events, read
	from obspy import read_inventory
	from obspy import Stream, Trace
	from obspy import UTCDateTime
	import os
	client = Client('SCEDC')

	os.environ['PATH'] += os.pathsep + '/usr/local/bin'
	os.chdir('/Users/vidale/Documents/PyCode/LAB/Spare')

#	select_event = 15
	chan_type = 'HNZ' # e.g., BHZ
	network_sel = 'CE'  # CE has four traces, but won't deconvolve
	#network_sel = 'CI,NP'
	min_lat = 33.75
	max_lat = 34.2
	min_lon = -118.5
	max_lon = -117.75
	start_buff = 0
	end_buff = 50
	st = Stream()

	if select_event > 15:
		fname_inv = 'LAB.QUAKEML2'
		LAB = read_events(fname_inv, format='QUAKEML')
		if select_event == 16:
			t = LAB[1].origins[0].time
		elif select_event == 17:
			t = LAB[5].origins[0].time
		elif select_event == 18:
			t = LAB[14].origins[0].time
		elif select_event == 19:
			t = LAB[13].origins[0].time
	else:
		fname_inv = 'LAB.QUAKEML'
		LAB = read_events(fname_inv, format='QUAKEML')
		t = LAB[select_event].origins[0].time

		#print('event:',LAB)
		#plt.style.use('ggplot')
		#plt.rcParams['figure.figsize'] = 12, 8
		#LAB.plot(projection = 'local', resolution = 'h')

		#%% Make inventory of all stations in box recording this channel

	print(str(t))
	s_t = t - start_buff
	e_t = t + end_buff

	inventory = client.get_stations(starttime = s_t, endtime = e_t,
						channel=chan_type, level='station', network=network_sel,
						minlatitude  = min_lat, maxlatitude  = max_lat,
						minlongitude = min_lon, maxlongitude = max_lon)

#	inventory = client.get_stations(starttime = s_t, endtime = e_t,
#						channel=chan_type, level='response', network=network_sel,
#						minlatitude  = min_lat, maxlatitude  = max_lat,
#						minlongitude = min_lon, maxlongitude = max_lon)

	#print(inventory)
	print('inventory has ' + str(len(inventory)) + ' networks recording data')
#	print(inventory)

	#for network in inventory:
	#	sta_cnt = 0
	#	for station in network:
	#		sta_cnt += sta_cnt
	#	print('Network ' + str(network) + ' has ' + str(sta_cnt) + ' stations to try')
	#inventory.plot(projection = 'local', resolution = 'h')  # not working

	#%% Check inventory of stations for traces at time of event
	cnt_try = 0
	cnt_got = 0
	for network in inventory:
		for station in network:
			print('trying ' + network.code + '  ' + station.code)
			if cnt_try % 20 == 0:
				print('Try ' + str(cnt_try) + ' got ' + str(cnt_got) + ' sgrams ' + str(len(st)))
			cnt_try += +1
			try:
				st += client.get_waveforms(network.code, station.code, location='*',channel=chan_type, starttime=s_t, endtime = e_t, attach_response=False)
				cnt_got += 1
			except:
				pass
'''

	print(str(cnt_try) + ' stations examined, ' + str(cnt_got) + ' have data, ' + str(len(st)) + ' traces extracted  ')
	fname = 're_event' + str(select_event) + '/event' + str(select_event) + 'Z_all.mseed'
	st.write(fname,format = 'MSEED')
	print('1  Wrote ' + str(len(st)) + ' stations')

	for tr in st:
		print('Station ' + tr.stats.station + ' channel ' + tr.stats.channel)

	tr = Trace()
	hnz_chosen = 0
	ehz_chosen = 0
	hhz_chosen = 0
	hlz_chosen = 0
	bhz_chosen = 0
	for tr in st:
		if tr.stats.channel == 'HNZ':
			hnz_chosen += 1
		if tr.stats.channel == 'EHZ':
			ehz_chosen += 1
		if tr.stats.channel == 'HHZ':
			hhz_chosen += 1
		if tr.stats.channel == 'HLZ':
			hlz_chosen += 1
		if tr.stats.channel == 'BHZ':
			bhz_chosen += 1

	print('Total channels ' + str(len(st)) + ' - HNZ, EHZ, HHZ, HLZ, BHZ have '
		   + str(hnz_chosen) + ' ' + str(ehz_chosen) + ' '
		   + str(hhz_chosen) + ' ' + str(hlz_chosen) + ' ' + str(bhz_chosen))

	for tr in st:
		if (tr.stats.network != 'CE') and (tr.stats.station != 'BVH') and (tr.stats.station != 'LAX'):
#		if tr.stats.network != 'CE' and not (tr.stats.station == 'BVH' and tr.stats.channel == 'EHZ'):
			print('Station ' + tr.stats.station + ' channel ' + tr.stats.channel)
			tr.remove_response(water_level=40, inventory=inventory, output='ACC')

	fname = 're_event' + str(select_event) + '/event' + str(select_event) + 'Z_decon.mseed'
	st.write(fname,format = 'MSEED')
	print('2  Wrote ' + str(len(st)) + ' stations')

	tr2 = Trace()
	st_chosen = Stream()
	hhz_chosen = 0
	ehz_chosen = 0
	hnz_chosen = 0
	hlz_chosen = 0
	bhz_chosen = 0
	for tr in st:
		if tr.stats.channel == 'HNZ':
			st_chosen += tr
			hnz_chosen += 1
		elif tr.stats.channel == 'EHZ':  # write EHZ if present and BHZ is not present
			skip = 0
			for tr2 in st:
				if tr2.stats.channel == 'HNZ' and tr2.stats.station == tr.stats.station:
					skip = 1
			if skip == 0:
				st_chosen += tr
				ehz_chosen += 1
		elif tr.stats.channel == 'HHZ':
			skip = 0
			for tr2 in st:
				if tr2.stats.channel == 'HNZ' and tr2.stats.station == tr.stats.station:
					skip = 1
				if tr2.stats.channel == 'EHZ' and tr2.stats.station == tr.stats.station:
					skip = 1
			if skip == 0:
				st_chosen += tr
				hhz_chosen += 1
		elif tr.stats.channel == 'HLZ':
			skip = 0
			for tr2 in st:
				if tr2.stats.channel == 'HNZ' and tr2.stats.station == tr.stats.station:
					skip = 1
				if tr2.stats.channel == 'EHZ' and tr2.stats.station == tr.stats.station:
					skip = 1
				if tr2.stats.channel == 'HHZ' and tr2.stats.station == tr.stats.station:
					skip = 1
			if skip == 0:
				st_chosen += tr
				hlz_chosen += 1
		elif tr.stats.channel == 'BHZ':
			skip = 0
			for tr2 in st:
				if tr2.stats.channel == 'HNZ' and tr2.stats.station == tr.stats.station:
					skip = 1
				if tr2.stats.channel == 'EHZ' and tr2.stats.station == tr.stats.station:
					skip = 1
				if tr2.stats.channel == 'HHZ' and tr2.stats.station == tr.stats.station:
					skip = 1
				if tr2.stats.channel == 'HLZ' and tr2.stats.station == tr.stats.station:
					skip = 1
			if skip == 0:
				st_chosen += tr
				bhz_chosen += 1

	for tr in st_chosen:
		print(tr.stats.station + '  ' + tr.stats.channel)

	print('Chosen - HNZ, EHZ, HHZ, HLZ, BHZ have ' + str(hnz_chosen) + ' ' + str(ehz_chosen) +
		   ' ' + str(hhz_chosen) + ' ' + str(hlz_chosen) + ' ' + str(bhz_chosen))
	print(str(len(st_chosen)) + ' traces in dataset')
	fname = 're_event' + str(select_event) + '/event' + str(select_event) + 'Z_chosen.mseed'
	st_chosen.write(fname,format = 'MSEED')
	print('3  Wrote ' + str(len(st_chosen)) + ' stations')
'''