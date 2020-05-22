#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 17:13:10 2019

@author: vidale
"""
def get_E(select_event=0):
	from obspy.clients.fdsn import Client
	import matplotlib.pyplot as plt
	from obspy import read_events, read
	from obspy import read_inventory
	from obspy import Stream, Trace
	from obspy import UTCDateTime
	import os
	client = Client('SCEDC')

	os.environ['PATH'] += os.pathsep + '/usr/local/bin'
	os.chdir('/Users/vidale/Documents/PyCode/LAB/Mseed')

#	select_event = 12
	chan_type = 'EHE,HHE,HNE,HLE,BHE' # e.g., BHE
	network_sel = 'CI,CE,NP'  # CE has four traces, but won't deconvolve
	#network_sel = 'CI,NP'
	min_lat = 33.75
	max_lat = 34.2
	min_lon = -118.5
	max_lon = -117.75
	start_buff = 50
	end_buff = 300
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

		t = LAB[select_event].origins[0].time
	print(str(t))
	s_t = t - start_buff
	e_t = t + end_buff

	inventory = client.get_stations(starttime = s_t, endtime = e_t,
						channel=chan_type, level='response', network=network_sel,
						minlatitude  = min_lat, maxlatitude  = max_lat,
						minlongitude = min_lon, maxlongitude = max_lon)

	#print(inventory)
	print('inventory has ' + str(len(inventory)) + ' networks recording data')
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
			if cnt_try % 20 == 0:
				print('Try ' + str(cnt_try) + ' got ' + str(cnt_got) + ' sgrams ' + str(len(st)))
			cnt_try += +1
			try:
				st += client.get_waveforms(network.code, station.code, location='*',channel=chan_type, starttime=s_t, endtime = e_t, attach_response=True)
				cnt_got += 1
			except:
				pass

	print(str(cnt_try) + ' stations examined, ' + str(cnt_got) + ' have data, ' + str(len(st)) + ' traces extracted  ')
	fname = 'event' + str(select_event) + '/event' + str(select_event) + 'E_all.mseed'
	st.write(fname,format = 'MSEED')
	#st=read(fname)

	for tr in st:
		print('Station ' + tr.stats.station + ' channel ' + tr.stats.channel)

	tr = Trace()
	hne_chosen = 0
	ehe_chosen = 0
	hhe_chosen = 0
	hle_chosen = 0
	bhe_chosen = 0
	for tr in st:
		if tr.stats.channel == 'HNE':
			hne_chosen += 1
		if tr.stats.channel == 'EHE':
			ehe_chosen += 1
		if tr.stats.channel == 'HHE':
			hhe_chosen += 1
		if tr.stats.channel == 'HLE':
			hle_chosen += 1
		if tr.stats.channel == 'BHE':
			bhe_chosen += 1

	print('Total channels ' + str(len(st)) + ' - HNE, EHE, HHE, HLE, BHE have '
		   + str(hne_chosen) + ' ' + str(ehe_chosen) + ' '
		   + str(hhe_chosen) + ' ' + str(hle_chosen) + ' ' + str(bhe_chosen))

	for tr in st:
		if (tr.stats.network != 'CE') and (tr.stats.station != 'BVH') and (tr.stats.station != 'LAX'):
#		if tr.stats.network != 'CE':
			tr.remove_response(water_level=40, inventory=inventory, output='ACC')

	fname = 'event' + str(select_event) + '/event' + str(select_event) + 'E_decon.mseed'
	st.write(fname,format = 'MSEED')

	'''
	st=read('event14E_decon.mseed')
	'''

	tr2 = Trace()
	st_chosen = Stream()
	hhe_chosen = 0
	ehe_chosen = 0
	hne_chosen = 0
	hle_chosen = 0
	bhe_chosen = 0
	for tr in st:
		if tr.stats.channel == 'HNE':
			st_chosen += tr
			hne_chosen += 1
		elif tr.stats.channel == 'EHE':  # write EHE if present and BHE is not present
			skip = 0
			for tr2 in st:
				if tr2.stats.channel == 'HNE' and tr2.stats.station == tr.stats.station:
					skip = 1
			if skip == 0:
				st_chosen += tr
				ehe_chosen += 1
		elif tr.stats.channel == 'HHE':
			skip = 0
			for tr2 in st:
				if tr2.stats.channel == 'HNE' and tr2.stats.station == tr.stats.station:
					skip = 1
				if tr2.stats.channel == 'EHE' and tr2.stats.station == tr.stats.station:
					skip = 1
			if skip == 0:
				st_chosen += tr
				hhe_chosen += 1
		elif tr.stats.channel == 'HLE':
			skip = 0
			for tr2 in st:
				if tr2.stats.channel == 'HNE' and tr2.stats.station == tr.stats.station:
					skip = 1
				if tr2.stats.channel == 'EHE' and tr2.stats.station == tr.stats.station:
					skip = 1
				if tr2.stats.channel == 'HHE' and tr2.stats.station == tr.stats.station:
					skip = 1
			if skip == 0:
				st_chosen += tr
				hle_chosen += 1
		elif tr.stats.channel == 'BHE':
			skip = 0
			for tr2 in st:
				if tr2.stats.channel == 'HNE' and tr2.stats.station == tr.stats.station:
					skip = 1
				if tr2.stats.channel == 'EHE' and tr2.stats.station == tr.stats.station:
					skip = 1
				if tr2.stats.channel == 'HHE' and tr2.stats.station == tr.stats.station:
					skip = 1
				if tr2.stats.channel == 'HLE' and tr2.stats.station == tr.stats.station:
					skip = 1
			if skip == 0:
				st_chosen += tr
				bhe_chosen += 1

	for tr in st_chosen:
		print(tr.stats.station + '  ' + tr.stats.channel)

	print('Chosen - HNE, EHE, HHE, HLE, BHE have ' + str(hne_chosen) + ' ' + str(ehe_chosen) +
		   ' ' + str(hhe_chosen) + ' ' + str(hle_chosen) + ' ' + str(bhe_chosen))
	print(str(len(st_chosen)) + ' traces in dataset')
	fname = 'event' + str(select_event) + '/event' + str(select_event) + 'E_chosen.mseed'
	st_chosen.write(fname,format = 'MSEED')
