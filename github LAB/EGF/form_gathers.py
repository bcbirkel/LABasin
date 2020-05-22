#!/usr/bin/env python
# input is set of hinet traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# This programs deals with a single event.
# John Vidale 2/2019

import os
from obspy import Stream
from obspy import read

sum = 1

os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/vidale/Documents/GitHub/LAB/EGF')
#from run_1_all_to_obspy    import run1convert
os.chdir('/Users/vidale/Documents/PyCode/LAB/Spare/SAC')

sta_test = ['BHP','BRE']
sta_all  = ['BHP','BRE','CAC','CBC','DEC','DJJ','DLA','FMP','FUL','GR2','GVR',
			'HLL','KIK','LAF','LCG','LGB','LLS','LTP','MIKB','MWC','OGC','OLI',
			'PASC','PDR','PSR','RHC2','RIO','RPV','RUS','SAN','SMF2','SRN',
			'STS','WLT','WNS','WTT2','USC']
comp     = ['ZZ','ZR','ZT','RZ','RR','RT','TZ','TR','TT']

print(str(len(sta_all)) + ' stations in sta_all ' + str(len(comp)) + ' components ')

#%% Choose appropriate station list
sta_list = sta_all  # short list for testing
#comp = ['ZZ']
#sta_list = sta_all

#%% read in all pairs with sta1 for comp1 component
for comp_one in comp:
		for sta in sta_list:
			st_all = Stream()
			counter = 0
			for sta1 in sta_list:
				if sta != sta1:  # no auto-correlation
					# form the strings with the file names
					if sum == 1:
						fname = comp_one + '_sum/' + sta  + '_' + sta1 + '_sum.sac'
					else:
						fname = comp_one + '/' + sta  + '_' + sta1 + '.sac'
					exists = os.path.isfile(fname)
					if exists:
						st = Stream()
						st = read(fname)
						counter += 1
	#				print('Station name ' + str(sta_all[0].stats.station))
						st_all.append(st[0])
			print('Station ' + sta + ' comp ' + comp_one + ', wrote out ' + str(counter) + '  ' + str(len(st_all)) + ' traces')
		#	print('First Station ' + str(sta_all[0].stats.station))
			if sum == 1:
				outname = '/Users/vidale/Documents/PyCode/LAB/Mseed/' + comp_one + '/' + sta + '_' + comp_one + '_sum.mseed'
			else:
				outname = '/Users/vidale/Documents/PyCode/LAB/Mseed/' + comp_one + '/' + sta + '_' + comp_one + '.mseed'
			st_all.write(outname,format = 'MSEED')