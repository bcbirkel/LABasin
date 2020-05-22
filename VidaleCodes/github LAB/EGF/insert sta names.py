#!/usr/bin/env python
# This programs deals with a all components of all stations.
# John Vidale 2/2019

import os
from obspy import Stream, Trace
from obspy import read

os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/vidale/Documents/GitHub/LAB/EGF')
#from run_1_all_to_obspy    import run1convert
os.chdir('/Users/vidale/Documents/PyCode/Spare/SAC')

sta_test = ['BHP','BRE']
sta_list  = ['BHP','BRE','CAC','CBC','DEC','DJJ','DLA','FMP','FUL','GR2','GVR',
			'HLL','KIK','LAF','LCG','LGB','LLS','LTP','MIKB','MWC','OGC','OLI',
			'PASC','PDR','PSR','RHC2','RIO','RPV','RUS','SAN','SMF2','SRN',
			'STS','WLT','WNS','WTT2','USC']
comp_list     = ['ZZ','ZR','ZT','RZ','RR','RT','TZ','TR','TT']

print(str(len(sta_list)) + ' stations in sta_list ' + str(len(comp_list)) + ' components ')

#%%  short list for testing
#sta_list = sta_test

#%% read in all pairs with sta1 for comp1 component
for comp in comp_list:
	for sta in sta_list:
		st_all = Stream()
		for sta1 in sta_list:
			if sta != sta1:  # no auto-correlation
				# form the strings with the file names
				fname1 = comp + '/' + sta  + '_' + sta1 + '.sac'
				fname2 = comp + '/' + sta1 + '_' + sta  + '.sac'
				fname1s = comp + '_new/' + sta  + '_' + sta1 + '.sac'
				fname2s = comp + '_new/' + sta1 + '_' + sta  + '.sac'
				exists1 = os.path.isfile(fname1)
				exists2 = os.path.isfile(fname2)

				if exists1 and exists2:
					st1 = Stream()
					st2 = Stream()
					st1 = read(fname1)
					st2 = read(fname2)
					st1[0].stats.station = sta1
					st2[0].stats.station = sta
#					print('Both files are present, ' + fname1 + ' and ' + fname2)
					st1.write(fname1s,format = 'sac')
					st2.write(fname2s,format = 'sac')
					fname3 = comp + '_sum/' + sta  + '_' + sta1 + '_sum.sac'
					fname4 = comp + '_sum/' + sta1 + '_' + sta  + '_sum.sac'
					for it in range(len(st1[0].data)):  # check points one at a time
						st1[0].data[it] = (st1[0].data[it] + st2[0].data[it])/2
						st2[0].data[it] =  st1[0].data[it]
					st1.write(fname3,format = 'sac')
					st2.write(fname4,format = 'sac')
				elif exists1 or exists2:
	#				print('File ' + fname1)
					print('In comp ' + comp + ', Only one file is present, ' + fname1 + ' or ' + fname2)
