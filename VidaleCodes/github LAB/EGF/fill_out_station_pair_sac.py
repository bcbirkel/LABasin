#!/usr/bin/env python
# This programs deals with a all components of all stations.
# John Vidale 2/2019

import os
from obspy import Stream, Trace
from obspy import read

os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/vidale/Documents/GitHub/LAB/EGF')
#from run_1_all_to_obspy    import run1convert
os.chdir('/Users/vidale/Documents/PyCode/LAB/Spare/SAC')

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
			if sta == sta1:  # no auto-correlation
				x = 1
				print('Sta and Sta1 are same ' + sta + '  ' + sta1)
			else:
				# form the strings with the file names
				fname1 = comp + '/' + sta  + '_' + sta1 + '.sac'
				fname2 = comp + '/' + sta1 + '_' + sta  + '.sac'
				exists1 = os.path.isfile(fname1)
				exists2 = os.path.isfile(fname2)

				st = Stream()
				if exists1 and exists2:
					print('Both files are present, ' + fname1 + ' and ' + fname2)
				elif exists1:
					st = read(fname1)
					st[0].data = st[0].data[::-1] # flip time series
					st.write(fname2,format = 'sac')
	#				print('File ' + fname1)
				elif exists2:
					st = read(fname2)
					st[0].data = st[0].data[::-1] # flip time series
					st.write(fname1,format = 'sac')
	#				st.data = st.data.reverse()
	#				print('File ' + fname2)
				else:
					print('Neither file is present, ' + fname1 + ' nor ' + fname2)
