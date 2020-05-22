#!/usr/bin/env python
# this program translate Ricardo's synthetic and data files to mseed
# John Vidale 6/2019

import os

#os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/bcbirkel/Documents/Research/LABasin/github LAB/Quakes')

from conv_LABsynasc_to_mseed     import con_synasc_mseed
from conv_LABdataasc_to_mseed    import con_dataasc_mseed

#%% define common parameters
event_no = ['9064568', '9093975', '9096972', '9140050', '9644101', '9703873', '9716853','9753489',
	'9818433', '10216101', '10275733', '10370141', '10399889', '10403777', '10410337', '10530013',
	'10541957', '10972299', '13692644', '14000376', '14116972', '14155260', '14239184', '14312160',
	'14383980', '14494128', '14601172', '14738436', '15237281', '15481673']
vmodel = ['cvmhn', 'cvmhy', 'cvms400-100', 'cvms426-223']
# 0 - Harvard, no GTL
# 1 - Harvard, with GTL
# 2 - SCEC, simpler
# 3 - SCEC, full model

delta = 0.1
taper = 0
taper_fac = 0.05
motion = 'acc'

#%% Run over all synthetic event files and both models
#for evno in event_no:
#	for vmod in vmodel:
#		con_synasc_mseed(event_no = evno, vmodel = vmod, delta = delta)

# individual test runs for debugging
con_synasc_mseed(event_no = event_no[27], vmodel = vmodel[1], delta = delta, motion = 'acc')
con_synasc_mseed(event_no = event_no[27], vmodel = vmodel[2], delta = delta, motion = 'acc')
con_synasc_mseed(event_no = event_no[27], vmodel = vmodel[0], delta = delta, motion = 'acc')

# Runs over all processed data files
#for evno in event_no:
#	con_dataasc_mseed(event_no = evno, delta = delta, motion = motion, taper = taper, taper_fac = taper_fac)

# individual test runs for debugging
#con_dataasc_mseed(event_no = event_no[27], delta = delta, motion = motion, taper = taper, taper_fac = taper_fac)

#for evno in event_no:
#	con_synasc_mseed(event_no = evno, vmodel = vmodel[0], delta = delta)
#for evno in event_no:
#	con_synasc_mseed(event_no = evno, vmodel = vmodel[1], delta = delta)
#for evno in event_no:
#	con_synasc_mseed(event_no = evno, vmodel = vmodel[2], delta = delta)
#for evno in event_no:
#	con_synasc_mseed(event_no = evno, vmodel = vmodel[3], delta = delta)

# returns to a convenient directory when done
os.chdir('/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo')