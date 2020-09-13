#!/usr/bin/env python
# this program runs python code to combine sac files into a single mseed file
# John Vidale 9/2020
# edit Brianna Birkel 9/13/20

import os

from combo_dir_sac import combo_dir_sac

events = ['lahabra_2014', 'beverlyhills_2001', 'chatsworth_2007', 'chinohills_2008', 'inglewood_2009']
evids = ['15481673','9703873','14312160','14383980','10410337']

for ii in range(5):
    event = events[ii]

    dir_in = '/Users/bcbirkel/Documents/GitHub/LABasin/CompiledEvents/' + event + '/allData/'
    labl = 'Z'  # rotated data
    combo_dir_sac(dir_in = dir_in,labl = labl)
    labl = 'N'  # rotated data
    combo_dir_sac(dir_in = dir_in,labl = labl)
    labl = 'E'  # rotated data
    combo_dir_sac(dir_in = dir_in,labl = labl)
    labl = 'R'  # rotated data
    combo_dir_sac(dir_in = dir_in,labl = labl)
    labl = 'T'  # rotated data
    combo_dir_sac(dir_in = dir_in,labl = labl)
    
    dir_in = '/Users/bcbirkel/Documents/GitHub/LABasin/CompiledEvents/' + event + '/GravesSyn/CVM-S4/'
    labl = 'Z'  # rotated data
    combo_dir_sac(dir_in = dir_in,labl = labl)
    labl = 'N'  # rotated data
    combo_dir_sac(dir_in = dir_in,labl = labl)
    labl = 'E'  # rotated data
    combo_dir_sac(dir_in = dir_in,labl = labl)
    labl = 'R'  # rotated data
    combo_dir_sac(dir_in = dir_in,labl = labl)
    labl = 'T'  # rotated data
    combo_dir_sac(dir_in = dir_in,labl = labl)
    
    dir_in = '/Users/bcbirkel/Documents/GitHub/LABasin/CompiledEvents/' + event + '/GravesSyn/CVM-H/'
    labl = 'Z'  # rotated data
    combo_dir_sac(dir_in = dir_in,labl = labl)
    labl = 'N'  # rotated data
    combo_dir_sac(dir_in = dir_in,labl = labl)
    labl = 'E'  # rotated data
    combo_dir_sac(dir_in = dir_in,labl = labl)
    labl = 'R'  # rotated data
    combo_dir_sac(dir_in = dir_in,labl = labl)
    labl = 'T'  # rotated data
    combo_dir_sac(dir_in = dir_in,labl = labl)
