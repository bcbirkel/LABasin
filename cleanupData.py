#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 13:47:23 2022

@author: bcbirkel
"""

import shutil
from os import listdir

for ii in range(4):
 
    event_no = ii
    events = ['lahabra_2014', 'chatsworth_2007', 'chinohills_2008', 'inglewood_2009', 'beverlyhills_2001']
    event = events[event_no]
    dataPath = './CompiledEvents/' + event + '/IRIS/'
    outPath = './CompiledEvents/' + event + '/allData/'
    
    if event_no == 0:
        event_title = 'lahabra'
        event_lat = 33.9325
        event_lon = -117.9158
    elif event_no == 1:
        event_title = 'beverlyhills'
        event_lat = 34.0541
        event_lon = -118.3929
    elif event_no == 2:
        event_title = 'chinohills'
        event_lat = 33.9465
        event_lon = -117.7667
    elif event_no == 3:
        event_title = 'inglewood'
        event_lat = 33.9377
        event_lon = -118.3357
    elif event_no == 4:
        event_title = 'chatsworth'
        event_lat = 34.2983
        event_lon = -118.6255
        
    files = listdir(dataPath)
    keep = []
    
    bands = ['B','L','V','H']
    # instrs = ['H','L']
    instrs =['N']
    orients = ['N','E','Z']
    
    for f in files:
        for band in bands:
            for instr in instrs:
                for orient in orients:
                    code = '.' + band + instr + orient
                    if code in f:
                        keep.append(f)


    for file in keep:
        shutil.copyfile(dataPath + file, outPath + file)