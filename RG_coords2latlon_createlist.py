#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 15:57:47 2022

@author: bcbirkel
"""
import pickle

latlonmod = []
xymod = []
model_trans = []

# pull in latlon to xy coordinate transformation info
modelcoord = './model_coords_sc01-h0.100'
coord = open(modelcoord, 'r')
pairs = coord.readlines()
for line in pairs:
    split_line = line.split()
    lon = float(split_line[0])
    lat = float(split_line[1])
    xmod = int(split_line[2])
    ymod = int(split_line[3])
    latlonmod.append([lat,lon])
    xymod.append([xmod,ymod])
    model_trans.append([xmod,ymod,lat,lon])

print("coordinates stripped") #check

# %% ##### SORT LAT,LON TO FIT ORDERING OF DATA FROM SIMULATION #####

xymod_sorted = []

# skips to every 10th point in coord transform file to match sim grid
model_trans_skip = model_trans[::10]

#pulls out only useful points from coord transform file
for pt in model_trans_skip:
    if pt[0]%10 == 0 and pt[1]%10 == 0:
        xymod_sorted.append(pt)

print("lat lon indices sorted")

with open("RG_latlon", "wb") as f:   #Pickling
     pickle.dump(xymod_sorted, f)
     

# CHECK

with open("RG_latlon", "rb") as f:   # Unpickling
    test = pickle.load(f)