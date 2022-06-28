#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 12:08:45 2022

@author: bcbirkel
"""

import pickle
import numpy as np
import struct
import time


mat = ('planedisplacements.0')

with open(mat, mode='rb+') as f:
    fileContent = f.read()
    # Invoke the above lazy file read data function.
    # header = f.read(60)
    tshead = struct.unpack("ddddddddddiidi", fileContent[:100])
    print("Header: " + str(tshead))
    
    Aux = tshead[0:7]
    LengthX = tshead[8]
    LengthY = tshead[9]   
    PointsDip = tshead[10]  
    PointsStrike = tshead[11]
    DeltaT = tshead[12] 
    Steps = tshead[13] 
    TimeSteps = round(200/(Steps*DeltaT))
#%%
# define variables
# fsize = nx*ny*3*nt # number of points
# spapts = nx*ny*3
# fpts = nx*ny # number of grid points

t0 = time.time()
ts_array = np.zeros((PointsStrike,TimeSteps), dtype=float)
count = 0

def read_file_in_chunks(file_object, chunk_size=int(PointsStrike)*3):
    while True:
        data = []
        # Just read chunk_size size data.
        f = file_object.read(chunk_size*8)
        
        if not f:
            # Break the loop.
            break
        
        data = struct.unpack("d"*chunk_size, f[0:chunk_size*8])

        yield data
        
# Open the big data file.
with open(mat, mode='rb+') as f:
    fileContent = f.read()
    # Invoke the above lazy file read data function.
    # header = f.read(60)
    tshead = struct.unpack("ddddddddddiidi", fileContent[:100])
    print("Header: " + str(tshead))
    for tslice in read_file_in_chunks(f):
        # Process the piece of data
        ts_array[:, count] = tslice
        count += 1
        if count % 1000 == 0:
            t1 = time.time()
            print("Elapsed time on data read: " + str(t1-t0) + " seconds. " + str(count) + " of 8000 timeslices done.")

# %%
text = np.fromfile(mat)

for line in mat:
        line = line.strip()
        # Skip comments
        if line.startswith("#") or line.startswith("%"):
            pieces = line.split()[1:]
            # Write header
            if len(pieces) >= 10:
                print("# her header: # %s %s %s %s\n" %
                                 (pieces[0], pieces[1], pieces[2], pieces[3]))

# Aux=  fread(Plane1, 8,'double')
# BoxCorners=reshape(Aux,4,2)
# LengthY=  fread(Plane1, 1,'double')
# LengthX=  fread(Plane1, 1,'double')
# PointsDip = fread(Plane1, 1,'int')
# PointsStrike    = fread(Plane1, 1,'int')
# DeltaT  = fread(Plane1, 1,'double')
# Steps      = fread(Plane1, 1,'int')
# TimeSteps = round(200/(Steps*DeltaT)+1)
# ncount=0;
# rcount=0;

# Aux = text[0:7]
# BoxCorners=...
# [[0, 475000], [0, 760000], [0, 475000], [760000, 0]]
# LengthY=  760000
# LengthX=  475000
# PointsDip = 1
# PointsStrike    = 132
# DeltaT  =  0.0500
# Steps      = 5
# TimeSteps = round(200/(Steps*DeltaT)+1)
# ncount=0;
# rcount=0;

# SynS1 = np.zeros(TimeSteps,PointsStrike)

# for n in range(TimeSteps):
   
#     for r in range(PointsStrike):
        
#         SynS1[n][r] = fread(Plane1, 1,'double');
#         SynS1[n][r] = fread(Plane1, 1,'double');
#         SynS1[n][r] = fread(Plane1, 1,'double');
#         rcount = rcount + 1;
#         riter = r;
#     end
#     ncount = ncount + 1;
#     niter = n;
# end
# fclose(Plane1);