#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 09:51:39 2020

@author: bcbirkel
"""

#when called, should have sta_code_js, sta_seq_js, sta_code_rt = get_stations()

def get_stations():
    
    import pandas as pd
    import os
    
    #Stewart syn - pull stations into list, mapped with station code and 
    print('check')
    os.chdir('/Users/bcbirkel/Documents/Research/LABasin/stewart_synthetics')
    
    csv_file = 'stationtable11252019ngaw2iriscsmed_forquery120819.csv'
    
    df = pd.read_csv(csv_file)
    #sta_id = pd.Series(df.station_id_no_nga, dtype=float)
    sta_code_js = df.station_id_no_nga
    sta_code_js = pd.Series.tolist(sta_code_js)
    sta_seq_no = df.station_sequence_no
    sta_seq_no = pd.Series.tolist(sta_seq_no)
    
    print(type(sta_code_js))
    #print(sta_id)
    
    #Taborda syn - pull stations
    
    os.chdir('/Users/bcbirkel/Documents/Research/LABasin/PyCode LAB/ricardo')
    
    csv_file_2 = 'ricardo_full-station-list.csv'
    
    df = pd.read_csv(csv_file_2)
    #sta_id = pd.Series(df.station_id_no_nga, dtype=float)
    sta_code_rt = df.Code
    sta_code_rt = pd.Series.tolist(sta_code_rt)
    
    update_sta_code_js = []
    update_sta_seq_no_js = []
    update_sta_code_rt = []
    
    #for sta in len(sta_code_rt):
    #    sta_code_rt = sta_code_rt.append(int(sta))
    
    for sta_ii in range(len(sta_code_rt)):
        for sta_jj in range(len(sta_code_js)):
            if sta_code_rt[sta_ii] == sta_code_js[sta_jj]:
                #print('check')
                update_sta_code_js.append(sta_code_js[sta_jj])
                update_sta_seq_no_js.append(sta_seq_no[sta_jj])
                update_sta_code_rt.append(sta_code_rt[sta_ii])
    print(update_sta_code_rt)
    
    #sc = open('sta_crossover.txt', 'w')
    #sc.write(update_sta_code_rt)
    #sc.close()
    
    os.chdir('/Users/bcbirkel/Documents/Research/LABasin/stewart_synthetics')
    return update_sta_code_js, update_sta_seq_no_js, update_sta_code_rt
                