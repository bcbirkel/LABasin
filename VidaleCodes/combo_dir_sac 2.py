#!/usr/bin/env python
# This programs reads all SAC files in a directory, then writes out an mseed file.
# John Vidale 9/2020

def combo_dir_sac(dir_in, labl):
    import os
    from obspy import read, Trace, Stream
    import glob

    os.chdir(dir_in)
    file_list = glob.glob(labl + '*')
    print('Number of files ' + str(len(file_list)) + '; first -  ' + file_list[0])
    st = Stream()
    for sgrams in file_list:
        new_one = read(sgrams)
        st += new_one
        # print('Channel is ' + str(new_one[0].stats.network))
    cnt = len(st)
    print('Number of output traces ' + str(cnt))
    if cnt > 0:
        st.write(labl + '.mseed', format='MSEED')
