#This program is designed to find all of the matches in two different fits
#flat files, finding the addresses where the objects in 1 match the objects in 2.
#rec_match and rec_match_srt only work well with hash arrays (they compare one array)
#so you will need to use mk_hsh_arr first to make PMF hash arrays for each fits
#flat file.

#Neither rec_match function should be used if you are just trying to find one or
#two matches. These are for large numbers of matches to be made.

#Code and ideas for rec_match lifted from https://www.followthesheep.com/?p=1366.

from astropy.io import fits
import numpy as np
import tmark
import progressBar as pb

#This is the very slow version. It's brute force. It can take a while. Not to be
#used with very large arrays of objects. Still needs PMF hash arrays.
def rec_match(rec1,rec2):
    tmark.tm('Starting Record Match')
    bool_rec1 = np.in1d(rec1,rec2)
    rec1_adr = np.arange(len(rec1))
    rec1_adr = rec1_adr[bool_rec1]

    tmark.tm('Finding Address Matches in rec2')
    rec2_adr = np.array([np.argwhere(rec2 == rec1[x]) for x in rec1_adr]).flatten()

    return rec1_adr,rec2_adr

#This works much faster than rec_match and will work with non-unique values.
#I'm not sure what it does with non-uniques though.
def rec_match_srt(rec1,rec2):
    rec1a = np.argsort(rec1)
    rec2a = np.argsort(rec2)
    tmark.tm('Starting Searchsorted')
    sort_left_rec1 = rec1[rec1a].searchsorted(rec2[rec2a],side='left')
    sort_right_rec1 = rec1[rec1a].searchsorted(rec2[rec2a],side='right')
    sort_left_rec2 = rec2[rec2a].searchsorted(rec1[rec1a],side='left')
    sort_right_rec2 = rec2[rec2a].searchsorted(rec1[rec1a],side='right')

    rec2_adr = (sort_right_rec1 - sort_left_rec1 > 0).nonzero()[0]
    rec1_adr = (sort_right_rec2 - sort_left_rec2 > 0).nonzero()[0]

    return rec1a[rec1_adr],rec2a[rec2_adr]

#Makes a PMF hash array. Useful outside the rec_match context. Necessary for
#rec_match and rec_match_srt.
def mk_hsh_arr(rec_arr):
    tmark.tm('Creating Hash Array')
    num_rec = len(rec_arr)
    rec_hsh = np.chararray(num_rec,itemsize=16,unicode=True)
    for i in range(num_rec):
        pt = rec_arr['PLATE'][i]
        mt = rec_arr['MJD'][i]
        ft = rec_arr['FIBERID'][i]
        hst = '{0:05d}-{1:05d}-{2:04d}'.format(pt,mt,ft)
        rec_hsh[i] = hst
        pb.pbar(i,num_rec)

    return rec_hsh
