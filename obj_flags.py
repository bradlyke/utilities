from astropy.io import fits
from pydl.pydlutils.sdss import sdss_flagname
from pydl.pydlutils.sdss import sdss_flagval
import numpy as np

def flags(infile,plate_in,fiber_in):
    pt = plate_in
    ft = fiber_in
    w1 = np.where((infile['PLATE']==pt)&(infile['FIBERID']==ft))[0]

    bt1_flags = sdss_flagname('BOSS_TARGET1',infile['BOSS_TARGET1'][w1])
    et0_flags = sdss_flagname('EBOSS_TARGET0',infile['EBOSS_TARGET0'][w1])
    et1_flags = sdss_flagname('EBOSS_TARGET1',infile['EBOSS_TARGET1'][w1])
    et2_flags = sdss_flagname('EBOSS_TARGET2',infile['EBOSS_TARGET2'][w1])
    at1_flags = sdss_flagname('ANCILLARY_TARGET1',infile['ANCILLARY_TARGET1'][w1])
    at2_flags = sdss_flagname('ANCILLARY_TARGET2',infile['ANCILLARY_TARGET2'][w1])

    print('\n')
    print('Object Flags')
    print('------------')
    print('BOSS_TARGET1: {}'.format(bt1_flags))
    print('EBOSS_TARGET0: {}'.format(et0_flags))
    print('EBOSS_TARGET1: {}'.format(et1_flags))
    print('EBOSS_TARGET2: {}'.format(et2_flags))
    print('ANCILLARY_TARGET1: {}'.format(at1_flags))
    print('ANCILLARY_TARGET2: {}'.format(at2_flags))
    print('\n')
