from astropy.io import fits
import numpy as np
from pydl.pydlutils.sdss import sdss_flagname
from pydl.pydlutils.sdss import sdss_flagval
import tmark
import progressBar as pb
import glob
from astropy import units as u
from astropy.coordinates import SkyCoord
from sciCon import mks

#Makes either a PMF hash array or a PM hash array. Parameter include_fiber
#controls this. Returns a 1-D array of strings.
def mk_hash(inrec,include_fiber=True,include_mjd=True,mname='MJD',fibname='FIBERID'):
    tmark.tm('Creating Hash Array')
    num_rec = len(inrec)
    len_plate = 5
    if ((include_fiber==True)&(include_mjd==True)):
        str_size = len_plate + 11
        rec_hsh = np.chararray(num_rec,itemsize=str_size,unicode=True)
        rec_hsh = np.array(['%05d-%05d-%04d'%(pt,mt,ft) for pt,mt,ft in zip(
                            inrec['PLATE'],inrec[mname],inrec[fibname])])
    elif ((include_fiber==False)&(include_mjd==True)):
        str_size = len_plate + 6
        rec_hsh = np.chararray(num_rec,itemsize=str_size,unicode=True)
        rec_hsh = np.array(['%05d-%05d'%(pt,mt) for pt,mt in zip(
                            inrec['PLATE'],inrec[mname])])
    elif ((include_fiber==True)&(include_mjd==False)):
        str_size = len_plate + 5
        rec_hsh = np.chararray(num_rec,itemsize=str_size,unicode=True)
        rec_hsh = np.array(['%05d-%04d'%(pt,ft) for pt,ft in zip(
                            inrec['PLATE'],inrec[fibname])])
    return rec_hsh

#This function is designed to load a FITS table into a numpy structured array,
#preserving the column names and data types (case-wise).
def file_load(infile):
    cat_file = fits.open(infile)[1].data
    num_rec = len(cat_file)
    column_names = np.array(cat_file.columns.names)
    column_forms = np.array(cat_file.columns.formats,dtype='U10')
    num_cols = len(column_names)

    for i in range(num_cols):
        if 'A' in column_forms[i]:
            column_forms[i] = 'U{}'.format(column_forms[i][0:-1])
        elif 'B' in column_forms[i]:
            num_check = column_forms[i][0:-1]
            if num_check == '':
                column_forms[i] = 'B'
            else:
                column_forms[i] = '{}B'.format(num_check)
        elif 'D' in column_forms[i]:
            num_check = column_forms[i][0:-1]
            if num_check == '':
                column_forms[i] = 'float64'
            else:
                column_forms[i] = '{}float64'.format(num_check)
        elif 'E' in column_forms[i]:
            num_check = column_forms[i][0:-1]
            if num_check == '':
                column_forms[i] = 'float32'
            else:
                column_forms[i] = '{}float32'.format(num_check)
        elif 'K' in column_forms[i]:
            num_check = column_forms[i][0:-1]
            if num_check == '':
                column_forms[i] = 'int64'
            else:
                column_forms[i] = '{}int64'.format(num_check)
        elif 'J' in column_forms[i]:
            num_check = column_forms[i][0:-1]
            if num_check == '':
                column_forms[i] = 'int32'
            else:
                column_forms[i] = '{}int32'.format(num_check)
        elif 'I' in column_forms[i]:
            num_check = column_forms[i][0:-1]
            if num_check == '':
                column_forms[i] = 'int16'
            else:
                column_forms[i] = '{}int16'.format(num_check)
    #Load the drfile into a numpy structured array
    struct_temp = np.zeros(num_rec,dtype={'names':column_names,'formats':column_forms})
    for cname in column_names:
        struct_temp[cname] = cat_file[cname]

    return struct_temp

#This program is designed to write out a fits file even if a file already exists
#with that name in the output destination. It will attempt to write the file with
#the given file name first. If that fails, it will append an _00 to the file name
#before the .fits and attempt to write the file. It will iterate through numbers
#until it finds a two-digit appendix that works.

#Call the function in your program AFTER you have created the BinTableHDU. The name
#of that HDU object is the in_rec. The in_name is the name of the file you want
#without the .fits at the end.
#This essentially replaces the last xxx.writeto(file_name) that is used to write
#fits files out.
def fet(in_rec,in_name,quiet=False):
    #rename the BinTableHDU as dof for internal purposes.
    dof = in_rec
    #Glob is used to find the list of files that exists. This is only for attempting
    #to write the file exactly as the name given.
    check_str = '{}*.fits'.format(in_name)
    nm_list = glob.glob(check_str)
    #Version number initialization. Used if the original file name already exists.
    ver_num = 0
    #If the file doesn't already exist, write it out.
    if len(nm_list) == 0:
        out_name = '{}.fits'.format(in_name)
        dof.writeto(out_name)
    #If the original desired file name exists, use try syntax to attempt to write
    #out the file with _00. Failing that (IOError), iterate the version number
    #and then try again.
    else:
        while True:
            try:
                out_name = '{}_{:02d}.fits'.format(in_name,ver_num)
                dof.writeto(out_name)
                break
            except IOError:
                ver_num += 1

    #Tell me the final file name so I have a record of what was output.
    if quiet==False:
        print('\nFile Written out as: {}'.format(out_name))

    return out_name

#This program is designed to find all of the matches in two different fits
#flat files, finding the addresses where the objects in 1 match the objects in 2.
#rec_match_srt only works well with hash arrays (they compare one array)
#so you will need to use mk_hsh_arr first to make PMF hash arrays for each fits
#flat file.
def rec_match_srt(rec1,rec2,verbose=True):
    rec1a = np.argsort(rec1)
    rec2a = np.argsort(rec2)
    if verbose==True:
        tmark.tm('Starting Searchsorted')

    sort_left_rec1 = rec1[rec1a].searchsorted(rec2[rec2a],side='left')
    sort_right_rec1 = rec1[rec1a].searchsorted(rec2[rec2a],side='right')
    sort_left_rec2 = rec2[rec2a].searchsorted(rec1[rec1a],side='left')
    sort_right_rec2 = rec2[rec2a].searchsorted(rec1[rec1a],side='right')

    rec2_adr = (sort_right_rec1 - sort_left_rec1 > 0).nonzero()[0]
    rec1_adr = (sort_right_rec2 - sort_left_rec2 > 0).nonzero()[0]

    return rec1a[rec1_adr],rec2a[rec2_adr]

#This program finds the object that you want in the file, then tells you all of
#the targeting flags that object has, by name.
def flags(infile,plate_in,fiber_in):
    pt = plate_in
    ft = fiber_in
    w1 = np.where((infile['PLATE']==pt)&(infile['FIBERID']==ft))[0]
    unstr = '--UNUSED--'

    try:
        bt1_flags = sdss_flagname('BOSS_TARGET1',infile['BOSS_TARGET1'][w1])
        if infile['BOSS_TARGET1'][w1] == -1:
            bt1_flags = unstr
    except (ValueError, KeyError):
        bt1_flags = unstr
    try:
        et0_flags = sdss_flagname('EBOSS_TARGET0',infile['EBOSS_TARGET0'][w1])
        if infile['EBOSS_TARGET0'][w1] == -1:
            et0_flags = unstr
    except (ValueError, KeyError):
        et0_flags = unstr
    try:
        et1_flags = sdss_flagname('EBOSS_TARGET1',infile['EBOSS_TARGET1'][w1])
        if infile['EBOSS_TARGET1'][w1] == -1:
            et1_flags = unstr
    except (ValueError, KeyError):
        et1_flags = unstr
    try:
        et2_flags = sdss_flagname('EBOSS_TARGET2',infile['EBOSS_TARGET2'][w1])
        if infile['EBOSS_TARGET2'][w1] == -1:
            et2_flags = unstr
    except (ValueError, KeyError):
        et2_flags = unstr
    try:
        at1_flags = sdss_flagname('ANCILLARY_TARGET1',infile['ANCILLARY_TARGET1'][w1])
        if infile['ANCILLARY_TARGET1'][w1] == -1:
            at1_flags = unstr
    except (ValueError, KeyError):
        at1_flags = unstr
    try:
        at2_flags = sdss_flagname('ANCILLARY_TARGET2',infile['ANCILLARY_TARGET2'][w1])
        if infile['ANCILLARY_TARGET2'][w1] == -1:
            at2_flags = unstr
    except (ValueError, KeyError):
        at2_flags = unstr

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

#This is for printing out a more human-readable version of a single row in a
#FITS record. The rec_adr needs to be used if input masked array has more than
#one "row" of data. Defaults to 0, same as calling inrec[w][0] or inrec[w[0]].
def rec_view(inrec,rec_adr=0):
    cnames = np.array(inrec.dtype.names)
    max_char = 0
    for i in range(len(cnames)):
        temp_length = len(cnames[i])
        if temp_length > max_char:
            max_char = temp_length
    for i in range(len(cnames)):
        out_str = '{0:{1}} | {2}'.format(cnames[i],max_char,inrec[cnames[i]][rec_adr])
        print(out_str)
    print('\n')

#This is for getting the index for a record matching a plate,mjd,fiber.
#I use this functionality enough, I was tired of typing the where command.
def find_rec(inrec,plate,mjd,fiberid,pname='PLATE',mname='MJD',fibname='FIBERID'):
    pt,mt,ft = plate,mjd,fiberid
    wpmf = np.where((inrec[pname]==pt)&(inrec[mname]==mt)&(inrec[fibname]==ft))[0]
    return wpmf

#The next two are for converting from sexagessimal to decimal RA/DEC or back.
#Returns the desired in human-readable format.
def dec2sex(rad,decd):
    c = SkyCoord(ra=rad*u.degree,dec=decd*u.degree,frame='icrs')
    c_hms = c.ra.hms
    c_dms = c.dec.dms
    sig = '+'
    if decd < 0:
        sig = '-'
    RAstr = 'RA - {0:02d}:{1:02d}:{2:05.2f}'.format(int(c_hms.h),int(c_hms.m),c_hms.s)
    DECstr = 'DEC - {0}{1:02d}:{2:02d}:{3:04.1f}'.format(sig,int(c_dms.d),int(c_dms.m),c_dms.s)
    out_str = '{0} | {1}'.format(RAstr,DECstr)
    print(out_str)

def sex2dec(ras,decs):
    c = SkyCoord(ras,decs,frame='icrs')
    c_ra = c.ra.degree
    c_dec = c.dec.degree
    out_str = 'RA - {0:.4f} | DEC - {1:.4f}'.format(c_ra,c_dec)
    print(out_str)

#The next two are for finding redshift error in km/s. p controls whether it
#prints the delta_v or returns the value (for importing).
def zdiff(z1,zt,p=0):
    ckms = mks.c / 1000
    zd = (ckms*np.abs(z1-zt))/(1+zt)
    if p!=0:
        print('dv = {:4d} km/s'.format(int(zd)))
    else:
        return zd

#This calculates the redshift error based on the wavelength shown as the line
#center (lam1) and what actually corresponds to the line center (lamt). For
#testing whether a misidentified line center corresponds to a delta_v > 3000 km/s.
#The p paramter controls printing vs. returning like above.
def lzdiff(lam1,lamt,p=0):
    ckms = mks.c / 1000
    ld = (ckms*np.abs(lam1-lamt))/lamt
    if p!=0:
        print('dv = {:4d} km/s'.format(int(ld)))
    else:
        return ld

def cat_set(cat1,cat2):
    cathash1,cathash2 = mk_hash(cat1),mk_hash(cat2)
    cat1args,cat2args = rec_match_srt(cathash1,cathash2)
    return cat1args,cat2args
