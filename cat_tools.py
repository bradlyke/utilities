from astropy.io import fits
import numpy as np
import tmark
import progressBar as pb
import glob

#Makes either a PMF hash array or a PM hash array. Parameter include_fiber
#controls this. Returns a 1-D array of strings.
def mk_hash(input_rec,include_fiber=True):
    tmark.tm('Creating Hash Array')
    num_rec = len(input_rec)
    max_plate = str(np.amax(input_rec['PLATE']))
    len_plate = len(max_plate)
    if include_fiber==True:
        str_size = len_plate + 11
        rec_hsh = np.chararray(num_rec,itemsize=str_size,unicode=True)
    else:
        str_size = len_plate + 6
        rec_hsh = np.chararray(num_rec,itemsize=str_size,unicode=True)
    for i in range(num_rec):
        pt = input_rec['PLATE'][i]
        mt = input_rec['MJD'][i]
        if include_fiber == True:
            ft = input_rec['FIBERID'][i]
            hst = '{0:0{1}d}-{2:05d}-{3:04d}'.format(pt,len_plate,mt,ft)
        else:
            hst = '{0:0{1}d}-{2:05d}'.format(pt,len_plate,mt)
        rec_hsh[i] = hst
        pb.pbar(i,num_rec)

    return rec_hsh

#This function is designed to load a FITS table into a numpy structured array,
#preserving the column names and data types (case-wise).
def file_load(infile):
    cat_file = fits.open(infile)[1].data
    num_rec = len(cat_file)
    column_names = np.array(cat_file.columns.names)
    column_forms = np.array(cat_file.columns.formats,dtype='U9')
    num_cols = len(column_names)

    for i in range(num_cols):
        if 'A' in column_forms[i]:
            column_forms[i] = 'U{}'.format(column_forms[i][0:-1])
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
def fet(in_rec,in_name):
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
    print('\nFile Written out as: {}'.format(out_name))

#This program is designed to find all of the matches in two different fits
#flat files, finding the addresses where the objects in 1 match the objects in 2.
#rec_match_srt only works well with hash arrays (they compare one array)
#so you will need to use mk_hsh_arr first to make PMF hash arrays for each fits
#flat file.
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
