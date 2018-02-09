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
#For ease, use the following command to import to your program:
#  import fits_writer_error_trap as fet
from astropy.io import fits
import glob

def nm_up(in_rec,in_name):
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
