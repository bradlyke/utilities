################################################################################
#                                                                              #
#  SCRIPT: SDSS_NAME Generator                                                 #
#  AUTHOR: Brad Lyke                                                           #
#                                                                              #
#  PURPOSE: This script will generate an SDSS_NAME, without the J from a given #
#           RA/DEC that are in decimal degree form. Designed for array of      #
#           RAs and DECs.                                                      #
#                                                                              #
#  VARIABLES: ra  - the RA of the object in decimal form.                      #
#             dec - the DEC of the object in decimal form.                     #
#                                                                              #
#  INPUT: Designed to be imported into another function and called.            #
#         Import with the following:                                           #
#           import sdss_name as sdn                                            #
#         Use the funtion with the following:                                  #
#           <name_variable> = sdn.s_name(ra_in,dec_in)                         #
#  OUTPUT: Returns just the SDSS_NAME with the correct sign and no J. Truncates#
#          decimal places in seconds without rounding.                         #
#                                                                              #
# Note: REQUIRES PYTHON 3.5 or later.                                          #
################################################################################

#These are the only two packages needed for converting to sexegessimal coords.
from astropy.coordinates import SkyCoord as sc
from astropy import units as u
import numpy as np

def s_name(rat,det):
    #Find the number of records
    num_obj = len(rat)
    #Convert the Decimal RA/DEC to Sexegessimal.
    c = sc(ra = rat*u.degree, dec = det*u.degree)
    c_hms = c.ra.hms
    c_dms = c.dec.dms
    #Placeholder arrays for the sign in the name, and for the names themselves.
    sig = np.chararray(num_obj,itemsize=1,unicode=True)
    sname_array = np.chararray(num_obj,itemsize=18,unicode=True)
    #To get the correct sign in the SDSS_NAME check the leading sign on the
    #sexegessimal dec.
    #Check the degree sign--
    wdp = np.where(c_dms.d > 0.0)[0]
    sig[wdp] = '+'
    wdm = np.where(c_dms.d < 0.0)[0]
    sig[wdm] = '-'
    #If the degree is 0, check the minutes sign
    wd0 = np.where(c_dms.d == 0.0)[0]
    wmp = np.where(c_dms.m[wd0] > 0.0)[0]
    sig[wd0[wmp]] = '+'
    wmn = np.where(c_dms.m[wd0] < 0.0)[0]
    sig[wd0[wmn]] = '-'
    #If both the degree and minutes are 0, check the seconds sign. If seconds
    #are 0, then assign a +. This should only happen if DEC is 00:00:00.0
    wm0 = np.where(c_dms.m[wd0] == 0.0)[0]
    wsp = np.where(c_dms.s[wd0[wm0]] >= 0.0)[0]
    sig[wd0[wm0[wsp]]] = '+'
    wsm = np.where(c_dms.s[wd0[wm0]] < 0.0)[0]
    sig[wd0[wm0[wsm]]] = '-'

    #Need to truncate the decimal places in seconds without rounding.
    #Converting to a string and then cutting string characters does this.
    #First convert the seconds to a string for both RA and DEC.
    c_htemp = np.array(['%07.4f'%X for X in abs(c_hms.s)])
    c_dtemp = np.array(['%07.4f'%Y for Y in abs(c_dms.s)])
    #Cut down the number of characters in the string. The :-2 means all but the
    #last two. For an SDSS name, the RA seconds have 2 decimal places, the DEC
    #seconds have 1 decimal place.
    c_hsec = np.array([obj[:-2] for obj in c_htemp])
    c_dsec = np.array([obj[:-3] for obj in c_dtemp])
    #Create the string for the RA part.
    snr = np.array(['%02d%02d%s'%(Hh,mm,ss) for Hh,mm,ss in zip(abs(c_hms.h),abs(c_hms.m),c_hsec)])
    #Create the string for the DEC part.
    snd = np.array(['%02d%02d%s'%(Dd,mm,ss) for Dd,mm,ss in zip(abs(c_dms.d),abs(c_dms.m),c_dsec)])

    #Create the total string, with sign, and without J for updating fields. If
    #a J is needed, then add it in the parent function calling this one.
    sname_array = np.array(['%s%s%s'%(rs,sigs,ds) for rs,sigs,ds in zip(snr,sig,snd)])

    #Output the final SDSS name.
    return sname_array
