################################################################################
#                                                                              #
#  SCRIPT: SDSS_NAME Generator                                                 #
#  AUTHOR: Brad Lyke                                                           #
#                                                                              #
#  PURPOSE: This script will generate an SDSS_NAME, without the J from a given #
#           RA/DEC that are in decimal degree form.                            #
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

def s_name(ra,dec):
    #Convert the Decimal RA/DEC to Sexegessimal.
    c = sc(ra = ra*u.degree, dec = dec*u.degree)
    c_hms = c.ra.hms
    c_dms = c.dec.dms

    #To get the correct sign in the SDSS_NAME check the leading sign on the
    #sexegessimal dec.
    if c_dms[0] > 0:
        sig = '+'
    elif c_dms[0] < 0:
        sig = '-'
    elif c_dms[0] == 0:
        if c_dms[1] > 0:
            sig = '+'
        elif c_dms[1] < 0:
            sig = '-'
        elif c_dms[1] == 0:
            if c_dms[2] > 0:
                sig = '+'
            elif c_dms[2] < 0:
                sig = '-'
            else:
                sig = '+'

    #Need to truncate the decimal places in seconds without rounding. Converting
    #to a string and then cutting string characters does this.
    c_htemp = '{0:07.4f}'.format(abs(c_hms[2][0]))
    c_dtemp = '{0:07.4f}'.format(abs(c_dms[2][0]))
    #Cut down the number of characters in the string. The :-2 means all but the
    #last two.
    c_hsec = c_htemp[:-2]
    c_dsec = c_dtemp[:-3]
    #Create the string for the RA part.
    snr = '{0:02d}{1:02d}{2}'.format(abs(int(c_hms[0])),abs(int(c_hms[1])),c_hsec)
    #Create the string for the DEC part.
    snd = '{0:02d}{1:02d}{2}'.format(abs(int(c_dms[0])),abs(int(c_dms[1])),c_dsec)

    #Create the total string, with sign, and without J for updating fields. If
    #a J is needed, then add it in the parent function calling this one.
    sntest = '{0}{1}{2}'.format(snr,sig,snd)

    #Output the final sdss_name.
    return sntest
