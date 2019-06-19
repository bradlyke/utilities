from astropy.io import fits
import numpy as np
import fits_writer_error_trap as fet
import progressBar as pb
import tmark
import sys

#USES FLAT COSMOLOGY WITH COSMOLOGICAL CONSTANT: BENCHMARK MODEL
#Define a number of used constants, taken 20180314 from Ned Wright's page.
#These are my values for DR15Q

h0 = 69.6 #km/s/Mpc
omegaM = 0.286
#omegaL = 0.714
omegaR = 9.0e-5
omegaL = 1.0 - omegaR - omegaM
'''
#These are the values Isabelle used in DR14Q
h0 = 67.8 #Hubble Constant
omegaM = 0.308 #Matter density parameter
omegaR = 9.0e-5 #radiation density parameter
omegaL = 0.692 #Vacuum density parameter
'''
#Common constants for both me and Isabelle.
#alphaV = -0.5
c = 299792.458 #speed of light in km/s

#Calculates the f(z) for finding proper distance using integral form.
def z_func(z_t):
    b1 = omegaR * (1 + z_t)**4 #Radiation term
    b2 = omegaM * (1 + z_t)**3 #Matter term
    b3 = omegaL              #Cosmological Constant (vac) term
    denom = h0 * (b1 + b2 + b3)**(0.5) #denominator of proper distance
    dP_temp = c  * denom**(-1.0) #differential proper distance for current z_now
    return dP_temp

#This version of the proper distance uses a simple Reimann sum with small steps.
def prop_dist_reim(z_p):
    num_steps = 100000 #10^6 is too many, 10^3 is too few.
    #Using num_steps = 10000 and omegaR = 9.0e-5 makes dP agree with Ned Wright's
    #values to 6 significant figures (dP and dL actually).
    dz = z_p / num_steps #step size
    z_now = 0 #initialized redshift
    dP = 0 #initialized proper distance
    #Calculate the proper distance using simple Reimann sum with very small
    #steps. dP = sum((c * dz)/(h0*(Or(1+zi)^4 + Om(1+zi)^3 + Ol)^0.5)) over zi
    for i in range(num_steps):
        z_now = i * dz #current redshift value
        dP_temp = z_func(z_now) #differential proper distance for current z_now
        dP += dP_temp #running total for proper distance
    dP = dz * dP
    return dP

#This version of the proper distance is calculated using the trapezoid rule.
def prop_dist_trap(z_p):
    num_steps = 100000 #number of steps from 0 to z_p.
    dz = z_p / num_steps #step size.
    p_now = 0 #initlized the placeholder for partial proper distance
    dP0T = z_func(0) #Evaluate the start point
    dPNT = z_func(z_p) #Evaluate the end point
    dpT = dP0T + dPNT #Add these together to initialize dP
    #Adds middle terms in trapezoid method: Sum_i (2 * f(z_i))
    for i in range(1,num_steps):
        z_0 = i * dz #Current z_i
        p_now = 2 * z_func(z_0) #differential proper distance
        dpT = dpT + p_now #cumulative sum
    dPT = (dz / 2.0) * dpT #Last step in trapezoid method.
    return dPT #Return the proper distance

#This version of the proper distance is calculated using Simpson's Rule.
#During testing it was found that trapezoid and simpson's method match to 6
#significant figures.
def prop_dist_simp(z_p):
    num_steps = 100000 #number of steps from 0 to z_p.
    dz = z_p / num_steps #step size.
    z_now = 0 #initlized the placeholder for partial proper distance
    dP0S = z_func(0) #Evaluate the start point
    dPNS = z_func(z_p) #Evaluate the end point
    dpS = dP0S + dPNS #Add these together to initialize dP
    #Adds the middle terms in the Simpson's method.
    for i in range(1,num_steps):
        z_0 = i * dz #current z_i
        #If we are on an even value for i, use 2*f(z)
        if (i % 2 == 0):
            p_now = 2 * z_func(z_0)
        #If we are on an odd value for i, use 4*f(z)
        else:
            p_now = 4 * z_func(z_0)
        dpS = dpS + p_now#Cumulative sum
    dPS = (dz / 3.0) * dpS #Last step of Simpson's Method.
    return dPS

#Calculate the Luminosity distance from proper distance and redshift.
def lum_dist(dp,z_l):
    #dp must in input as Mpc
    return dp * (1+z_l) #Will output in Mpc

#Calculate the distance modulus using luminosity distance.
def dist_mod(dl):
    #dL must be input as Mpc
    return 5*np.log10(dl) + 25 #unitless

#K-Correction for redshift changing slope of continuum due to cosmological reddening.
#def k_corr(z_k):
    #alphaV is defined above.
    #return -2.5*(1+alphaV)*np.log10(1+z_k) #unitless

def k_corr(z_k):
    k_table = np.loadtxt('k_corr_tab.dat',dtype='f4',delimiter=',') #Load k-correction table
    diff_arr = np.absolute(z_k - k_table[:,0])
    min_diff = np.amin(diff_arr)
    wmin = np.where(diff_arr == min_diff)[0]
    kC = k_table[wmin[0],1]
    return kC

#Find the absolute magnitude given the previously calculated values and data
#taken from the fits file.
def abs_Mag(apmag,ext,dm,kC):
    #apmag is the apparent magnitude, ext is galactic extinction
    #dm is distance modulus from dist_mod, and kC is k correction from k_corr.
    return apmag - ext - dm - kC #unitless

#Primary program for reading the FITS file and calculating redshift for all records.
#ifile is the input fits file. norml is for normalizing all of the records to
#redshift of 2.0. Not really used.
def get_mags(ifile,norml=False):
    #norml means normalized all to the same redshift (z=2) if true.
    #If norml is false, uses the best redshift from the catalogue ('Z').
    drfile = fits.open(ifile)[1].data #Load the FITS file.
    num_rec = len(drfile) #Find the number of records.
    magarr = np.zeros(num_rec,dtype='f4') #Create an array to hold absolute magnitude
    magarr[:] = 999 #replace all of them with 999. If 999 appears afterwards I missed that record.
    tmark.tm('Starting Absolute Magnitude Calculations')
    #Iterate through all objects.
    for i in range(num_rec):
        #Define the redshift for the object
        rec_z = drfile['Z'][i] #Current working redshift.
        if rec_z < 0:
            magarr[i] = 0 #This catches blazars or other objects with a bad Z.
            continue
        if (norml==False):
            z_temp = rec_z
        else:
            z_temp = 2.0
        appMag = drfile['PSFMAG'][i,3] #apparent PSF magnitude in i-band
        A_i = drfile['EXTINCTION'][i,3] #galactic extinction in magnitudes (i-band)
        Dp = prop_dist_simp(z_temp) #find the proper distance using simpson's method
        Dl = lum_dist(Dp,z_temp) #find the luminosity distance
        Dm = dist_mod(Dl)
        #Kc = k_corr(z_temp) #This is for K-correction using record Z.
        Kc = k_corr(2.0) #This is the standard way of finding K-correction.
        absMag = abs_Mag(appMag,A_i,Dm,Kc) #get the Absolute magnitude in i-band.
        magarr[i] = absMag #Save it to placeholder array

        pb.pbar(i,num_rec)
        #The following two lines are for testing purposes. Comment if you don't want
        #to print interim values.
        #test_str = 'Z: {0:6.4f} | Dp: {1:6.1f} | Dl: {2:7.1f} | Mag: {3:9.5f}'.format(z_temp,Dp,Dl,absMag)
        #print(test_str)

    #User feedback for when the program completes.
    print('\n')
    print('Absolute Magnitudes Calculated')
    #The following writes out the new fits record with the absolute magnitude column.

    #Need to place the new column in the right place, between PSFMAGERR and EXTINCTION
    colnames = np.array(drfile.columns.names) #Make an array of all current column names.
    wpsf = np.where(colnames=='psfmagerr')[0] #Find column address of PSFMAGERR
    psfmagerr_adr = wpsf[0] + 1 #Move one beyond that. If you don't have this, you remove PSFMAGERR column
    data_col1 = drfile.columns[0:psfmagerr_adr] #Old original columns.
    data_col2 = drfile.columns[psfmagerr_adr:] #Columns to shift after absMag column.
    magCol = fits.ColDefs([fits.Column(name='m_i',format='E',array=magarr)]) #New column definition.

    #Put together the HDU and write out the file using my fits_writer_error_trap function.
    data_out = fits.BinTableHDU.from_columns(data_col1 + magCol + data_col2)
    ofile = '../../data/newq/DRNQ_v1_8_wisemags'
    fet.nm_up(data_out,ofile)

get_mags(sys.argv[1]) #Used so I can run this from command line.
