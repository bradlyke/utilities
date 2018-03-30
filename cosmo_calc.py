"""
A cosmological calculator that can find various ages and times based on a given
redshift. Can also calculate and absolute i-band magnitude based on an apparent
magnitude, redshift, and galactic extinction value. Requires the k-correction
table in the same folder: k_corr_tab.dat.

use in another program via
import cosmo_calc as csm

Parameters
----------
z_in : :class:'float'
	The redshift of the object you want to calculate stuff for.
amag : :class:'float'
	OPTIONAL - the apparent magnitude (i-band) of an object.
    (used for mag_test)
a_i : :class:'float'
	OPTIONAL - the i-band galactic extinction in magnitudes for an object.
    (used for mag_test)

Originally designed for a terminal 96 characters wide.

Returns
----------
--red_report--
:class:'str'
	Outputs a report similar to Ned Wright's catalogue for a given redshift.

--mag_test--
:class:'float'
    Returns the absolute i-band magnitude of the object.
"""
import numpy as np

#USES FLAT COSMOLOGY WITH COSMOLOGICAL CONSTANT: BENCHMARK MODEL
#Define a number of used constants, taken 20180314 from Ned Wright's page.
#These are my values for DR15Q

h0 = 69.6 #Hubble Constant km/s/Mpc
omegaM = 0.286 #Matter density parameter
#omegaL = 0.714 #Vacuum density parameter defined by Ned Wright.
omegaR = 9.0e-5 #Taken from Barbara Ryden's Introduction to Cosmology, 2ed.
omegaL = 1.0 - omegaR - omegaM #Vacuum density parameter defined using omegaR
c = 299792.458 #speed of light in km/s

#Calculates the f(z) for finding proper distance using integral form.
def z_func(z_t):
    b1 = omegaR * (1 + z_t)**4 #Radiation term
    b2 = omegaM * (1 + z_t)**3 #Matter term
    b3 = omegaL              #Cosmological Constant (vac) term
    denom = h0 * (b1 + b2 + b3)**(0.5) #denominator of proper distance
    dP_temp = c  * denom**(-1.0) #differential proper distance for current z_now
    return dP_temp

#This version of the comoving distance is calculated using Simpson's Rule.
#During testing it was found that trapezoid and simpson's method match to 6
#significant figures, but Simpon's Method is O(h^4), while trapezoid is O(h^2).
def comove_dist(z_p):
    num_steps = 100000 #number of steps from 0 to z_p.
    dz = z_p / num_steps #step size.
    p_temp = 0 #initlized the placeholder for differential comoving distance
    pa = z_func(0) #Evaluate the start point
    pb = z_func(z_p) #Evaluate the end point
    dp = pa + pb #Add these together to initialize d
    #Adds the middle terms in the Simpson's method.
    for i in range(1,num_steps):
        z_0 = i * dz #current z_i
        #If we are on an even value for i, use 2*f(z)
        if (i % 2 == 0):
            p_temp = 2 * z_func(z_0)
        #If we are on an odd value for i, use 4*f(z)
        else:
            p_temp = 4 * z_func(z_0)
        dp += p_temp#Cumulative sum
    dPS = (dz / 3.0) * dp #Last step of Simpson's Method.
    return dPS #Returns the comoving distance in Mpc.

def comove_vol(z_v):
    rad = comove_dist(z_v)
    vol = ((4.0/3) * np.pi * rad**(3.0))/(1.0e9)
    return vol

#Lighttravel time. If you want to find the current age of the universe, use a z=1090
#z_l is the lower boundary. If you want age of the universe at a specific redshift
#z_l is that redshift, z_p is 1089.
#If you want to find light travel time to a redshift, z_l is 0, z_p is that redshift.
def lookback_time(z_l,z_p):
    num_steps = 1000000 #number of steps from 0 to z_p.
    dz = z_p / num_steps #step size.
    t_temp = 0 #initlized the placeholder for differential comoving distance
    ta = (1+z_l)**(-1.0) * z_func(z_l) #Evaluate the start point
    tb = (1+z_p)**(-1.0) * z_func(z_p) #Evaluate the end point
    dt = ta + tb #Add these together to initialize d
    #Adds the middle terms in the Simpson's method.
    for i in range(1,num_steps):
        z_0 = z_l + (i * dz) #current z_i
        #If we are on an even value for i, use 2*f(z)
        if (i % 2 == 0):
            t_temp = 2 * (1+z_0)**(-1.0) * z_func(z_0)
        #If we are on an odd value for i, use 4*f(z)
        else:
            t_temp = 4 * (1+z_0)**(-1.0) * z_func(z_0)
        dt += t_temp#Cumulative sum
    lT = (dz / 3.0) * dt #Last step of Simpson's Method.
    lT = lT * 978.3386176 / c #Corrects for h0
    return lT #Returns the lookback time in Gyr.

#Calculate the Luminosity distance from comoving distance and redshift.
def lum_dist(dp,z_l):
    #dp must in input as Mpc
    return dp * (1+z_l) #Will output in Mpc

#Calculate the Angular distance from the comoving distance and redshift.
def ang_dist(dp,z_l):
    #dp must in input as Mpc
    return dp * (1+z_l)**(-1.0) #Will output in Mpc

#Calculate the distance modulus using luminosity distance.
def dist_mod(dl):
    #dL must be input as Mpc
    return 5*np.log10(dl) + 25 #unitless

#Table is for i-band only and is loaded from Richards et. al. (2006)
def k_corr(z_k):
    k_table = np.loadtxt('k_corr_tab.dat',dtype='f4',delimiter=',') #Load k-correction table
    diff_arr = np.absolute(z_k-k_table[:,0]) #Find the difference between our redshift and table.
    min_diff = np.amin(diff_arr) #Find the minimum difference
    wmin = np.where(diff_arr == min_diff)[0] #Find the address of the minimum distance.
    kC = k_table[wmin[0],1] #Find the k-correction value for that redshift (usually 2.0).
    return kC #Return the k-correction, unitless (it's a magnitude correction).

#Find the absolute magnitude given the previously calculated values and object
#apparent magnitude and galactic extinction.
#Only good for i-band.
def abs_Mag(apmag,ext,dm,kC):
    #apmag is the apparent magnitude, ext is galactic extinction
    #dm is distance modulus from dist_mod, and kC is k correction from k_corr.
    return apmag - ext - dm - kC #unitless

#For getting an absolute magnitude from the program. Example for a QSO with
#redshift of 2.355, apparent magnitude of 20.56, and extinction of 0.045:
#mag_test(2.355,20.56,0.045)
def mag_test(z_test,amag,a_i):
    Dp = comove_dist(z_test)
    Dl = lum_dist(Dp,z_test)
    Dm = dist_mod(Dl)
    Kc = k_corr(2.0)
    Mag = abs_Mag(amag,a_i,Dm,Kc)
    tstr1 = ' '*3 +'Z'+' '*3+'|'+' '*4+'Dp'+' '*6+'|'+' '*5+'Dl'+' '*5+'|'+' '*4+'Mag'
    brkr = '-'*44
    tstr = '{0:6.4f} | {1:6.1f} Mpc | {2:7.1f} Mpc| {3:9.5f}'.format(z_test,Dp,Dl,Mag)
    print(tstr1)
    print(brkr)
    print(tstr)

#This will generated a report on the universe and calculated distances/ages
#for a given redshift, then output to the terminal. Also reports the redshift
#given, density parameters used, and the Hubble constant.
def red_report(z_in):
    uni_age = np.round(lookback_time(0,1089),1) #Gyr
    age_z = np.round(lookback_time(z_in,1089),1) #Gyr
    age_light = np.round(lookback_time(0,z_in),1) #Gyr
    coradial = np.round(comove_dist(z_in),1) #In Mpc
    covol = np.round(comove_vol(z_in),2) #In Gpc^3
    size_ang = np.round(ang_dist(coradial,z_in),1) #In Mpc
    lumdist = np.round(lum_dist(coradial,z_in),1)
    print('\n')
    print('Redshift: {}, Wm: {}, Wr: {}, Wl: {:0.3f}, H0: {}'.format(z_in,omegaM,omegaR,omegaL,h0))
    print('---------------------------------------------------------')
    print('Age of Universe:   {:7.1f} Gyr'.format(uni_age))
    print('Age at Redshift:   {:7.1f} Gyr'.format(age_z))
    print('Light Travel Time: {:7.1f} Gyr'.format(age_light))
    print('Comoving Distance: {:7.1f} Mpc'.format(coradial))
    print('Comoving Volume:   {:7.2f} Gpc^3'.format(covol))
    print('Angular Size Dist: {:7.1f} Mpc'.format(size_ang))
    print('Luminosity Dist:   {:7.1f} Mpc'.format(lumdist))
