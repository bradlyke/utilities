################################################################################
#                                                                              #
#  SCRIPT: Visual Inspection Calculator                                        #
#  AUTHOR: Brad Lyke                                                           #
#                                                                              #
#  PURPOSE: This script has 3 subroutines useful while completing visual       #
#           inspections. The program has a 'active menu,' so it will not exit  #
#           until it is told to or encounters bad input. The subroutines are:  #
#           1) Z Calc - calculates redshift based on a line name and center    #
#           2) E Calc - calculates likely observed line center based on a      #
#              a given name and redshift. Useful when looking for specific     #
#              emission lines.                                                 #
#           3) Ratio Calc - calculates the ratio between line centers. These   #
#              should be unique to line pairs. It will report a bad pair       #
#              if the ratio > 2.5 or the difference between a known ratio and  #
#              the calculated one is greater than 0.01.                        #
#                                                                              #
#  VARIABLES: None                                                             #
#                                                                              #
#  INPUT: Runs from the terminal. Will ask for user input along the way.       #
#         Run with the following command:                                      #
#         python brad_calc.py                                                  #
#  OUTPUT: Prints desired calculated value to the terminal in pretty pretty    #
#          #-outlined boxes. Each subroutine returns its own output.           #
#                                                                              #
# Note: REQUIRES PYTHON 3.5 or later.                                          #
################################################################################
import numpy as np
from colorama import Fore
from colorama import Back
from colorama import Style

################################################################################
#                                                                              #
#  FUNCTION: EM Line Matcher                                                   #
#  PURPOSE: Stores and matches (for recall) the line name with its rest        #
#           wavelength. Wavelengths taken from SDSS line table.                #
#                                                                              #
#  ACCEPTS: Line name in all caps. Input call will format correctly.           #
#  RETURNS: Formatted line name (for terminal) and rest wavelength.            #
#  INPUTS:  None                                                               #
#  OUTPUTS: None                                                               #
#  OUTPUTS To: Nowhere                                                         #
#                                                                              #
################################################################################
def em_line(l_name):
    #Initialize values
    w_out = 0  #The output wavelength
    l_name_out = l_name #The output name.
    #This should probably be done with a dict. But I don't know those yet.
    if l_name == 'LYA': #Lyman-Alpha line. Our favorite.
        w_out = 1215.24
        l_name_out = 'Ly a'
    elif ((l_name == 'OIV')|(l_name == 'SIV')): #Complex is two blended lines.
        w_out = 1399.8
        l_name_out = 'S IV + O IV'
    elif l_name == 'CIV': #Carbon IV. If MgII isn't visible, this is the goto line.
        w_out = 1549.48
        l_name_out = 'C IV'
    elif l_name == 'HEII': #This shows up occassionally I suppose.
        w_out = 1640.4
        l_name_out = 'He II'
    elif l_name == 'CIII': #Shows up more often. Useful for ratio_calc.
        w_out = 1908.734
        l_name_out = 'C III'
    elif l_name == 'CII': #Usually pretty weak.
        w_out = 2326.0
        l_name_out = 'C II'
    elif l_name == 'MGII': #Our most favoritest. Shows up strong at low Z.
        w_out = 2800
        l_name_out = 'Mg II'
    elif l_name == 'LYB': #Occassionally shows up blue of LyA.
        w_out = 1033
        l_name_out = 'Ly b'
    elif l_name == 'OIII': #The red peak of the double horns.
        w_out = 5007
        l_name_out = '[O III]'
    elif l_name == 'HB': #Shows up with OIII. Sometimes a redshift can be found.
        w_out = 4862
        l_name_out = 'H b'
    elif l_name == 'OII': #More for e_calc. Makes it findable for flagging during VI.
        w_out = 3728
        l_name_out = '[O II]'
    elif l_name == 'HG': #H-gamma that can sometimes appear with strength at low z.
        w_out = 4341.68
        l_name_out = 'H g'
    elif l_name == 'HD': #Included for completeness.
        w_out = 4102.89
        l_name_out = 'H d'
    elif l_name == 'HA':
        w_out = 6563
        l_name_out = 'H a'


    return l_name_out,w_out #Send the formatted name and rest wavelength out into the world.

################################################################################
#                                                                              #
#  FUNCTION: Redshift Calculator                                               #
#  PURPOSE: Calculates the redshift for a spectrum based on a line name and a  #
#           line's observed wavelength. User can specify quick select or type  #
#           in a name. It will have to match one of the lines in em_line().    #
#                                                                              #
#  ACCEPTS: Quick Select(chooses a common line), or type the name in all       #
#           lowercase with no spaces. It will capitalize it for the user.      #
#           Also requires the observed wavelength in both cases.               #
#           Example name input: 'lya' for Lyman Alpha, or 'hb' for H-beta.     #
#           Note: roman numerals are used. example: 'mgii'                     #
#  RETURNS: The name of the line used (for checking) and the redshift to 5     #
#           decimal places. Does so in a pretty pretty #-lined box.            #
#  INPUTS:  None                                                               #
#  OUTPUTS: None                                                               #
#  OUTPUTS To: Nowhere                                                         #
#                                                                              #
################################################################################
def z_calc():
    #Formatting and pretty stuff to make visual cues obvious.
    spcrline = '-'*10
    print('\n{0}Z CALC{0}'.format(spcrline))
    #Quick select menu. LyA, CIV, and MgII are used most often for redshift.
    #If the user chooses something other than 1,2,3 it will prompt to type name.
    uch = input('[1] Quick Select\n[2] Type Name?')
    if uch == '1':
        print('\n---Quick Select---')
        print('[1] Ly A\n[2] C IV\n[3] Mg II\n[4] C III')
        qchl = input('Selection: ')
        #Selection uses strings so that an input letter doesn't crash but instead
        #prompts the user to type the name.
        if qchl == '1':
            l_name_in = 'lya'
        elif qchl == '2':
            l_name_in = 'civ'
        elif qchl == '3':
            l_name_in = 'mgii'
        elif qchl == '4':
            l_name_in = 'ciii'
        else:
            print('Not a valid option, try typing the name.') #You didn't choose a given option.
            l_name_in = input('\nTransition Name: ')
    else:
        #TYPE THE NAME IN ALL LOWERCASE WITHOUT SPACES. EXAMPLE: 'lya','mgii'
        l_name_in = input('\nTransition Name: ') #If they choose to type the name. Rare, but used.

    #Either name input method requires they input an observed wavelength.
    l_o = float(input('Observed Lambda: '))
    #Simple redshift calculation.
    l_in = l_name_in.upper()
    l_n_in, l_r = em_line(l_in)
    z = (l_o - l_r) / float(l_r)

    #Redshift printed for the user in a box.
    ident_str = '#     Line: {0:11s}     #'.format(l_n_in)
    out_str = '#     Redshift: {0:05.3f}       #'.format(z)
    breaker = '#'*29
    print('\n')
    print(breaker)
    print(ident_str)
    print(out_str)
    print(breaker)

    #End of program will return the user to the main menu.

################################################################################
#                                                                              #
#  FUNCTION: Emission Line Observable Calculator                               #
#  PURPOSE: Calculates the observed wavelength for an emission line given a    #
#           possible redshift. This will call the lines from em_line().        #
#           The name is input (or quick select) line z_calc(). Useful when     #
#           hunting for an emission line in a noisy spectra. Quick select here #
#           includes O II as that is a common line to hunt for.                #
#                                                                              #
#  ACCEPTS: Quick Select(chooses a common line), or type the em name in all    #
#           lowercase with no spaces. It will capitalize it for the user.      #
#           Also requires the possible redshift in both cases.                 #
#           Example name input: 'lya' for Lyman Alpha, or 'hb' for H-beta.     #
#           Note: roman numerals are used. example: 'mgii'                     #
#  RETURNS: The name of the line used (for checking), the input redshift (for  #
#           checking), and the observable wavelength at that z-value (L_Obs).  #
#           Does so in a pretty pretty #-lined box.                            #
#  INPUTS:  None                                                               #
#  OUTPUTS: None                                                               #
#  OUTPUTS To: Nowhere                                                         #
#                                                                              #
################################################################################
def e_calc():
    #Formatting stuff and boxing strings.
    spcrline = '-'*10
    print('\n{0}E CALC{0}'.format(spcrline))
    #Ask the user for a quick-select line name, or to type it.
    uche = print('[1] Quick Select\n[2] Type Name? ')
    uche = input('Selection: ')
    #For the quick-select, lists the 4 most common lines you hunt for.
    if uche == '1':
        print('\n---Quick Select---')
        print('[1] Ly A\n[2] C IV\n[3] Mg II\n[4] O II\n[5] H b\n[6] C III')
        qchle = input('Selection: ')
        #Selection uses strings so that an input letter doesn't crash but instead
        #prompts the user to type the name.
        if qchle == '1':
            l_name_in = 'lya'
        elif qchle == '2':
            l_name_in = 'civ'
        elif qchle == '3':
            l_name_in = 'mgii'
        elif qchle == '4':
            l_name_in = 'oii'
        elif qchle == '5':
            l_name_in = 'hb'
        elif qchle == '6':
            l_name_in = 'ciii'
        else:
            #If the user did not input a valid number, ask for input.
            print('Not a valid option, try typing the name.')
            l_name_in = input('\nTransition Name: ')
    else:
        #TYPE THE NAME IN ALL LOWERCASE WITHOUT SPACES. EXAMPLE: 'lya','mgii'
        l_name_in = input('\nTransition Name: ')

    #Input the possible redshift here. Common mistake is to accidentally type the
    #'observed wavelength' (like z_calc()). If this is done it will output "Out
    #of Range" most likely.
    l_in = l_name_in.upper()
    l_z = float(input('\nPossible Redshift: '))
    l_n_in, l_r = em_line(l_in)
    l_o = int(float(l_r) * float(l_z + 1))
    if ((l_o < 3600)|(l_o > 9800)):
        l_o = -1

    #Name, redshift used, and observable wavelenth printed in a box for the user.
    #If the observable wavelength is outside the 3600-9800 Angstrom range, it
    #won't be visible in our spectra, so report this to the user.
    ident_str = '#     Line: {0:11s}      #'.format(l_n_in)
    out_str = '#     Redshift: {0:05.3f}        #'.format(l_z)
    if l_o != -1:
        out_str2 = '#     L_Obs: {0:5d} ang       #'.format(l_o)
    else:
        out_str2 = '#     L_Obs: Out of Range    #'
    breaker = '#'*30
    print('\n')
    print(breaker)
    print(ident_str)
    print(out_str)
    print(out_str2)
    print(breaker)

    #End of program will return the user to the main menu.

################################################################################
#                                                                              #
#  FUNCTION: Emission Line Ratio Calculator                                    #
#  PURPOSE: Calculates the ratio between emission line centers. These ratios   #
#           are unique and redshift independent (assuming z is all systemic).  #
#           The redder line is "Line 1 Center" and the bluer line is "Line 2   #
#           Center". This ensures the ratio is always greater than 1. Only     #
#           common ratios are used with the redder line being Mg II, C III, or #
#           C IV.                                                              #
#                                                                              #
#  ACCEPTS: No quick select here. Input the line centers in integer Angstroms  #
#           or at least both in the same units) for line 1 and 2. Use the red, #
#           or right-most line as Line 1 Center. Use the blue, or left line as #
#           Line 2 Center.                                                     #
#           Note: Line centers are input as integers, but are converted to     #
#                 floats for safety.                                           #
#  RETURNS: The likely line pair, if the ratio is close enough to a real ratio.#
#           Some lines will see non-systemic redshifts, so line ratios won't   #
#           exactly match rest wavelength ratios. This subroutine will NOT     #
#           return a suggested pair unless the ratio is less than 2.5 OR       #
#           the closest ratio different to the calculated ratio is less than   #
#           0.01 (absolute difference).                                        #
#           If the ratio IS good, it will report which line is red and which   #
#           blue. It will also show the calculated ratio and the absolute      #
#           difference between it and the closest known ratio.                 #
#           Does so in a pretty pretty #-lined box.                            #
#  INPUTS:  None                                                               #
#  OUTPUTS: None                                                               #
#  OUTPUTS To: Nowhere                                                         #
#                                                                              #
################################################################################
def ratio_calc():
    #Formatting stuff for the pretty boxes and lines for visual cues.
    mg2check = 0
    spcrline = '-'*10
    numbox = '#'*53
    pline = ' '*16
    pline2 = ' '*21
    boxstr = '#'*50
    print('\n{0}RATIO CALC{0}'.format(spcrline))

    #Ask for the ratios. No error trapping here. You put in letters, it crashes.
    #Remember: Red = Line 1, Blue = Line 2.
    l1c = float(input('Line 1 Center: '))
    l2c = float(input('Line 2 Center: '))

    #Calculate the ratio.
    if l1c >= l2c:
        ratc = l1c / l2c
    else:
        ratc = l2c / l1c

    #Define the known ratios with the names of the red/blue lines. Uses a numpy
    #structured array so I can call columns by a name and can store strings and
    #floats in the same array.
    ratarr = np.zeros(20,dtype=[('LINE1_NAME','U7'),('LINE2_NAME','U7'),('RATIO','f4')])
    ratarr['LINE1_NAME'] = np.array(['Mg II','Mg II','Mg II','Mg II','C III','C III','C III','C IV','C IV',
                            'O III','O III','O III','H b','H b','H b','H g','H g','O II','H a','Ly B'])
    ratarr['LINE2_NAME'] = np.array(['C III','C IV','SIV+OIV','Ly A','C IV','SIV+OIV','Ly A', 'SIV+OIV','Ly A',
                            'Mg II','C III','H g','Mg II','C III','H g','Mg II','C III','Mg II','H b','Ly A'])
    ratarr['RATIO'] = np.array([1.4675,1.806,2.00,2.3026,1.23097,1.36286,1.56908,1.10714,1.27467,
                            1.78857,2.6237,1.1535,1.7364,2.547,1.1198,1.5506,2.2746,1.33143,1.34985,1.17642])

    #Find the smallest absolute difference between the calculated and known ratios.
    #The smallest difference should be the right ratio (if it's close enough).
    rdiff = np.absolute(ratc - ratarr['RATIO'])
    mindiff = np.amin(rdiff)
    ml_val = np.where(rdiff == np.amin(rdiff))[0]
    if ((ratarr['LINE1_NAME'][ml_val[0]] == 'Mg II')|(ratarr['LINE2_NAME'][ml_val[0]] == 'Mg II')):
        mg2check = 1
        name_arr = np.array([ratarr['LINE1_NAME'][ml_val[0]],ratarr['LINE2_NAME'][ml_val[0]]])
        wname = np.where(name_arr == 'Mg II')[0]
        if wname[0] == 0:
            l_obs = l1c
        elif wname[0] == 1:
            l_obs = l2c
        l_rest = 2800
        z_out = (l_obs - l_rest) / float(l_rest)
        l_name = 'Mg II'
    #If the difference is too large, or the ratio is larger than one we would select
    #warn the user and return to the main menu.
    if ((ratc >= 2.5)|(mindiff > 0.01)):
        print('Calculated ratio: R = {0:01.4f}'.format(ratc))
        print('\n')
        print(boxstr)
        print('#  Not a recognized pair, try a different pair.  #')
        print(boxstr)
    #If the ratio is good, return the most likely pair Red first (left most line
    #center used), Blue second. Also return the calculated ratio from user input
    #line centers and the difference between tha ratio and the closest one found.
    #Does it in a box.
    elif mg2check == 1:
        print('\n')
        print(numbox)
        print('#  Most likely line pair:'+' '*27+'#')
        print('#  Red Line: {0} | Blue Line: {1}{2}#'.format(ratarr['LINE1_NAME'][ml_val[0]],ratarr['LINE2_NAME'][ml_val[0]],pline))
        print('#  Calculated ratio: R = {0:01.4f}{1}#'.format(ratc,pline2))
        print('#  Difference between calculated and table: {0:01.4f}  #'.format(np.amin(rdiff)))
        print(numbox)
        print('\n')
        extra_str = 'USING THE INPUT LAMBDA FOR MGII'
        ident_str = '#     Line: {0:11s}     #'.format(l_name)
        wave_str =  '#     Wavelength: {0:04d}      #'.format(int(l_obs))
        out_str = '#     Redshift: {0:05.3f}       #'.format(z_out)
        breaker = '#'*29
        print(extra_str)
        print(breaker)
        print(ident_str)
        print(wave_str)
        print(out_str)
        print(breaker)
    else:
        print('\n')
        print(numbox)
        print('#  Most likely line pair:'+' '*27+'#')
        print('#  Red Line: {0} | Blue Line: {1}{2}#'.format(ratarr['LINE1_NAME'][ml_val[0]],ratarr['LINE2_NAME'][ml_val[0]],pline))
        print('#  Calculated ratio: R = {0:01.4f}{1}#'.format(ratc,pline2))
        print('#  Difference between calculated and table: {0:01.4f}  #'.format(np.amin(rdiff)))
        print(numbox)

    #End of program will return the user to the main menu.

################################################################################
#                                                                              #
#  FUNCTION: Broad Absorption Lines Check for CIV                              #
#  PURPOSE: Calculates the width of the trough blueward of CIV to see if it    #
#           meets the 2000 km/s width criteria.                                #
#                                                                              #
#  ACCEPTS: Wavelength inputs in Angstroms for the left and right side of the  #
#           trough and the systemic redshift (z) of the QSO.                   #
#  RETURNS: The rest-frame width in km/s of the trough and a BAL flag          #
#           recommendation.                                                    #
#  INPUTS:  None                                                               #
#  OUTPUTS: None                                                               #
#  OUTPUTS To: Nowhere                                                         #
#                                                                              #
################################################################################
def bal_check():
    #Formatting stuff for the pretty boxes and lines for visual cues.
    bal_flag = 'NOT BAL'
    spcrline = '-'*10
    numbox = '#'*53
    pline = ' '*16
    pline2 = ' '*21
    boxstr = '#'*50
    print('\n{0}BAL CHECK{0}'.format(spcrline))

    #Ask for the trough edges. No error trapping here. You put in letters, it crashes.
    #Left/right doesn't matter, it's taking an absolute difference.
    l1e = float(input('Trough Edge 1 (Ang): '))
    l2e = float(input('Trough Edge 2 (Ang): '))
    zsys = float(input('QSO Redshift: '))

    #Calculate z values for the two edges and the difference.
    l1ez = (l1e - 1550) / 1550
    l2ez = (l2e - 1550) / 1550
    z_diff = np.abs(l1ez - l2ez)

    #Calculate the velocity difference.
    v_diff = z_diff * 2.99792458e5 #multiplied by speed of light in km/s

    #Now move it to the rest frame:
    v_rest = int(v_diff / (1 + zsys))
    if v_rest >= 2000:
        bal_flag = 'BAL (2)'
    elif ((v_rest>=1000)&(v_rest<2000)):
        bal_flag = 'Possible BAL (20)'
    else:
        bal_flag = 'NOT BAL'

    #If the velocity in rest frame is above 2000 km/s, call it a BAL. If it's
    #between 1000 and 2000, give it the BAL? flag (20). Otherwise, not a BAL.
    print('\n')
    print(numbox)
    print('#  Velocity Difference (rest): {0:5d} km/s           #'.format(v_rest))
    print('#' + ' '*51 + '#')
    print('#      BAL? {:17}                       #'.format(bal_flag))
    print(numbox)

    #End of program will return the user to the main menu.

################################################################################
#                                                                              #
#  FUNCTION: Emission Line Observable Calculator 2                             #
#  PURPOSE: Calculates the observed wavelengths for a number of emission lines #
#           from the given redshift. This will call the lines from em_line().  #
#           The name is input (or quick select) line z_calc(). Useful when     #
#           hunting for an emission line in a noisy spectra.                   #
#                                                                              #
#  ACCEPTS: Quick Select(chooses a common line), or type the em name in all    #
#           lowercase with no spaces. It will capitalize it for the user.      #
#           Also requires the possible redshift in both cases.                 #
#           Example name input: 'lya' for Lyman Alpha, or 'hb' for H-beta.     #
#           Note: roman numerals are used. example: 'mgii'                     #
#  RETURNS: The name of the line used (for checking), the input redshift (for  #
#           checking), and the observable wavelength at that z-value (L_Obs).  #
#           Does so in a pretty pretty #-lined box.                            #
#  INPUTS:  None                                                               #
#  OUTPUTS: None                                                               #
#  OUTPUTS To: Nowhere                                                         #
#                                                                              #
################################################################################
def e_calc2():
    #Formatting stuff and boxing strings.
    spcrline = '-'*10
    print('\n{0}E CALC2{0}'.format(spcrline))
    #Ask the user for a possible redshift:
    uche = float(input('Possible Redshift: '))
    l_name_arr = np.array(['lya','siv','civ','ciii','mgii','oii','hb','oiii'])
    breaker = '#'*30
    out_str = '#     Redshift: {0:05.3f}        #'.format(uche)
    print('\n')
    print(breaker)
    print(out_str)
    for lna in l_name_arr:
        l_name_in = lna
        l_in = l_name_in.upper()
        l_n_in, l_r = em_line(l_in)
        l_o = int(float(l_r) * float(uche + 1))
        if l_name_in == 'siv':
            l_n_in = 'S/O IV'
        if l_o > 10000:
            padstr = ' '*12
        else:
            padstr = ' '*13
        if ((l_o >= 9800)|(l_o<=3600)):
            line_str = '#  {}{:7}: {:4d}{}{}#'.format(Fore.RED,l_n_in,l_o,Style.RESET_ALL,padstr)
        elif ((l_o>3600)&(l_o<=7200)):
            line_str = '#  {}{}{:7}: {:4d}{}{}#'.format(Fore.GREEN,Style.BRIGHT,l_n_in,l_o,Style.RESET_ALL,padstr)
        else:
            line_str = '#  {}{:7}: {:4d}{}{}#'.format(Fore.YELLOW,l_n_in,l_o,Style.RESET_ALL,padstr)
        print(line_str)
    print(breaker)

    #End of program will return the user to the main menu.

def v_diff():
    #Formatting and pretty stuff to make visual cues obvious.
    spcrline = '-'*10
    print('\n{0}DELTA V{0}\n'.format(spcrline))

    ckms = 299792.458
    z1 = float(input('Z_1: '))
    zt = float(input('Z_true: '))
    dv = (ckms*np.abs(z1-zt))/(1+zt)

    #delta_v printed for the user in a box.
    ident_str = '#     Z_1:    {:5.3f}         #'.format(z1)
    ident_str2 = '#     Z_true: {:5.3f}         #'.format(zt)
    out_str = '#     del_v: {:6,d} km/s   #'.format(int(dv))
    breaker = '#'*29
    print('\n')
    print(breaker)
    print(ident_str)
    print(ident_str2)
    print(out_str)
    print(breaker)

    #End of program will return the user to the main menu.


################################################################################
#                                                                              #
#  FUNCTION: Main Menu (active)                                                #
#  PURPOSE: Creates a main menu for a user to select a subroutine. Active menu #
#           means that when a subroutine completes, it returns to this menu    #
#           and only exits on user selection (Quit) or if it crashes. I got    #
#           tired of rerunning the program from the terminal for every         #
#           spectrum I looked at. This will keep it open until I'm done with   #
#           spectra for that session.                                          #
#                                                                              #
#  ACCEPTS: A number defined by the menu. Uses a string instead of an int so   #
#           that a letter can be used to quit instead of crash.                #
#  RETURNS: The subroutine called by that selection. See subroutine definitions#
#           above for what each one does.                                      #
#  INPUTS:  None                                                               #
#  OUTPUTS: None                                                               #
#  OUTPUTS To: Nowhere                                                         #
#                                                                              #
################################################################################
def main_menu():
    #This while loop is what keeps the menu open between subroutine calls.
    while True:
        #Formatting stuff for the box and lines seen. Startstr is used to visually
        #separate the end of one subroutine call from a menu restart so it's easier
        #to see the subroutine outputs.
        startstr = '!'*50
        boxstr = '#'*50
        linestr = '-'*16
        spcrstr = ' '*9
        byespcr = ' '*22
        print('\n')
        print(startstr)
        print(startstr)
        print('\n')
        print(boxstr)
        print("#{0} Brad's EZ Redshift Calculator{0}#".format(spcrstr))
        print(boxstr)
        print('\n')
        print(linestr)
        print('Which Subprogram')
        print(linestr)
        #Ask the user for their selection. Numbers are used, but letters are
        #added in case users can't follow instructions.
        print('[1] Z Calc\n[2] E Calc\n[3] Ratio Calc\n[4] BAL Check\n[5] E Calc 2\n[6] Delta V\n[7] Quit')
        print(linestr)
        n = input('Selection: ')
        if ((n == '1')|(n == 'Z')|(n == 'z')):
            z_calc()
        elif ((n == '2')|(n == 'E')|(n == 'e')):
            e_calc()
        elif ((n == '3')|(n == 'R')|(n == 'r')):
            ratio_calc()
        elif ((n == '4')|(n == 'B')|(n == 'b')):
            bal_check()
        elif ((n == '5')|(n == 'E2')|(n == 'e2')):
            e_calc2()
        elif ((n== '6')|(n == 'D')|(n == 'd')):
            v_diff()
        else:
            print('\n')
            print(boxstr)
            print('#{0}Bye!{0}#'.format(byespcr))
            print(boxstr)
            print('\n')
            break

#Calls the main menu so I can launch this from the terminal command.
main_menu()
