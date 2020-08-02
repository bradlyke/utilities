class SciCons:
    def __init__(self,usys):
        self.sec = self.get_cons('ext_defs')['sec'] #Number of Cs-133 transitions per second
        self.yrd = self.get_cons('ext_defs')['yrd'] #Number of days in Julian year
        self.yrs = self.get_cons('ext_defs')['yrs'] #Number of sec in a Julian year
        self.c = self.get_cons(usys)['c'] #Speed of light (exact)
        self.G = self.get_cons(usys)['G'] #Gravitational constant
        self.Na = self.get_cons('ext_defs')['Na'] #Avogadro's constant (exact)
        self.R = self.get_cons(usys)['R'] #Ideal gas constant (exact)
        self.e0 = self.get_cons(usys)['e0'] #Vacuum permissivity
        self.pi = self.get_cons('ext_defs')['pi'] #For getting the defined value of pi used
        self.mu0 = self.get_cons(usys)['mu0'] #Vacuum permeability
        self.h = self.get_cons(usys)['h'] #J s (exact)
        self.hbar = self.h / (2 * self.pi) #reduced plank constant, J s / rad
        self.kb = self.get_cons(usys)['kb'] #Boltzmann constant (exact)
        self.sb = self.get_cons(usys)['sb'] #Stefan-Boltzmann constant
        self.e = self.get_cons('ext_defs')['e'] #Elementary charge (exact)
        self.me = self.get_cons(usys)['me'] #Electron mass
        self.mp = self.get_cons(usys)['mp'] #Proton mass
        self.mn = self.get_cons(usys)['mn'] #Neutron mass
        self.amu = self.get_cons(usys)['amu'] #Atomic mass constant (same as Da, or u)

        #Atomic hydrogen mass in proper units
        if usys == 'nat':
            #Nat units: amu
            self.mH = self.get_cons('ext_defs')['mH']
        else:
            #MKS units: kg | CGS units: g
            self.mH = self.get_cons('ext_defs')['mH'] * self.amu

        self.a0 = self.get_cons(usys)['a0'] #Bohr radius
        self.muB = self.get_cons(usys)['muB'] #Bohr magneton
        self.Msun = self.get_cons(usys)['Msun'] #Solar mass
        self.Rsun = self.get_cons(usys)['Rsun'] #Solar radius
        self.Lsun = self.get_cons(usys)['Lsun'] # Solar luminosity
        self.Tsun = self.get_cons('ext_defs')['Tsun'] #Effective solar temperature
        self.Me = self.get_cons(usys)['Me'] #Earth mass
        self.Re = self.get_cons(usys)['Re'] #Earth mean Radius
        self.Mj = self.get_cons(usys)['Mj'] #Jupiter mass
        self.Rj = self.get_cons(usys)['Rj'] #Jupiter radius
        self.Rc = self.get_cons(usys)['Rc'] #Rydberg constant at infinity

        #Radiation density constant in appropriate units
        if usys == 'nat':
            #Nat units: eV/cm^3/K^4
            self.a = (4 * self.sb) / self.get_cons('cgs')['c']
        else:
            #MKS units: J / m^3 / K^4 | CGS units: erg/m^3/K^4
            self.a = (4 * self.sb) / self.c

        self.eV = self.get_cons(usys)['eV'] #eV--J relation, see get_cons for units
        self.au = self.get_cons(usys)['au'] #Astronomical unit (exact)
        self.ly = self.yrs* self.c #Lightyear
        self.pc = self.get_cons(usys)['pc'] #Parsec
        self.wb = self.get_cons(usys)['wb'] #Wein's constant in m K
        self.alpha = self.get_cons('ext_defs')['alpha'] #Fine-structure constant
        if usys == 'nat':
            self.MsunE = self.Msun / self.Me #Solar mass in Earth masses
            self.RsunE = self.Rsun / self.Re #Solar radius in Earth radii
            self.MsunJ = self.Msun / self.Mj #Solar mass in Jupiter masses
            self.RsunJ = self.Rsun / self.Rj #Solar radius in Jupiter radii
            self.MeS = self.Me / self.Msun #Earth mass in Solar masses
            self.ReS = self.Re / self.Rsun #Earth radius in Solar radii
            self.MeJ = self.Me / self.Mj #Earth mass in Jupiter masses
            self.ReJ = self.Re / self.Rj #Earth radius in Jupiter radii
            self.MjS = self.Mj / self.Msun #Jupiter mass in Solar masses
            self.RjS = self.Rj / self.Rsun #Jupiter radius in Solar radii
            self.MjE = self.Mj / self.Me #Jupiter mass in Earth masses
            self.RjE = self.Rj / self.Re #Jupiter radius in Earth radii

            self.aupc = self.au / self.pc #1 au in parsec
            self.auly = self.au / self.ly #1 au in ly ---EXACT
            self.lyau = self.ly / self.au #1 ly in au ---EXACT
            self.lypc = self.ly / self.pc #1 ly in pc
            self.pcau = self.pc / self.au #1 parsec in au
            self.pcly = self.pc / self.ly #1 parsec in ly

        self.unit_defs = self.get_defs(usys)

    def get_cons(self,usys):
        ext_defs = {
            "sec" : 9.192631770e9, #1 second in number of Cs-133 hyperfine transitions
            "yrd" : 365.25, # days
            "yrs" : 31557600.0, # seconds
            "Na" : 6.02214076e23, # unitless
            "pi" : 3.141592653589793,
            "e" : 1.602176634e-19, # coulombs
            "mH" : 1.00782503223, # amu
            "Tsun" : 5777, # K
            "alpha" : 7.2973525693e-3 # unitless
        }

        mksv = {
            "c" : 2.99792458e8, # m/s
            "G" : 6.67430e-11, # N m^2/kg^2
            "R" : 8.31446261815324, # J/mol/K
            "e0" : 8.8541878128e-12, # F/m
            "mu0" : 1.25663706212e-6, # H/m
            "h" : 6.62607015e-34, # J s
            "kb" : 1.380649e-23, # J/K
            "sb" : 5.670374419e-8, # W/m^2/K^4
            "me" : 9.1093837015e-31, # kg
            "mp" : 1.67262192369e-27, # kg
            "mn" : 1.67492749804e-27, # kg
            "amu" : 1.66053906660e-27, # kg
            "a0" : 5.29177210903e-11, # m
            "muB" : 9.2740100783e-24, # J/tesla
            "Msun" : 1.98855e30, # kg
            "Rsun" : 6.96392e8, # m
            "Lsun" : 3.828e26, # W
            "Me" : 5.97237e24, # kg
            "Re" : 6.371e6, # m
            "Mj" : 1.8982e27, # kg
            "Rj" : 6.9911e7, # m
            "Rc" : 1.0973731568160e7, # 1/m
            "eV" : 1.602176634e-19, # 1 eV in joules
            "au" : 1.49597870700e11, # m
            "pc" : 3.0856775814671916e16, # m
            "wb" : 2.897771955e-3 # m K
        }

        cgsv = {
            "c" : mksv['c'] * 1e2, # cm/s
            "G" : 6.67430e-13, # N cm^2/g^2
            "R" : mksv['R'] * 1e7, # erg/mol/K
            "e0" : 8.8541878128e-14, # F/cm
            "mu0" : 1.25663706212e-8, # H/cm
            "h" : 6.62607015e-27, # erg s
            "kb" : 1.380649e-16, # erg/K
            "sb" : 5.670374419e-5, # erg/s/cm^2/K^4
            "me" : 9.1093837015e-28, # g
            "mp" : 1.67262192369e-24, # g
            "mn" : 1.67492749804e-24, # g
            "amu" : 1.66053906660e-24, # g
            "a0" : 5.29177210903e-9, # cm
            "muB" : 9.2740100783e-21, # erg/gauss
            "Msun" : 1.98855e33, # g
            "Rsun" : 6.96392e10, # cm
            "Lsun" : 3.828e33, # erg/s
            "Me" : 5.97237e27, #  g
            "Re" : 6.371e8, # cm
            "Mj" : 1.8982e30, # g
            "Rj" : 6.9911e9, # cm
            "Rc" : 1.0973731568160e5, # 1/cm
            "eV" : 1.602176634e-12, # 1 eV in erg
            "au" : 1.49597870700e13, # cm
            "pc" : 3.0856775814671916e18, # cm
            "wb" : 0.2897771955 # cm K
        }

        natv = {
            "c" : mksv['c'] * 1e10, # Ang/s
            "G" : 4.30091e-3, # pc (km/s)^2 per Solar mass
            "R" : mksv['R'] / ext_defs['e'], # eV/mol/K
            "e0" : cgsv['e0'], # F/cm
            "mu0" : cgsv['mu0'], # H/cm
            "h" : mksv['h'] / ext_defs['e'], # eV s
            "kb" : mksv['kb'] / ext_defs['e'], # eV/K
            "sb" : (mksv['sb'] * 1e-4) / ext_defs['e'], # eV/s/cm^2/K^4
            "me" : 0.51099895000e6, # eV/c^2
            "mp" : 938.27208816e6, # eV/c^2
            "mn" : 939.56542052e6, # eV/c^2
            "amu" : (mksv['amu']*mksv['c']**(2))/ext_defs['e'], # eV/c^2
            "a0" : 0.529177210903, # Ang
            "muB" : 5.7883818060e-5, # eV/T
            "Msun" : 1.98855e30, # kg
            "Rsun" : 6.96392e8, # m
            "Lsun" : 3.828e26, # W
            "Me" : 5.97237e24, # kg
            "Re" : 6.371e6, # m
            "Mj" : 1.8982e27, # kg
            "Rj" : 6.9911e7, # m
            "Rc" : mksv['Rc'] * 1e-10, # 1/Ang
            "eV" : 1 / mksv['eV'], # 1 J in electron volts
            "au" : 1.49597870700e11, # m
            "pc" : 3.0856775814671916e16, # m
            "wb" : mksv['wb'] * 1e10 # Ang K
        }

        unit_systems = {
            "mks" : mksv,
            "cgs" : cgsv,
            "nat" : natv,
            "ext_defs" : ext_defs
        }

        return unit_systems[usys]

    def get_defs(self,usys):
        mksu = {
            "sec" : "Number of Cs-133 transitions per second",
            "yrd" : "Number of days in a Julian year",
            "yrs" : "Number of sec in a Julian year",
            "c" : "Speed of light, m/s (exact)",
            "G" : 'Gravitational constant, N m^2/kg^2',
            "Na" : "Avogadro's constant, unitless (exact)",
            "R": "Ideal gas constant, J/mol/K (exact)",
            "e0" : "Vacuum permissivity, F/m",
            "mu0" : "Vacuum permeability, H/m",
            "h" : "Planck constant, J s (exact)",
            "hbar" : "Reduced plank constant, J s/rad",
            "kb" : "Boltzmann constant, J/K (exact)",
            "sb" : "Stefan-Boltzmann constant, W/m^2/K^4",
            "e" : "Elementary charge, coulombs (exact)",
            "me" : "Electron mass, kg",
            "mp" : "Proton mass, kg",
            "mn" : "Neutron mass, kg",
            "amu" : "Atomic mass constant, kg (same as Da, or u)",
            "mH" : "Atomic hydrogen mass (H I), kg",
            "a0" : "Bohr radius, m",
            "muB" : "Bohr magneton, J/tesla",
            "Msun" : "Solar mass, kg",
            "Rsun" : "Solar radius, m",
            "Lsun" : "Solar luminosity, W",
            "Tsun" : "Effective solar temperature, K",
            "Me" : "Earth mass, kg",
            "Re" : "Earth mean Radius, m",
            "Mj" : "Jupiter mass, kg",
            "Rj" : "Jupiter radius, m",
            "Rc" : "Rydberg constant at infinity, 1/m",
            "a" : "Radiation density constant, J / m^3 / K^4",
            "eV" : "1 electron volt, J (exact)",
            "au" : "Astronomical unit, m (exact)",
            "ly" : "Lightyear, m (exact)",
            "pc" : "Parsec, m",
            "wb" : "Wein's constant, m K",
            "alpha" : "Fine-structure constant, unitless"
        }

        cgsu = {
            "sec" : "Number of Cs-133 transitions per second",
            "yrd" : "Number of days in a Julian year",
            "yrs" : "Number of sec in a Julian year",
            "c" : "Speed of light, cm/s (exact)",
            "G" : 'Gravitational constant, N cm^2/g^2',
            "Na" : "Avogadro's constant, unitless (exact)",
            "R": "Ideal gas constant, erg/mol/K (exact)",
            "e0" : "Vacuum permissivity, F/cm",
            "mu0" : "Vacuum permeability, H/cm",
            "h" : "Planck constant, erg s (exact)",
            "hbar" : "Reduced plank constant, erg s/rad",
            "kb" : "Boltzmann constant, erg/K (exact)",
            "sb" : "Stefan-Boltzmann constant, erg/s/cm^2/K^4",
            "e" : "Elementary charge, coulombs (exact)",
            "me" : "Electron mass, g",
            "mp" : "Proton mass, g",
            "mn" : "Neutron mass, g",
            "amu" : "Atomic mass constant, g (same as Da, or u)",
            "mH" : "Atomic hydrogen mass (H I), g",
            "a0" : "Bohr radius, cm",
            "muB" : "Bohr magneton, erg/gauss",
            "Msun" : "Solar mass, g",
            "Rsun" : "Solar radius, cm",
            "Lsun" : "Solar luminosity, erg/s",
            "Tsun" : "Effective solar temperature, K",
            "Me" : "Earth mass, g",
            "Re" : "Earth mean Radius, cm",
            "Mj" : "Jupiter mass, g",
            "Rj" : "Jupiter radius, cm",
            "Rc" : "Rydberg constant at infinity, 1/cm",
            "a" : "Radiation density constant, erg/m^3/ K^4",
            "eV" : "1 electron volt, erg (exact)",
            "au" : "Astronomical unit, cm (exact)",
            "ly" : "Lightyear, cm (exact)",
            "pc" : "Parsec, cm",
            "wb" : "Wein's constant, cm K",
            "alpha" : "Fine-structure constant, unitless"
        }

        natu = {
            "sec" : "Number of Cs-133 transitions per second",
            "yrd" : "Number of days in a Julian year",
            "yrs" : "Number of sec in a Julian year",
            "c" : "Speed of light, Ang/s (exact)",
            "G" : 'Gravitational constant, pc (km/s)^2 per Solar mass',
            "Na" : "Avogadro's constant, unitless (exact)",
            "R": "Ideal gas constant, J/mol K (exact)",
            "h" : "Planck constant, eV s (exact)",
            "hbar" : "Reduced plank constant, eV s/rad",
            "kb" : "Boltzmann constant, eV/K (exact)",
            "sb" : "Stefan-Boltzmann constant, eV/s/cm^2/K^4",
            "e" : "Elementary charge, coulombs (exact)",
            "me" : "Electron mass, eV/c^2 (not MeV)",
            "mp" : "Proton mass, eV/c^2",
            "mn" : "Neutron mass, eV/c^2",
            "amu" : "Atomic mass constant, eV/c^2",
            "mH" : "Atomic hydrogen mass (H I), eV/c^2",
            "a0" : "Bohr radius, Ang",
            "muB" : "Bohr magneton, eV/tesla",
            "Msun" : "Solar mass, kg",
            "Rsun" : "Solar radius, m",
            "Lsun" : "Solar luminosity, W",
            "Tsun" : "Effective solar temperature, K",
            "Me" : "Earth mass, kg",
            "Re" : "Earth mean Radius, m",
            "Mj" : "Jupiter mass, kg",
            "Rj" : "Jupiter radius, m",
            "Rc" : "Rydberg constant at infinity, 1/Ang",
            "a" : "Radiation density constant, eV/cm^3/K^4",
            "eV" : "1 joule in eV (exact)",
            "au" : "Astronomical unit, m (exact)",
            "ly" : "Lightyear, m (exact)",
            "pc" : "Parsec, m",
            "wb" : "Wein's constant, Ang K",
            "alpha" : "Fine-structure constant, unitless",
            "MsunE" : "Solar mass in Earth masses",
            "RsunE" : "Solar radius in Earth radii",
            "MsunJ" : "Solar mass in Jupiter masses",
            "RsunJ" : "Solar radius in Jupiter radii",
            "MeS" : "Earth mass in Solar masses",
            "ReS" : "Earth radius in Solar radii",
            "MeJ" : "Earth mass in Jupiter masses",
            "ReJ" : "Earth radius in Jupiter radii",
            "MjS" : "Jupiter mass in Solar masses",
            "RjS" : "Jupiter radius in Solar radii",
            "MjE" : "Jupiter mass in Earth masses",
            "RjE" : "Jupiter radius in Earth radii",
            "aupc" : "Astronomical unit, pc",
            "auly" : "Astronomical unit, ly (exact)",
            "lyau" : "Lightyear, au (exact)",
            "lypc" : "Lightyear, pc",
            "pcau" : "Parsec, au",
            "pcly" : "Parsec, ly"
        }

        unit_systems = {
            "mks" : mksu,
            "cgs" : cgsu,
            "nat" : natu,
        }

        unit_defs = unit_systems[usys]

        return unit_defs

    def defs(self):
        print()
        print(61*'-')
        for keyval in self.unit_defs.items():
            kdisp = keyval[0]
            if kdisp == 'yrs':
                print('     {:6s}: {:45s}'.format(keyval[0],keyval[1]))
                print(61*'-')
            else:
                #print('|     {:6s}: {:45s} |'.format(keyval[0],keyval[1]))
                print('     {:6s}: {:45s}'.format(keyval[0],keyval[1]))
        print(61*'-')

# Import SciCons so you can load a set of constants in a given unit systems. Ex.
#   [1]: from sciCon_v2 import SciCons
#   [2]: mks = SciCons('mks')
#   [3]: print(mks.c)
#   299792458.0
if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='A set of scientific constants for MKS, CGS, or "Natural" units.')
    parser.add_argument('u', default='mks', choices = ['mks','cgs','nat'],
                        help="The system of units to use.")
    args = parser.parse_args()

    usys_name = SciCons(args.u)
    usys_name.defs()
