#Import the tools necessary to work with spectra files.
import cat_tools as ct #This is a tool I wrote, so you'll need the script in your path.
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as ticker
matplotlib.rc('text',usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
matplotlib.rc('xtick',direction='in',top=True)
matplotlib.rc('xtick.major',size=8)
matplotlib.rc('xtick.minor',visible=True,size=5)
matplotlib.rc('ytick',direction='in',right=True)
matplotlib.rc('ytick.major',size=8)
matplotlib.rc('ytick.minor',visible=True,size=5)
matplotlib.rc('legend',framealpha=0)

#A spectrum is a class. The object has 3 attributes right now: self, boxcar,
#plot_spec.
class Spectrum:
    #initialize the object by loading the file. Needs a file name for infile.
    def __init__(self,infile):
        self.infile = infile #Load the filename as a thing for later use in plot saving.
        self.data = ct.file_load(infile) #Load the file as a numpy structured array.
        self.loglam = self.data['loglam'] #SDSS stores wavelengths as log10 values.
        self.lam = 10**self.loglam #Convert this to a decimal wavelength in Angstroms.
        self.flux = self.data['flux'] #Load the flux. Absolute, units are in plot.
        self.ivar = self.data['ivar'] #Load the inverse variance. Needed for plotting and boxcar.
        self.ferr = self.data['ivar']**(-0.5)

    #We want to be able to smooth it fairly easily. Smooth_pct = 10 is good for some stuff.
    def boxcar(self,flux_arr,flux_var,smooth_pct,weight=False):
        box_size = smooth_pct
        if box_size%2 != 0:
            box_size += 1
        num_dpoints = len(flux_arr)
        box_flux = np.zeros(num_dpoints,dtype='f8')
        box_err = np.zeros(num_dpoints,dtype='f8')
        box_flux[:] = flux_arr[:]
        box_err[:] = flux_var[:]

        for i in range(num_dpoints):
            #Define range of values to smooth
            if i < int(box_size/2):
                lower = 0
            else:
                lower = -int(box_size/2) + i
            if ((int(box_size/2) + i) > num_dpoints):
                upper = num_dpoints
            else:
                upper = int(box_size/2) + i
            #Smooth the values depending on smoothing type
            if weight==True:
                noise_temp = np.sqrt(box_err[lower:upper])
                signal_temp = box_flux[lower:upper]
                flux_temp = np.average(signal_temp,weights=noise_temp)
                ivar_temp = (np.average((signal_temp-flux_temp)**2, weights=noise_temp))**(-1.0)
            else:
                signal_temp = box_flux[lower:upper]
                flux_temp = np.median(signal_temp)
                ivar_temp = (np.average((signal_temp-flux_temp)**2))**(-1.0)

            box_flux[i] = flux_temp
            box_err[i] = ivar_temp

        return box_flux,box_err


    def plot_spec(self, plot_size, rest_twin=False, spec_redshift=0.0,
                    smooth=False, smooth_box=10, weighted_avg=False,
                    err=False, sky=False, scale_sky=False,
                    save=False, file_type='png'):
        wobs = np.where((self.lam>=3700)&(self.lam<=10000))[0]
        flux_range = np.amax(self.flux[wobs]) - np.amin(self.flux[wobs])
        flux_pad = float(flux_range)/10
        y_lower = np.amin(self.flux[wobs]) - flux_pad
        y_upper = np.amax(self.flux[wobs]) + flux_pad
        x_lower = np.amin(self.lam)
        x_upper = np.amax(self.lam)

        def round_updown(tval):
            rup = tval + (-tval %100)
            rdo = tval + (-tval % -100)
            rmid = int((rup+rdo)/2)
            if tval >= rmid:
                y = rup
            else:
                y = rdo
            return y

        def rest_tick_step_gen(obs_min,obs_max,z_in):
            if z_in < 0.5:
                return 2000
            else:
                rest_max = obs_max / (1 + z_in)
                rest_min = obs_min / (1 + z_in)
                step_size = round_updown(int((rest_max - rest_min)/5))
                return step_size

        rest_ticks = rest_tick_step_gen(x_lower,x_upper,spec_redshift)

        if plot_size=='s':
            matplotlib.rc('font',size=11)
            fig,ax = plt.subplots(figsize=(5,4))
        elif plot_size=='l':
            matplotlib.rc('font',size=11)
            fig,ax = plt.subplots(figsize=(10,4))
        else:
            matplotlib.rc('font',size=24)
            fig,ax = plt.subplots(figsize=(15,12))

        if rest_twin==True:
            axT = ax.twiny()
        if smooth==True:
            self.bflux,self.berr = self.boxcar(self.flux,self.ivar,smooth_box,weight=weighted_avg)
            ax.plot(self.lam,self.flux,color='0.70',linewidth=0.8)
            ax.plot(self.lam,self.bflux,color='black',linewidth=0.6)
        else:
            ax.plot(self.lam,self.flux,color='black')
        if sky==True:
            if scale_sky==True:
                self.sky = self.data['sky'] / 10.0
            else:
                self.sky = self.data['sky']
            ax.plot(self.lam,self.sky,color='green',linewidth=0.7,alpha=0.5)
        if err==True:
            ax.plot(self.lam,self.ferr,color='red',linewidth=0.6)
        ax.set_xlim((x_lower,x_upper))
        ax.set_xticks(np.arange(4000,11000,1000))
        ax.set_xlabel(r'\textbf{Observed Frame Wavelength (\AA)}')
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(100))
        if rest_twin==True:
            axT.set_xlim((x_lower/(1+spec_redshift),x_upper/(1+spec_redshift)))
            axT.set_xlabel(r'\textbf{Rest Frame Wavelength (\AA)}')
            axT.xaxis.set_major_locator(ticker.MultipleLocator(rest_ticks))
            axT.xaxis.set_minor_locator(ticker.AutoMinorLocator(4))
        ax.set_ylim((y_lower,y_upper))
        ax.set_ylabel(r'$f_{\lambda}$ ($10^{-17}$ \textbf{ergs s}$^{-1}$ \textbf{cm}$^{-2}$\,\textbf{\AA}$^{-1}$)')

        name_str = r'\textbf{SDSS J003713.64+241121.5, z = %.3f}'%spec_redshift
        ax.text(0.6,0.95,name_str,transform=ax.transAxes, verticalalignment='top',
                    bbox=dict(boxstyle='square,pad=0.2',fc='magenta',alpha=0.0))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))

        plt.tight_layout()
        if save==True:
            if file_type=='eps':
                fname_out = self.infile.replace('.fits','_lrg.eps')
                fname_out = fname_out.replace('/data/','/plots/')
                fig.savefig(fname_out,bbox_inches='tight',pad_inches=0.03,format='eps')
                plt.close()
            else:
                fname_out = self.infile.replace('.fits','_lrg.png')
                fname_out = fname_out.replace('/data/','/plots/')
                fig.savefig(fname_out,bbox_inches='tight',pad_inches=0.03,format='png')
                plt.close()
        else:
            plt.tight_layout()
            plt.show(block=False)

#This is mainly for testing the functions within the class.
if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Plot one of the spectra in Lyke et al. 2020')
    parser.add_argument('t', choices=['s','l','p'], help='Plot size: (s)mall, (l)arge, (p)oster')
    parser.add_argument('-d', '--default_params', action='store_true',
                        help='Use a default parameter set for testing')
    parser.add_argument('-r', '--rest_top', action='store_true',
                        help='Mark top x-axis with rest frame wavelengths')
    parser.add_argument('-z', '--redshift', type=float, default=0.0, metavar='',
                        help='Redshift of spectrum for rest frame shift')
    parser.add_argument('-m', '--smooth', action='store_true',
                        help='Smooth spectrum source flux')
    parser.add_argument('-n', '--smooth_box', type=int, default=10, metavar='',
                        help='Smoothing window size in pixels, must be even')
    parser.add_argument('-w', '--weighted', action='store_true',
                        help='Use an error-weighted smoothing')
    parser.add_argument('-e', '--error', action='store_true',
                        help='Include the error spectrum in red')
    parser.add_argument('-k', '--sky', action='store_true',
                        help='Include the sky spectrum in green')
    parser.add_argument('-x', '--scale_sky', action='store_true',
                        help='Scale the sky flux to the source flux')
    parser.add_argument('-s', '--save', action='store_true',
                        help='Save the plot, without display')
    parser.add_argument('-f', '--file_type', choices=['png','eps'],default='png',
                        metavar='', help="File type for saving plot: 'eps' or 'png'")

    args = parser.parse_args()

    spec = Spectrum('../dr16q/data/spec-7672-57339-0394.fits')
    args.redshift = 3.303
    '''
    def plot_spec(self, plot_size, rest_twin=False, spec_redshift=1.0,
                    smooth=False, smooth_box=10, weighted_avg=False,
                    err=False, sky=False, scale_sky=False,
                    save=False, file_type='png')
    '''
    if args.default_params==True:
        args.rest_top = True
        args.smooth = True
        args.smooth_box = 10
        args.weighted = True
        args.error = True
    spec.plot_spec(args.t, rest_twin=args.rest_top, spec_redshift=args.redshift,
                    smooth=args.smooth, smooth_box=args.smooth_box,
                    weighted_avg=args.weighted, err=args.error, sky=args.sky,
                    scale_sky=args.scale_sky, save=args.save,
                    file_type=args.file_type)
