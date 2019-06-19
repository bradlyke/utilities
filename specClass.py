#Import the tools necessary to work with spectra files.
import cat_tools as ct #This is a tool I wrote, so you'll need the script in your path.
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as ticker
matplotlib.rc('text',usetex=True)

#A spectrum is a class. The object has 5 attributes right now: self, boxcar,
#paper_plot_small, plot_large, and poster_plot.
class spectrum:
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
    def boxcar(self,flux_arr,flux_var,smooth_pct,weighted=False):
        spc = smooth_pct
        num_dpoints = len(flux_arr)
        box_flux = np.zeros(num_dpoints,dtype='f8')
        box_err = np.zeros(num_dpoints,dtype='f8')

        for i in range(num_dpoints):
            #Define range of values to smooth
            if i < int(spc/2):
                lower = 0
            else:
                lower = -int(spc/2) + i
            if ((int(spc/2) + i) > num_dpoints):
                upper = num_dpoints
            else:
                upper = int(spc/2) + i
            #Smooth the values depending on smoothing type
            if weighted==True:
                noise_temp = np.sqrt(flux_var[lower:upper])
                signal_temp = flux_arr[lower:upper]
                flux_temp = np.average(signal_temp,weights=noise_temp)
                ivar_temp = (np.average((signal_temp-flux_temp)**2, weights=noise_temp))**(-1.0)
            else:
                flux_temp = np.median(flux_arr[lower:upper])
                ivar_temp = flux_var[i]

            box_flux[i] = flux_temp
            box_err[i] = ivar_temp

        return box_flux,box_err


    def paper_plot_small(self,smooth=False,err=False,sky=False,scale_sky=False,save=False,form='p'):
        wobs = np.where((self.lam>=3700)&(self.lam<=9000))[0]
        flux_range = np.amax(self.flux[wobs]) - np.amin(self.flux[wobs])
        flux_pad = float(flux_range)/10
        y_lower = np.amin(self.flux[wobs]) - flux_pad
        y_upper = np.amax(self.flux[wobs]) + flux_pad
        x_lower = np.amin(self.lam)
        x_upper = np.amax(self.lam)

        matplotlib.rc('font',size=15)
        fig1,ax1 = plt.subplots(figsize=(5,4))
        if smooth==True:
            self.bflux,self.berr = self.boxcar(self.flux,self.ivar,10,weighted=True)
            ax1.plot(self.lam,self.flux,color='0.70',linewidth=0.8)
            ax1.plot(self.lam,self.bflux,color='black',linewidth=0.6)
        else:
            ax1.plot(self.lam,self.flux,color='black')
        if sky==True:
            if scale_sky==True:
                self.sky = self.data['sky'] / 10.0
            else:
                self.sky = self.data['sky']
            ax1.plot(self.lam,self.sky,color='green',linewidth=0.7,alpha=0.5)
        if err==True:
            ax1.plot(self.lam,self.ferr,color='red',linewidth=0.6)
        ax1.set_xlim((x_lower,x_upper))
        ax1.set_ylim((y_lower,y_upper))
        ax1.set_xlabel(r'Wavelength (\AA)',fontsize=15)
        ax1.set_ylabel(r'$f_{\lambda}$ ($10^{-17}$ ergs s$^{-1}$ cm$^{-2}$\,\AA$^{-1}$)',fontsize=15)
        ax1.tick_params(axis='both',direction='in')
        ax1.tick_params(axis='both',which='minor',direction='in')
        ax1.tick_params(top=True,right=True)
        ax1.tick_params(which='minor',top=True,right=True)
        ax1.xaxis.set_minor_locator(ticker.MultipleLocator(100))
        ax1.yaxis.set_minor_locator(ticker.MultipleLocator(1))
        matplotlib.rc('font',size=15)

        plt.tight_layout()
        if save==True:
            if form=='e':
                fname_out = self.infile.replace('.fits','_sml.eps')
                fig1.savefig(fname_out,format='eps')
            else:
                fname_out = self.infile.replace('.fits','_sml.png')
                fig1.savefig(fname_out,format='png')
        else:
            plt.show()

    def plot_large(self,smooth=False,err=False,sky=False,scale_sky=False,save=False,form='p'):
        wobs = np.where((self.lam>=3700)&(self.lam<=9000))[0]
        flux_range = np.amax(self.flux[wobs]) - np.amin(self.flux[wobs])
        flux_pad = float(flux_range)/10
        y_lower = np.amin(self.flux[wobs]) - flux_pad
        y_upper = np.amax(self.flux[wobs]) + flux_pad
        x_lower = np.amin(self.lam)
        x_upper = np.amax(self.lam)

        matplotlib.rc('font',size=18)
        fig1,ax1 = plt.subplots(figsize=(10,8))
        if smooth==True:
            self.bflux,self.berr = self.boxcar(self.flux,self.ivar,10,weighted=True)
            ax1.plot(self.lam,self.flux,color='0.70',linewidth=0.8)
            ax1.plot(self.lam,self.bflux,color='black',linewidth=0.6)
        else:
            ax1.plot(self.lam,self.flux,color='black')
        if sky==True:
            if scale_sky==True:
                self.sky = self.data['sky'] / 10.0
            else:
                self.sky = self.data['sky']
            ax1.plot(self.lam,self.sky,color='green',linewidth=0.7,alpha=0.5)
        if err==True:
            ax1.plot(self.lam,self.ferr,color='red',linewidth=0.6)
        ax1.set_xlim((x_lower,x_upper))
        ax1.set_ylim((y_lower,y_upper))
        ax1.set_xlabel(r'Wavelength (\AA)',fontsize=18)
        ax1.set_ylabel(r'$f_{\lambda}$ ($10^{-17}$ ergs s$^{-1}$ cm$^{-2}$\,\AA$^{-1}$)',fontsize=18)
        ax1.tick_params(axis='both',direction='in')
        ax1.tick_params(axis='both',which='minor',direction='in')
        ax1.tick_params(top=True,right=True)
        ax1.tick_params(which='minor',top=True,right=True)
        ax1.xaxis.set_minor_locator(ticker.MultipleLocator(100))
        ax1.yaxis.set_minor_locator(ticker.MultipleLocator(1))
        matplotlib.rc('font',size=18)

        plt.tight_layout()
        if save==True:
            if form=='e':
                fname_out = self.infile.replace('.fits','_big.eps')
                fig1.savefig(fname_out,format='eps')
            else:
                fname_out = self.infile.replace('.fits','_big.png')
                fig1.savefig(fname_out,format='png')
        else:
            plt.show()

    def poster_plot(self,smooth=False,err=False,sky=False,scale_sky=False,save=False,form='p'):
        wobs = np.where((self.lam>=3700)&(self.lam<=9000))[0]
        flux_range = np.amax(self.flux[wobs]) - np.amin(self.flux[wobs])
        flux_pad = float(flux_range)/10
        y_lower = np.amin(self.flux[wobs]) - flux_pad
        y_upper = np.amax(self.flux[wobs]) + flux_pad
        x_lower = np.amin(self.lam)
        x_upper = np.amax(self.lam)

        matplotlib.rc('font',size=24)
        fig1,ax1 = plt.subplots(figsize=(15,12))
        if smooth==True:
            self.bflux,self.berr = self.boxcar(self.flux,self.ivar,10,weighted=True)
            ax1.plot(self.lam,self.flux,color='0.70',linewidth=0.8)
            ax1.plot(self.lam,self.bflux,color='black',linewidth=0.6)
        else:
            ax1.plot(self.lam,self.flux,color='black')
        if sky==True:
            if scale_sky==True:
                self.sky = self.data['sky'] / 10.0
            else:
                self.sky = self.data['sky']
            ax1.plot(self.lam,self.sky,color='green',linewidth=0.7,alpha=0.5)
        if err==True:
            ax1.plot(self.lam,self.ferr,color='red',linewidth=0.6)
        ax1.set_xlim((x_lower,x_upper))
        ax1.set_ylim((y_lower,y_upper))
        ax1.set_xlabel(r'Wavelength (\AA)',fontsize=24)
        ax1.set_ylabel(r'$f_{\lambda}$ ($10^{-17}$ ergs s$^{-1}$ cm$^{-2}$\,\AA$^{-1}$)',fontsize=24)
        ax1.tick_params(axis='both',direction='in')
        ax1.tick_params(axis='both',which='minor',direction='in')
        ax1.tick_params(top=True,right=True)
        ax1.tick_params(which='minor',top=True,right=True)
        ax1.xaxis.set_minor_locator(ticker.MultipleLocator(100))
        ax1.yaxis.set_minor_locator(ticker.MultipleLocator(1))
        matplotlib.rc('font',size=24)

        plt.tight_layout()
        if save==True:
            if form=='e':
                fname_out = self.infile.replace('.fits','_poster.eps')
                fig1.savefig(fname_out,format='eps')
            else:
                fname_out = self.infile.replace('.fits','_poster.png')
                fig1.savefig(fname_out,format='png')
        else:
            plt.show()

#This is mainly for testing the functions within the class.
#if __name__=='__main__':
    #spec = spectrum('spec-7379-56713-0938-dr16.fits')
    #spec.paper_plot_small(smooth=True,err=True,sky=True,scale_sky=True)
    #spec.plot_large(smooth=True,err=True,sky=True,scale_sky=True)
    #spec.poster_plot(smooth=True,err=True,sky=True,scale_sky=True)
