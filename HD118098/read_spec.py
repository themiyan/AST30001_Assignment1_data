"""
python snippet to plot the MOOG synthetic spectra with the observed spectra.
email issues to: tnanayak@astro.swin.edu.au
-Themiya 09/08/2016
"""



###############DEFINE PARAMETERS BELOW###############


#enter the file names below. if the files are not in the same directory give the relative paths
#enter the name of the star file
star_name ='HD118098.fits'
#enter the name of the synthetic spectra file
syn_name='CH.out2'

#enter the radial velocity of the star in km/s
#REMEMBER: moving towards the Earth is +, moving away is -
vr= -24.67



#define zoom_range if you want to zoom in the plot to inspect features
#should be defined in a python list as min and max of x first and then min and max of y
#eg: 
#zoom_range=[xmin,xmax,ymin,ymax] 
#where xmin, xmax, ymin, ymax are integer/float values defined by the user 
zoom_range=None

#########################################################

def main():

    syn_wave, syn_flux  = open_synthetic_spectra(syn_name)
    obs_wave, obs_flux = open_observed_spectra(star_name)
    
    
    obs_corrected_wave = redshift_correction(obs_wave,vr )
    
    plot_spectra(star_name, obs_wave, obs_flux, syn_wave, syn_flux, zoom_range=zoom_range)






def open_synthetic_spectra(file_name): 
    """
    Open the synthetic spectra. 
    The code will skip the first 3 lines of the file and read everything else. 
    If the last column of the file doesn't have the same number of columns as the 
    rest of the spectra (number of columns in the 4th row onwards),
    the last line will be skipped. 
    
    INPUT: 
        file_name: file name of the synthetic spectra with extension
    
    RETURN: 
        synthetic_wavelength, synthetic_flux: sythetic spectra wavelength, flux
    

    """
    
    
    line = open(file_name, "r").readlines()[2]

    params = np.asarray(s.split(line,' '))
    params = np.asarray(params[np.asarray(params)!=''], dtype='float')

    min_wave   = params[0]
    max_wave   = params[1]
    delta_wave = params[2]
    
    synthetic_flux       = np.genfromtxt(file_name, skip_header=3, invalid_raise=False)
    synthetic_flux       = 1-(synthetic_flux.flatten())
    synthetic_wavelength = np.arange(min_wave, max_wave, delta_wave)

    
    return synthetic_wavelength, synthetic_flux


def open_observed_spectra(file_name):
    """
    Open the observed spectra from fits format. 
    
    The wavelength is computed using fits header information. 

    wavelength = ((x_pixel_number + 1.0) - CRPIX1) * CDELT1 + CRVAL1 

    INPUT: 
        file_name: file name of the observed spectra with extension
    
    RETURN: 
        observed_wave, observed_flux: observed spectra wavelength, flux
    
    
    """
    
    observed_fits = fits.open(file_name)
    
    observed_flux   = observed_fits[0].data
    observed_header = observed_fits[0].header 
    
    observed_wave   = ((np.arange(len(observed_flux)) + 1.0) - observed_header['CRPIX1']) * observed_header['CDELT1'] + observed_header['CRVAL1'] 
    
    return observed_wave, observed_flux
    




def redshift_correction(wavelength, velocity):
    """
    Corrects the observed spectra for radial velocities
    
    INPUT: 
        wavelength: the observed wavelength of the spectra
        velocity: the radial velocity of the star in km/s
        
    RETURN: 
        wave_reast: the wavelength corrected for radial velocity of the star
    
    """
    
    c=3e5 #speed of light in Km/s
    
    wave_rest = wavelength/(1 + (velocity/c))
    
    return wave_rest
    
    
def plot_spectra(star_name,observed_wave, observed_flux, synthetic_wavelength, synthetic_flux, zoom_range=None):
    """
    Plot the observed and synthesized spectra.
    Saves the figure in the same directory as the code. 
    
    INPUT: 
        star_name: Name of the star. The file name extension will be removed.
        observed_wave, observed flux: observed wavelength and flux of the spectra
        synthetic_wavelength, synthetic_flux: wavelength and flux of the synthesized spectra
        
    OPTIONAL: 
        zoom_range: the zoom pange of the plot. if set to None will use the min and max parameters
                    from the synthetic spectra for the x axis. y axis is set for 0.05 to 1.05. 
                    defined as zoom_range=[xmin, xmax, ymin,ymax]
                    
    RETURN: None
    
    
    """
    
    star_name= s.split(star_name, sep='.')[0]
    
    fig, ax = plt.subplots(figsize=(10,4))

    ax.step(observed_wave, observed_flux, color='navy', label=r'Observed spectra')

    ax.step(synthetic_wavelength, synthetic_flux, color='deeppink', label=r'MOOG Synthetic spectra', alpha=0.5)

    plt.title(r''+str(star_name), fontsize=20)
    plt.xlabel(r'Wavelength (\AA)', fontsize=15)
    plt.ylabel(r'Flux (normalized)', fontsize=15)

    plt.legend(loc='best', fontsize='large')

    if zoom_range==None:
        plt.xlim(np.min(synthetic_wavelength), np.max(synthetic_wavelength))
        plt.ylim(0.05,1.05)
    else:
        plt.xlim(zoom_range[0], zoom_range[1])
        plt.ylim(zoom_range[2], zoom_range[3])
        
    plt.tight_layout()
    plt.savefig('obs_syn_comp_'+str(star_name)+'.png', dpi=400)
    
    print "Figure saved to disk" 

    plt.close()

    



if __name__ == '__main__':
    """Main entry point to the code"""
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import fits
    import string as s
    from matplotlib import rc
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)

			
    main()
    



