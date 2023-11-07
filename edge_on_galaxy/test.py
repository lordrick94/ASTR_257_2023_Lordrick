# Astropy imports
from astropy.io import fits

# Other imports
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt


path = '/home/lordrick/Documents/ASTRO_PHD_UCSC/Year_1_Astro_PHD/ASTR_257/'

# Bias r1014 - r1020
def get_bias():
    bias_fits = [fits.open(os.path.join(path,'Spectro_KAST','data',f'r{i}.fits'))[0].data for i in range(1014,1021)]

    bias_med = np.stack(bias_fits,axis=0)
    bias = np.nanmedian(bias_med,axis=0)

    return bias

def row_normalize(array_2d):
    """
    Normalizes each row of a 2d array by dividing each row by the median of that row.

    Parameters:
    -----------
    array_2d : 2d array
        Array to normalize
    """
    for row in range(len(array_2d)):
        new_row = array_2d[row]/np.median(array_2d[row])
        array_2d[row] = new_row
    return array_2d

def get_fits_data(n1,n2,norm=False):
    """
    Returns the data from a fits file.

    Parameters:
    -----------
    n1 : int
        Starting number of fits file.
    n2 : int
        Ending number of fits file.
    ax : int
        Axis along which to take the median when normalizing.
    norm : bool
        Normalize the data
    
    """

    bias = get_bias()

    if norm:
        hdus  = [fits.open(os.path.join(path,
                'Spectro_KAST','data',f'r{i}.fits'))[0].data 
                - bias for i in range(n1,n2+1)]
        
        # Normalize by dividing each row by the median of that row
        hdus_norm = [row_normalize(hdu) for hdu in hdus]

        hdus = np.stack(hdus_norm,axis=0)
        
        hdu_data = np.nanmedian(hdus_norm,axis=0)

    else:
        hdus  = [fits.open(os.path.join(path,
                'Spectro_KAST','data',f'r{i}.fits'))[0].data 
                 - bias for i in range(n1,n2+1)]
        hdus = np.stack(hdus,axis=0)
        hdu_data = np.nanmedian(hdus,axis=0)

    return hdu_data

# Arcs r1000 - r1002
arcs_med = get_fits_data(1000,1002)

# Flats r1003 - r1013
flats_med_norm = get_fits_data(1003,1013,norm=True)

# UGC 18876 r1022 - r1024
UGC_med = get_fits_data(1025,1030)

# Science files
UGC_final = (UGC_med/300)/(flats_med_norm/3)
arcs_final = arcs_med/(flats_med_norm)


def make_fits(data,filename):
    """Makes a fits file from data."""
    hdr = fits.getheader(os.path.join(path,'Spectro_KAST','data','r1022.fits'))
    hdu = fits.PrimaryHDU(data=data,header=hdr)
    hdu.writeto(filename,overwrite=True)

#make_fits(UGC_final,'UGC_18876.fits')
_ = make_fits(UGC_final,'UGC.fits')
_ = make_fits(arcs_final,'arcs.fits')
