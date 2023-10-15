from astropy.io import fits
import numpy as np
import os


data_dir = '/home/lordrick/Documents/ASTRO_PHD_UCSC/Year_1_Astro_PHD/ASTR_257/data_Eros_cluster'

def get_bias():
    bias_fits = [fits.open(os.path.join(data_dir,f'd{i}.fits'))[0].data for i in range(2000,2002+1)]
    bias_med = np.stack(bias_fits,axis=0)
    bias = np.nanmedian(bias_med,axis=0)
    return bias

#print(get_bias())

def make_fits(data,filename,outfile=None):
    """Makes a fits file from data."""
    hdr = fits.getheader(filename)
    hdr['HISTORY'] = 'Reduced using reduction.py'
    hdu = fits.PrimaryHDU(data=data,header=hdr)

    if outfile is None:
        outfile = filename.replace('.fits','_reduced.fits')
    hdu.writeto(outfile,overwrite=True)

def process_image(filename, sub_bias:bool=True, normalize:bool=False,per_time:bool=False,outfile=None):
    """Processes an image."""
    bias = get_bias()
    image = fits.open(filename)[0].data
    hdr = fits.getheader(filename)
    if sub_bias:
        image = image - bias
    
    if normalize:
        image = image/np.nanmedian(image)

    
    if per_time:
        try:
            exp_time = hdr['EXPTIME']
            image = image/exp_time
        except KeyError:
            print('No exposure time found in header.')
            Warning('No exposure time found in header, using default value exp_time=1.')
            exp_time = 1
            image = image/exp_time

    # Make fits file
    _ = make_fits(image,filename,outfile=outfile)
    print(f'Processed {filename}')
    
    return None

def process_all():
    # Process darks d2005.fits to d2005.fits
    for i in range(2005,2006):
        process_image(os.path.join(data_dir,f'd{i}.fits'),sub_bias=True,normalize=False,per_time=True,outfile=f'dark_d{i}.fits')

    # Process b_flats d2058,d2059,d2061
    for i in [2058,2059,2061]:
        process_image(os.path.join(data_dir,f'd{i}.fits'),sub_bias=True,normalize=True,per_time=False,outfile=f'b_flat_d{i}.fits')

    # Process v_flats d2062,d2063,d2064
    for i in [2062,2063,2064]:
        process_image(os.path.join(data_dir,f'd{i}.fits'),sub_bias=True,normalize=True,per_time=False,outfile=f'v_flat_d{i}.fits')

    # Process v_science images d2096 to d2098
    for i in range(2096,2098+1):
        process_image(os.path.join(data_dir,f'd{i}.fits'),sub_bias=True,normalize=False,per_time=True,outfile=f'v_science_d{i}.fits')

    # Process b_science images d2099 to d2103
    for i in range(2099,2103+1):
        process_image(os.path.join(data_dir,f'd{i}.fits'),sub_bias=True,normalize=False,per_time=True,outfile=f'b_science_d{i}.fits')

    # Process b_cluster images d2124 to d2128
    for i in range(2124,2128+1):
        process_image(os.path.join(data_dir,f'd{i}.fits'),sub_bias=True,normalize=False,per_time=True,outfile=f'b_cluster_d{i}.fits')

    # Process v_cluster images d2121 to d2123
    for i in range(2121,2123+1):
        process_image(os.path.join(data_dir,f'd{i}.fits'),sub_bias=True,normalize=False,per_time=True,outfile=f'v_cluster_d{i}.fits')
if __name__ == '__main__':
    process_all()


