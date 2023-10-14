from astropy.io import fits
import numpy as np
from process_image import make_fits

def get_data_median(data_files:list,frame_type:str='science'):
    """
    Gets the median of a list of fits files.
    Writes the median to a fits file.
    Parameters:
    ----------
    data_files: list
        List of fits files.
    frame_type: str
        Type of frame. Default is 'science'.
    Returns:
    -------
    med_data: np.ndarray
        Median of the data.
    
    """
    data = [fits.open(f)[0].data for f in data_files]
    data = np.stack(data,axis=0)
    out = frame_type + '_median.fits'
    med_data = np.nanmedian(data,axis=0)
    make_fits(med_data,data_files[0],outfile=out)
    return np.nanmedian(data,axis=0)

#Read all the files
import os
path = os.getcwd()
files = os.listdir(path)

def spec_files(init:str='spec',files:list=files):
    """Returns the files that are spectra."""
    my_files = [f for f in files if f.startswith('b_flat')]
    work = get_data_median(my_files,frame_type=init)
    done_str = f'Processed {len(my_files)} of {init} files.'
    print(done_str)
    return work

    inits = ['b_flat','v_flat']
    for init in inits:
        print(spec_files(init))


def reduce_science_image(files:list=files):
    """Reduces the science image."""
    b_files = [f for f in files if f.startswith('b_science')]
    v_files = [f for f in files if f.startswith('v_science')]

    b_flat_med = spec_files(init='b_flat')
    v_flat_med = spec_files(init='v_flat')


    for f in b_files:
        data = fits.open(f)[0].data
        red_data = data/b_flat_med
        make_fits(red_data,f,outfile=f.replace('.fits','_reduced.fits'))
        print(f'Processed {f}')
    for f in v_files:
        data = fits.open(f)[0].data
        red_data = data/v_flat_med
        make_fits(red_data,f,outfile=f.replace('.fits','_reduced.fits'))
        print(f'Processed {f}')

    return None

#List all files ending in reduced.fits
def get_final_images(files:list=files):
    reduced_files = [f for f in files if f.endswith('reduced.fits')]
    #get median of reduced b_science files
    b_reduced = [f for f in reduced_files if f.startswith('b_science')]
    b_med = get_data_median(b_reduced,frame_type='b_science')

    #get median of reduced v_science files
    v_reduced = [f for f in reduced_files if f.startswith('v_science')]
    v_med = get_data_median(v_reduced,frame_type='v_science')



if __name__ == '__main__':
    get_final_images()



