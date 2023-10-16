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
    my_files = [f for f in files if f.startswith(init)]
    work = get_data_median(my_files,frame_type=init)
    done_str = f'Processed {len(my_files)} of {init} files.'
    print(done_str)
    return work

# Median of b_flat files
b_med = spec_files(init='b_flat')

# Median of v_flat files

v_med = spec_files(init='v_flat')

# Replace all zeros with 1e-6
b_med[b_med==0] = 1e-6
v_med[v_med==0] = 1e-6

# Divide science files by median of b_flat and v_flat
def get_reduced_image(init:str='b_science',files:list=files,med_data=None):
    """Returns the reduced image."""
    my_files = [f for f in files if f.startswith(init)]
    print(my_files)
    init_data = [fits.open(f)[0].data/med_data for f in my_files]
    data = np.stack(init_data,axis=0)
    reduced = np.nanmedian(data,axis=0)
    out_filename = init + '_reduced.fits'
    out = os.path.join(path,'reduced_files',out_filename)
    make_fits(reduced,my_files[0],outfile=out)
    print(f'Reduced {init} files.')
    return None

if __name__ == '__main__':
    get_reduced_image(init='b_science',med_data=b_med)
    get_reduced_image(init='v_science',med_data=v_med)
    get_reduced_image(init='b_cluster',med_data=b_med)
    get_reduced_image(init='v_cluster',med_data=v_med)