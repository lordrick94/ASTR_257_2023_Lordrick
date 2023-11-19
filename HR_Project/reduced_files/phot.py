"""
This script is used to find the stars in the image and make a H-R diagram.

Algorithm:
----------
1. Grab the stars from the image and calculate their aperture sum.

"""

import numpy as np
import pandas as pd

from astropy.io import fits

from photutils.aperture import aperture_photometry, CircularAperture
from photutils.datasets import load_star_image
from photutils.detection import DAOStarFinder
from astropy.stats import mad_std


def grab_stars(file:str='b_science_reduced.fits'):
    """
    This function takes in a file and returns a dataframe with the x and y coordinates of the stars in the image.
    Uses photutils DAOStarFinder with aperture photometry to find the stars in the image.

    Parameters
    ----------
    file : str
        The name of the file to be used.
    Returns
    -------
    df : pandas dataframe
        A dataframe with the x, y coordinates and the aperture sum of the stars in the image.
    """
    F_h,bkg_factor,radius = 10.0, 2.5, 2.0
    hdu = fits.open(file)
    image = hdu[0].data

    image_1 = image -  np.nanmedian(image)

    bkg_sigma = mad_std(image_1)
    daofind = DAOStarFinder(fwhm=F_h, threshold= bkg_factor*bkg_sigma)  
    sources = daofind(image_1)
    for col in sources.colnames:  
        sources[col].info.format = '%.8g'  # for consistent table output 

    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))  
    apertures = CircularAperture(positions, r=radius)  
    phot_table = aperture_photometry(image, apertures)  
    for col in phot_table.colnames:  
        phot_table[col].info.format = '%.8g'  # for consistent table output
    df = phot_table.to_pandas()

    # Sort the dataframe by aperture sum
    df = df.sort_values(by=['aperture_sum'], ascending=False)
    df = df.reset_index(drop=True)
    return df

x,y = 200, 510
b_science_stars = grab_stars('b_science_reduced.fits')
df = b_science_stars
df_star = df[(df['xcenter']>x) & 
   (df['xcenter']<x+20) & (df['ycenter']>y) & (
       df['ycenter']<y+20)]

L_b_star_sum = df['aperture_sum'][0].astype(float)

m_ref = 10.382+1.959 #Known star
F_ref = L_b_star_sum
df_b_cluster = grab_stars('b_cluster_reduced.fits')

F_1 = df_b_cluster['aperture_sum'].values
m_B = m_ref -2.5*np.log10(F_1/F_ref)

df_b_cluster['m_B'] = m_B

df_b_cluster['m_V'] = 100

# Remove Nan Rows
df_b_cluster = df_b_cluster[df_b_cluster['m_B'].notna()]

v_science_stars = grab_stars('v_science_reduced.fits')
df = v_science_stars
df_star = df[(df['xcenter']>x) & 
   (df['xcenter']<x+20) & (df['ycenter']>y) & (
       df['ycenter']<y+20)]

L_v_star_sum = df['aperture_sum'][0].astype(float)


m_ref = 10.382 #Known star
F_ref = L_v_star_sum
df_v_cluster = grab_stars('v_cluster_reduced.fits')

F_1 = df_v_cluster['aperture_sum'].values
m_V = m_ref -2.5*np.log10(F_1/F_ref)
df_v_cluster['m_V'] = m_V
df_v_cluster['m_B'] = 100
# Remove Nan Rows
df_v_cluster = df_v_cluster[df_v_cluster['m_V'].notna()]


# Cross match the stars
df_b_cluster = df_b_cluster.reset_index(drop=True)
df_v_cluster = df_v_cluster.reset_index(drop=True)

df_b_cluster = df_b_cluster.rename(columns={'xcenter':'xcenter_b', 'ycenter':'ycenter_b'})
df_v_cluster = df_v_cluster.rename(columns={'xcenter':'xcenter_v', 'ycenter':'ycenter_v'})

a = df_b_cluster.copy()
b = df_v_cluster.copy()

tolerance = 2
# Cartesian product of both dataframes
a['key'] = 1
b['key'] = 1

cartesian = pd.merge(a, b, on='key').drop('key', axis=1)

# Compute the distance between points
cartesian['dist_x'] = (cartesian['xcenter_b'] - cartesian['xcenter_v']).abs()
cartesian['dist_y'] = (cartesian['ycenter_b'] - cartesian['ycenter_v']).abs()

# Filter using the tolerance
filtered = cartesian[(cartesian['dist_x'] <= tolerance) & (cartesian['dist_y'] <= tolerance)]

# If there are duplicates, keep the closest match (in terms of x-distance)
filtered = filtered.sort_values(by=['xcenter_b', 'ycenter_b', 'dist_x', 'dist_y']).drop_duplicates(subset=['xcenter_b', 'ycenter_b'])

# Drop unnecessary columns
filtered = filtered.drop(columns=['dist_x', 'dist_y'])

df_final = filtered[['m_B_x', 'm_V_y']]
print(df_final)

# Plot the H-R diagram
import matplotlib.pyplot as plt

# Data
x = df_final['m_B_x'] - df_final['m_V_y']
y = df_final['m_V_y']

# Create figure and axis
fig, ax = plt.subplots(figsize=(10, 6))

# Scatter plot with customization
ax.scatter(x, y, c='crimson', alpha=0.6, edgecolors='w', linewidth=0.5, s=50, marker="o")
ax.set_facecolor('lightgrey')

# Labels, title, and grid
ax.set_title('H-R Diagram for NGC7618', fontsize=16, fontweight='bold')
ax.set_xlabel('$m_B - m_V$', fontsize=14)
ax.set_ylabel('$m_V$', fontsize=14)
ax.set_ylim(ax.get_ylim()[::-1])
ax.grid(True, linestyle='--', which='both', linewidth=0.5, alpha=0.7)

# Tight layout
plt.tight_layout()

# Display the plot
plt.show()

#Save the plot
fig.savefig('hr_diagram.png', dpi=400)