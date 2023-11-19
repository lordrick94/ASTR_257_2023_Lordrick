import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils.aperture import aperture_photometry, CircularAperture
from photutils.detection import DAOStarFinder
from astropy.stats import mad_std

def load_fits_image(file):
    hdu = fits.open(file)
    image = hdu[0].data
    return image - np.nanmedian(image)


def load_iso(file='isochrone.dat'):
    iso = pd.read_csv(file, delim_whitespace=True, header=None) 
    iso.columns = iso.iloc[0]
    #Drop first row
    iso = iso.drop(iso.index[0])

    #Adjusting isochrone
    del_mag_x,del_mag_y = 11.9,12
    iso['Vmag'] = iso['Vmag'].astype(float) + del_mag_x
    iso['Bmag'] = iso['Bmag'].astype(float) + del_mag_y

    #Sort by Vmag
    iso = iso.sort_values(by=['Vmag'], ascending=False).reset_index(drop=True)

    return iso


def find_stars(image, F_h=15.0, bkg_factor=1):
    bkg_sigma = mad_std(image)
    daofind = DAOStarFinder(fwhm=F_h, threshold=bkg_factor*bkg_sigma)
    sources = daofind(image)
    for col in sources.colnames:
        sources[col].info.format = '%.8g'
    return sources


def photometry_data(image, sources, radius=15):
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r=radius)
    phot_table = aperture_photometry(image, apertures)
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'
    df = phot_table.to_pandas()
    return df.sort_values(by=['aperture_sum'], ascending=False).reset_index(drop=True)


def magnitude_data(df, m_ref, f_ref, f):
    df['m_' + f] = m_ref - 2.5 * np.log10(df['aperture_sum'].values / f_ref)
    return df


def plot_hr_diagram(df_b, df_v, iso_df=None):
    x_data = df_b['m_b'] - df_v['m_v']
    y_data = df_v['m_v']

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(x_data, y_data, c='crimson', alpha=0.6, edgecolors='w', linewidth=0.5, s=100)
    ax.set_facecolor('lightgrey')
    ax.set_title('H-R Diagram for NGC6819', fontsize=20, fontweight='bold')
    ax.set_xlabel('$m_B - m_V$', fontsize=20)
    ax.set_ylabel('$m_V$', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.set_xlim(0, 2.5)
    ax.grid(True, linestyle='--', which='both', linewidth=0.5, alpha=0.7)

    # Plot isochrone if provided
    labels = ['Age = 1 Gyr', 'Age = 2 Gyr']
    if iso_df is not None:
        # Loop over isochrones and labels

        for isochrone_df, label in zip(iso_df, labels):
            x_iso = isochrone_df['Bmag'].astype(float) - isochrone_df['Vmag'].astype(float)
            y_iso = isochrone_df['Vmag'].astype(float)
            ax.scatter(x_iso, y_iso,s=20,label=label)
    
    ax.legend()

    plt.tight_layout()
    plt.show()
    fig.savefig('hr_diagram.png', dpi=400)

def main():
    # Load images
    b_image = load_fits_image('b_science_reduced.fits')
    v_image = load_fits_image('v_science_reduced.fits')

    # Find stars in images
    b_sources = find_stars(b_image)
    v_sources = find_stars(v_image)

    # Photometry data
    df_b = photometry_data(b_image, b_sources)
    df_v = photometry_data(v_image, v_sources)

    # Known stars and reference values
    m_ref_b, F_ref_b = 12.341, df_b['aperture_sum'][0].astype(float)
    m_ref_v, F_ref_v = 10.382, df_v['aperture_sum'][0].astype(float)

    # Magnitude data
    df_b_cluster = magnitude_data(photometry_data(load_fits_image('b_cluster_reduced.fits'), b_sources), m_ref_b, F_ref_b, 'b')
    df_v_cluster = magnitude_data(photometry_data(load_fits_image('v_cluster_reduced.fits'), v_sources), m_ref_v, F_ref_v, 'v')

    #Load isochrones
    iso = [load_iso(f'age{i}.dat') for i in range(1,3)]
    
    plot_hr_diagram(df_b_cluster, df_v_cluster, iso)



if __name__ == "__main__":
    main()
