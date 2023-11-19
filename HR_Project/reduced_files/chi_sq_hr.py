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


def find_stars(image, F_h=15.0, bkg_factor=2):
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


def compute_chisq(observed_magnitudes, isochrone_magnitudes):
    """Compute the chi-squared value between the observed data and isochrone."""
    chisq = np.sum((observed_magnitudes - isochrone_magnitudes) ** 2)
    return chisq


def find_nearest(array, value):
    """Find the nearest value in an array and return the index."""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def fit_isochrone(df_b, df_v, isochrone_df):
    """Fit the isochrone to the data using chi-squared statistic."""
    x_data = df_b['m_b'] - df_v['m_v']
    y_data = df_v['m_v']

    x_iso = (isochrone_df['Bmag'].astype(float) - isochrone_df['Vmag'].astype(float)).values
    y_iso =( isochrone_df['Vmag'].astype(float)).values

    chisq_values = []
    for y in y_data:
        idx = find_nearest(y_iso, y)
        chisq_values.append((y - y_iso[idx]) ** 2)

    chisq = sum(chisq_values)
    return chisq


def plot_hr_diagram(df_b, df_v, isochrone_df):
    x_data = df_b['m_b'] - df_v['m_v']
    y_data = df_v['m_v']

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(x_data, y_data, c='crimson', alpha=0.6, edgecolors='w', linewidth=0.5, s=50)
    ax.set_facecolor('lightgrey')
    ax.set_title('H-R Diagram for NGC7618', fontsize=16, fontweight='bold')
    ax.set_xlabel('$m_B - m_V$', fontsize=14)
    ax.set_ylabel('$m_V$', fontsize=14)
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.grid(True, linestyle='--', which='both', linewidth=0.5, alpha=0.7)

    # Plot isochrone if provided
    if isochrone_df is not None:
        x_iso = isochrone_df['Bmag'].astype(float) - isochrone_df['Vmag'].astype(float)
        y_iso = isochrone_df['Vmag'].astype(float) + 13.5
        ax.plot(x_iso, y_iso, color='blue', label='Isochrone')
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

    iso = pd.read_csv('isochrone.dat', delim_whitespace=True, header=None) 
    iso.columns = iso.iloc[0]
    # Drop first row
    iso = iso.drop(iso.index[0])

    chisq = fit_isochrone(df_b_cluster, df_v_cluster, iso)
    print(f"Chi-squared value: {chisq}")

    plot_hr_diagram(df_b_cluster, df_v_cluster, iso)


if __name__ == "__main__":
    main()
