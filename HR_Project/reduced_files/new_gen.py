import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils import DAOStarFinder, CircularAperture, aperture_photometry
from astropy.visualization import simple_norm

# Load image data
def load_image(filename):
    with fits.open(filename) as hdul:
        data = hdul[0].data
    return data

b_image_starfield = load_image("b_science_reduced.fits")
v_image_starfield = load_image("v_science_reduced.fits")
b_image_cluster = load_image("b_cluster_reduced.fits")
v_image_cluster = load_image("v_cluster_reduced.fits")

# Detect stars
def detect_stars(data):
    mean, median, std = np.mean(data), np.nanmedian(data), np.std(data)
    daofind = DAOStarFinder(fwhm=10.0, threshold=2.5*std, 
                            sharplo=0.2, sharphi=1.0, 
                            roundlo=-1.0, roundhi=1.0)
    sources = daofind(data - median)
    return sources


sources_starfield = detect_stars(b_image_starfield)
print(sources_starfield)
# # Define apertures and perform photometry
# apertures_starfield = CircularAperture([sources_starfield['xcentroid'], sources_starfield['ycentroid']], r=6.)
# rawflux_b_starfield = aperture_photometry(b_image_starfield, apertures_starfield)
# rawflux_v_starfield = aperture_photometry(v_image_starfield, apertures_starfield)

# sources_cluster = detect_stars(b_image_cluster)
# apertures_cluster = CircularAperture([sources_cluster['xcentroid'], sources_cluster['ycentroid']], r=6.)
# rawflux_b_cluster = aperture_photometry(b_image_cluster, apertures_cluster)
# rawflux_v_cluster = aperture_photometry(v_image_cluster, apertures_cluster)

# # Calculate magnitude using the reference star
# ref_flux_b = rawflux_b_starfield['aperture_sum'][0]  # Assuming the reference star is the first detected
# ref_flux_v = rawflux_v_starfield['aperture_sum'][0]
# m_b_ref = 12.1
# m_v_ref = 10.832

# def calculate_magnitude(rawflux, ref_flux, m_ref):
#     return -2.5 * np.log10(rawflux / ref_flux) + m_ref

# m_b_cluster = calculate_magnitude(rawflux_b_cluster['aperture_sum'], ref_flux_b, m_b_ref)
# m_v_cluster = calculate_magnitude(rawflux_v_cluster['aperture_sum'], ref_flux_v, m_v_ref)

# # Plot HR Diagram
# colors = m_b_cluster - m_v_cluster
# plt.scatter(colors, m_v_cluster, color='blue', s=50)
# plt.gca().invert_yaxis()
# plt.xlabel('B-V')
# plt.ylabel('V')
# plt.title('HR Diagram')
# plt.show()

# # Plot detected stars on the images
# def plot_stars_on_image(data, apertures):
#     norm = simple_norm(data, 'sqrt', percent=99)
#     plt.imshow(data, norm=norm, cmap='Greys_r', origin='lower')
#     apertures.plot(color='red', lw=1.5)
#     plt.show()

# plot_stars_on_image(b_image_starfield, apertures_starfield)
# plot_stars_on_image(b_image_cluster, apertures_cluster)
