from data_processing import prepare_catalog
from companion_finder import find_companions
import numpy as np

# Example paths (user must modify)
optap_path = "data/CEERS_optap.fits"
photoz_path = "data/CEERS_photoz.fits"

filtered_catalog, red_mask, blue_mask, F160W_mag, F444W_mag, color = prepare_catalog(
    optap_path,
    photoz_path
)

z_all = np.array(filtered_catalog['zphot'])

# Dummy binning (can be expanded)
z_edges = np.linspace(np.min(z_all), np.max(z_all), 10)

red_results = find_companions(
    filtered_catalog,
    red_mask,
    z_all,
    z_edges,
    F160W_mag,
    F444W_mag,
    color,
    "red"
)

print("Red galaxies with companions:", len(red_results))
