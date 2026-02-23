import numpy as np
from astropy.table import Table, join

def prepare_catalog(optap_path, photoz_path):

    catalog = Table.read(optap_path)
    redshift_catalog = Table.read(photoz_path)

    merged_catalog = join(catalog, redshift_catalog, keys='ID', join_type='inner')

    valid_mask = (
        (merged_catalog['f_f444w'] / merged_catalog['e_f444w'] >= 3) &
        (merged_catalog['f_f160w'] > 0) &
        (merged_catalog['f_f444w'] > 0) &
        (merged_catalog['e_f160w'] != -99) &
        (merged_catalog['e_f444w'] != -99) &
        (merged_catalog['flag_1'] < 400) &
        (merged_catalog['flag_2'] < 400) &
        (merged_catalog['flag_2'] >= 0)
    )

    filtered_catalog = merged_catalog[valid_mask]

    F160W_flux = filtered_catalog['f_f160w']
    F444W_flux = filtered_catalog['f_f444w']

    F160W_mag = -2.5 * np.log10(F160W_flux) + 23.9
    F444W_mag = -2.5 * np.log10(F444W_flux) + 23.9

    color = F160W_mag - F444W_mag

    red_mask = (color > 2.0)
    blue_mask = (color <= 2.0)

    return filtered_catalog, red_mask, blue_mask, F160W_mag, F444W_mag, color
