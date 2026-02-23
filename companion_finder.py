import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18 as cosmo
from astropy.table import Table
from scipy.optimize import brentq

def comoving_kpc(z):
    return cosmo.comoving_distance(z).value * 1000.0

def find_companions(filtered_catalog, source_mask, z_all, z_edges, F160W_mag, F444W_mag, color, type_label):

    coords_all = SkyCoord(
        ra=filtered_catalog['RA_1'],
        dec=filtered_catalog['DEC_1'],
        unit='deg'
    )

    sources = filtered_catalog[source_mask]
    coords_sources = SkyCoord(
        ra=sources['RA_1'],
        dec=sources['DEC_1'],
        unit='deg'
    )

    r_min_mpc = 5.0 / 1000.0
    r_max_mpc = 50.0 / 1000.0

    with_companions = Table(
        names=('ID', 'zphot', 'n_companions'),
        dtype=('i8', 'f8', 'i8')
    )

    for i, source in enumerate(sources):

        theta = coords_sources[i].separation(coords_all).to(u.rad).value
        d_ang = cosmo.angular_diameter_distance(source['zphot']).to('Mpc').value

        physical_sep = d_ang * theta

        neighbors = np.where(
            (physical_sep >= r_min_mpc) &
            (physical_sep <= r_max_mpc) &
            (filtered_catalog['ID'] != source['ID'])
        )[0]

        if len(neighbors) > 0:
            with_companions.add_row([
                source['ID'],
                source['zphot'],
                len(neighbors)
            ])

    return with_companions
