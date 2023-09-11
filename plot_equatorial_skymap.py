#!/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
from astropy.io import ascii
from astropy.coordinates import SkyCoord
import astropy.units as u
from psrqpy import QueryATNF

plt.rcParams['mathtext.fontset'] = 'dejavuserif'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 12


def parse_args():
    parser = argparse.ArgumentParser(usage='Usage: %(prog)s [options]',
        description='Make an equatorial skymap showing pulsars along the ecliptic. ' +\
            'User-specified beta/flux limits are for plotting purposes only.')
    
    parser.add_argument('-b', '--beta_limit',
                  action='store', type=float, dest='beta_lim', default=15.0,
                  help='Absolute value of ecliptic latitude limit [default: %(default)s]')
    parser.add_argument('-S', '--flux_limit',
                  action='store', type=float, dest='flux_lim', default=100.0,
                  help='Flux density limit for bright pulsars [default: %(default)s mJy]')
    parser.add_argument('--psr_ecliptic_file',
                  action='store', type=str, dest='fn_ecl', default='pulsars_near_ecliptic.csv',
                  help='CSV file where the first column lists the pulsars along the ecliptic [default: %(default)s]')
    parser.add_argument('--psr_bright_file',
                  action='store', type=str, dest='fn_bright', default='bright_pulsars_near_ecliptic.mrt',
                  help='Machine readable table (MRT) file listing the bright pulsars along the ecliptic [default: %(default)s]')
    parser.add_argument('-o', '--outfile',
                  action='store', type=str, dest='fn_outfile', default='sunpath_skymap.png',
                  help='Name of the output plot file [default: %(default)s]')
     
    args = parser.parse_args()

    return args


def ecliptic2equatorial(l, b):
    """Convert ecliptic coordinates to equatorial coordinates (J2000).

    Parameters
    ----------
    lambda : float
        Ecliptic longitude in degrees.
    beta : float
        Ecliptic latitude in degrees.

    Returns
    -------
    ra : float
        Right ascension in degrees.
    dec : float
        Declination in degrees.
    """
    # Create a SkyCoord object with ecliptic coordinates
    ecliptic_coord = SkyCoord(lon=l * u.degree, lat=b * u.degree, frame='barycentrictrueecliptic')

    # Transform to equatorial coordinates
    equatorial_coord = ecliptic_coord.transform_to('icrs')

    # Extract right ascension (ra) and declination (dec)
    ra = equatorial_coord.ra.degree
    dec = equatorial_coord.dec.degree

    return ra, dec


def map_coordinates(l, b):
    """Map list of ecliptic coordinates to sorted list of equatorial coordinates.

    Parameters
    ----------
    lambda : float
        Ecliptic longitude in degrees.
    beta : float
        Ecliptic latitude in degrees.

    Returns
    -------
    ra : float
        Right ascension in degrees.
    dec : float
        Declination in degrees.
    """
    # Convert coordinates
    ra, dec = ecliptic2equatorial(l, b)

    # Add to numpy array
    coords = np.empty(shape=(len(dec),2))
    for i, (ra_i, dec_i) in enumerate(zip(ra, dec)):
        coords[i][0] = ra_i
        coords[i][1] = dec_i

    # Sort by right ascension
    sorted_indices = np.argsort(coords[:,0])
    sorted_coords = coords[sorted_indices]

    # Separate into two arrays
    ra = sorted_coords[:,0]
    dec = sorted_coords[:,1]

    return ra, dec


def plot_equatorial_skymap(query_all, query_ecliptic, query_ecliptic_bright, 
                           beta_lim=15.0, flux_lim=100.0, filename='sunpath_skymap.png'):
    """Plot an equatorial skymap showing the sunpath and nearby pulsars.

    Parameters
    ----------
    query_all : `pandas.DataFrame`
        A psrqpy query for all pulsars.
    query_ecliptic : `pandas.DataFrame`
        A psrqpy query for pulsars within beta_lim of the ecliptic.
    query_ecliptic_bright : `pandas.DataFrame`
        A psrqpy query for pulsars within beta_lim of the ecliptic and with
        flux densities less than flux_lim.
    beta_lim : `float`, optional
        Ecliptic latitude limit in degrees. |br| Default: 15.0.
    flux_lim : `float`, optional
        Flux density limit in mJy. |br| Default: 100.0.
    filename : `string`, optional
        Name of the figure plot. |br| Default: 'sunpath_skymap.png'.
    """
    ras_all = np.array(query_all.table['RAJD'])
    decs_all = np.array(query_all.table['DECJD'])
    ras_ecl = np.array(query_ecliptic.table['RAJD'])
    decs_ecl = np.array(query_ecliptic.table['DECJD'])
    ras_ecl_bright = np.array(query_ecliptic_bright.table['RAJD'])
    decs_ecl_bright = np.array(query_ecliptic_bright.table['DECJD'])

    lambda_patch = np.linspace(0,360,500)
    beta_patch_centre = np.full_like(lambda_patch, 0.)
    beta_patch_top = np.full_like(lambda_patch, abs(beta_lim))
    beta_patch_bottom = np.full_like(lambda_patch, -abs(beta_lim))

    centre_coords = map_coordinates(lambda_patch, beta_patch_centre)
    top_coords = map_coordinates(lambda_patch, beta_patch_top)
    bottom_coords = map_coordinates(lambda_patch, beta_patch_bottom)

    _, ax = plt.subplots(figsize=(8, 6))
    psr1 = ax.scatter(ras_all, decs_all, s=4, linewidths=0, marker='.', color='grey', zorder=2)
    psr2 = ax.scatter(ras_ecl, decs_ecl, s=7, linewidths=0, marker='.', color='k', zorder=2)
    psr3 = ax.scatter(ras_ecl_bright, decs_ecl_bright, s=20, marker='*', color='crimson', zorder=3)
    ecl1, = ax.plot(*centre_coords, color='grey', linestyle='-', linewidth=1, zorder=1)
    ax.plot(*top_coords, color='grey', linestyle='--', linewidth=0.8, zorder=1)
    ecl2, = ax.plot(*bottom_coords, color='grey', linestyle='--', linewidth=0.8, zorder=1)

    ax.set_xlabel('Right ascension (deg)')
    ax.set_ylabel('Declination (deg)')
    ax.set_xlim(0, 360)
    ax.set_ylim(-90,90)
    ax.set_yticks(np.arange(-90, 91, 30))
    ax.set_xticks(np.arange(0, 361, 60))
    ax.grid(ls=':', alpha=0.5)
    ax.legend(
        [
            (psr1, psr2),
            psr3,
            ecl1,
            ecl2
        ],
        [
            'All known pulsars',
            f'Pulsars with $S_{{150}}>{flux_lim:.0f}\,\mathrm{{mJy}}$',
            r'$\beta=0\,\mathrm{deg}$',
            f'$|\\beta|<{beta_lim:.0f}\,\mathrm{{deg}}$'
        ],
        numpoints=1,
        handler_map={tuple: HandlerTuple(ndivide=None)},
        loc='upper right'
    )

    plt.savefig(filename, dpi=300, bbox_inches='tight')


def main():
    args = parse_args()

    # Read files
    ecliptic_pulsars = list(np.loadtxt(args.fn_ecl, dtype=str, unpack=True, usecols=0, skiprows=1))
    bright_ecliptic_pulsars_tab = ascii.read(args.fn_bright, format='mrt')

    # Convert astropy table column to list
    bright_ecliptic_pulsars = list(bright_ecliptic_pulsars_tab['JNAME'])

    # Query the ATNF catalogue
    query_all = QueryATNF(params=['PSRJ', 'RAJD', 'DECJD'], checkupdate=True)
    query_ecliptic = QueryATNF(params=['PSRJ', 'RAJD', 'DECJD'], psrs=ecliptic_pulsars)
    query_ecliptic_bright = QueryATNF(params=['PSRJ', 'RAJD', 'DECJD'], psrs=bright_ecliptic_pulsars)
    
    plot_equatorial_skymap(query_all, query_ecliptic, query_ecliptic_bright, args.beta_lim, args.flux_lim, args.fn_outfile)


if __name__ == '__main__':
    main()