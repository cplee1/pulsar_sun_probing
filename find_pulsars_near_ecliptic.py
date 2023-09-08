#!/bin/env python

import argparse
import numpy as np
from psrqpy import QueryATNF


def parse_args():
    parser = argparse.ArgumentParser(usage='Usage: %(prog)s [options]',
        description='Find pulsars within a given distance of the ecliptic.')
    
    parser.add_argument('-d', '--distance',
                  action='store', type=float, dest='dist', default=15.0,
                  help='Maximum distance from ecliptic in degrees [default: %(default)s]')
    parser.add_argument('-F', '--filename',
                  action='store', type=str, dest='filename', default='ecliptic_pulsars.csv',
                  help='Name of the file to write to [default: %(default)s]')
    
    args = parser.parse_args()

    if args.dist < 0 or args.dist > 90:
        parser.error('invalid distance given')

    return args


def equatorial2ecliptic(ra, dec):
    """Convert equatorial coordinates to ecliptic coordinates.

    Parameters
    ----------
    ra : float
        Right ascension in degrees.
    dec : float
        Declination in degrees.

    Returns
    -------
    lambda : float
        Ecliptic longitude in degrees.
    beta : float
        Ecliptic latitude in degrees.
    """
    ra = np.radians(ra)
    dec = np.radians(dec)
    epsilon = np.radians(23.439281)
    l = np.arctan2(np.sin(ra)*np.cos(epsilon) + np.tan(dec)*np.sin(epsilon), np.cos(ra))
    b = np.arcsin(np.sin(dec)*np.cos(epsilon) - np.cos(dec)*np.sin(epsilon)*np.sin(ra))
    return np.degrees(l), np.degrees(b)


def main():
    args = parse_args()

    # Query the ATNF pulsar catalogue
    query = QueryATNF(params=['PSRJ', 'RAJD', 'DECJD'], checkupdate=True)
    
    # Loop through catalogue to find pulsars near ecliptic
    ecliptic_pulsars = []
    for pulsar in np.array(query.table['PSRJ']):
        ra = np.array(query[pulsar]['RAJD'])[0]
        dec = np.array(query[pulsar]['DECJD'])[0]
        l, b = equatorial2ecliptic(ra, dec)
        print(pulsar, l, b)
        if abs(b) < abs(args.dist):
            ecliptic_pulsars.append([pulsar, abs(b)])

    # Write list of pulsars to file
    with open(args.filename, 'w') as f:
        f.write('# PSRJ |beta|\n')
        for pulsar in ecliptic_pulsars:
            f.write(f'{pulsar[0]} {round(pulsar[1],4)}\n')


if __name__ == '__main__':
    main()