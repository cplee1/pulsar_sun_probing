#!/bin/env python

import os
import csv
import argparse
import numpy as np
from psrqpy import QueryATNF
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(usage='Usage: %(prog)s [options]',
        description='Find pulsars within a given distance of the ecliptic.')
    
    parser.add_argument('-d', '--distance',
                  action='store', type=float, dest='dist', default=15.0,
                  help='Maximum distance from ecliptic in degrees [default: %(default)s]')
    parser.add_argument('--lookup',
                  action='store', type=str, dest='fn_lookup', default='ecliptic_coordinate_lookup.csv',
                  help='File name of lookup table (will generate on first execution) [default: %(default)s]')
    parser.add_argument('-o', '--outfile',
                  action='store', type=str, dest='fn_outfile', default='pulsars_near_ecliptic.csv',
                  help='Name of the file to write the list of pulsars to [default: %(default)s]')
    
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


def get_ecliptic_coordinates_from_ATNF(filename):
    """Create a lookup table of ecliptic coordinates for pulsars in the ATNF pulsar catalogue.

    Parameters
    ----------
    filename : string
        The name of file to save the lookup table.
    """
    # Query the ATNF pulsar catalogue
    query = QueryATNF(params=['PSRJ', 'RAJD', 'DECJD'], checkupdate=True)

    # Loop through catalogue to find pulsars near ecliptic
    pulsars_lb = []
    for pulsar in tqdm(np.array(query.table['PSRJ'])):
        ra = np.array(query[pulsar]['RAJD'])[0]
        dec = np.array(query[pulsar]['DECJD'])[0]
        l, b = equatorial2ecliptic(ra, dec)
        pulsars_lb.append([pulsar, l, b])

    # Write list of pulsars to file
    with open(filename, 'w') as f:
        spamwriter = csv.writer(f, delimiter=' ')
        spamwriter.writerow(['PSRJ', 'lambda', 'beta'])
        for pulsar in pulsars_lb:
            spamwriter.writerow([pulsar[0], round(pulsar[1],4), round(pulsar[2],4)])


def main():
    args = parse_args()

    # If a lookup table is not found, create one
    if not os.path.isfile(args.fn_lookup):
        get_ecliptic_coordinates_from_ATNF(args.fn_lookup)

    # Find pulsars within specified ecliptic lattitude range
    with open(args.fn_lookup, 'r') as f:
        spamreader = csv.reader(f, delimiter=' ')
        # Skip the header
        next(spamreader)
        near_ecliptic = []
        for row in spamreader:
            pulsar, l, b = row
            if abs(float(b)) < abs(args.dist):
                near_ecliptic.append([pulsar, abs(float(b))])

    # Write list of pulsars near ecliptic to file
    with open(args.fn_outfile, 'w') as f:
        spamwriter = csv.writer(f, delimiter=' ')
        spamwriter.writerow(['PSRJ', '|beta|'])
        for pulsar in near_ecliptic:
            spamwriter.writerow([pulsar[0], round(pulsar[1],4)])


if __name__ == '__main__':
    main()