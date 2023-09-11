#!/usr/bin/env python

import textwrap
import argparse
import numpy as np
from astropy.table import QTable
import astropy.units as u
from psrqpy import QueryATNF
from pulsar_spectra.spectral_fit import find_best_spectral_fit, estimate_flux_density
from pulsar_spectra.catalogue import collect_catalogue_fluxes


def parse_args():
    parser = argparse.ArgumentParser(usage='Usage: %(prog)s [options]',
        description='Find the brightest pulsars in a given observation.')
    
    parser.add_argument('-q', '--quiet',
                  action='store_false', dest='verbose', default=True,
                  help="Do not print status messages to stdout")
    parser.add_argument('-d', '--dmlim',
                  action='store', type=float, dest='dmlim', default=None,
                  help='DM limit in pc/cc [default: %(default)s]')
    parser.add_argument('-f', '--obsfreq',
                  action='store', type=float, dest='obsfreq', default=154.24,
                  help='Frequency of observation in MHz [default: %(default)s MHz]')
    parser.add_argument('-l', '--lowfreqs',
                  action='store_true', dest='lowfreqs', default=False,
                  help='Must have low frequency (< 600 MHz) flux data available? [default: %(default)s]')
    parser.add_argument('-m', '--minimum_flux',
                  action='store', type=float, dest='min_flux', default=100.0,
                  help='Minimum desired flux density in mJy [default: %(default)s mJy]')
    parser.add_argument('-p', '--pulsars',
                  action='store', type=str, dest='fn_psrs', default='pulsars_near_ecliptic.csv',
                  help='Name of the file with a list of pulsars [default: %(default)s]')
    parser.add_argument('-o', '--outfile',
                  action='store', type=str, dest='fn_outfile', default='bright_pulsars_near_ecliptic.mrt',
                  help='Name of the file to write the output data to [default: %(default)s]')
    parser.add_argument('--plot_best',
                  action='store_true', dest='plot_best', default=False,
                  help="Plot spectra for modelled pulsars, showing the best-fitting model.")
    
    args = parser.parse_args()

    if type(args.dmlim) is not type(None):
        if args.dmlim <= 0:
            parser.error('invalid DM given')

    if args.obsfreq <= 0:
        parser.error('invalid frequency given')

    if args.min_flux <= 0:
        parser.error('invalid minimum flux density given')

    return args


def get_flux(cat_dict, pulsar, interpfreq, lowfreqs=False, plot_best=False):
    """Estimate the flux density of a pulsar based on the best-fit spectral model.

    Parameters
    ----------
    cat_dict : `dict`
        pulsar_spectra catalogue dictionary.
    pulsar : `string`
        J name of pulsar to model.
    interpfreq : `float`
        Frequency to estimate the flux density at in MHz.
    lowfreqs : `boolean`, optional
        Only consider pulsars with flux density measurements available < 600 MHz. |br| Default: False.

    Returns
    -------
    fitted_flux : `float`
        Estimate of the flux density in mJy.
    fitted_flux_err : `float`
        Uncertainty in the estimated flux density in mJy.
    """
    if pulsar not in cat_dict:
        print(f'{pulsar} not in catalogue')
        return None, None
    
    freqs, bands, fluxs, flux_errs, refs = cat_dict[pulsar]
    
    if len(freqs) < 3:
        print(f'Not enough flux measurements for {pulsar}')
        return None, None
    
    if lowfreqs:
        if np.min(freqs) > 600:
            print(f'No low frequency flux measurements for {pulsar}')
            return None, None
    
    model, m, _, _, _ = find_best_spectral_fit(pulsar, freqs, bands, fluxs, flux_errs, refs, plot_best=plot_best)
    
    if type(model) is type(None):
        print(f'No best fit model for {pulsar}')
        return None, None
    
    fitted_flux, fitted_flux_err = estimate_flux_density(interpfreq, model, m)
    
    return fitted_flux, fitted_flux_err


def print_results_table(query, flux_matrix, args):
    """Print pulsars and flux densities using an astropy QTable.

    Parameters
    ----------
    query : `pandas.DataFrame`
        A psrqpy query object.
    flux_matrix : `list`
        A list of tuples in the form (Jname, fitted_flux, fitted_flux_err).
    """
    jnames = []
    p0s = []
    dms = []
    fluxs = []
    flux_errs = []
    rajs = []
    decjs = []
    for flux_tup in flux_matrix:
        if type(flux_tup[1]) is not type(None) and type(flux_tup[2]) is not type(None):
            jnames.append(flux_tup[0])
            p0s.append(query.get_pulsar(flux_tup[0])['P0'][0]*1e3 * u.ms)
            dms.append(query.get_pulsar(flux_tup[0])['DM'][0] * u.pc / u.cm**3)
            fluxs.append(flux_tup[1] * u.mJy)
            flux_errs.append(flux_tup[2] * u.mJy)
            rajs.append(query.get_pulsar(flux_tup[0])['RAJ'][0])
            decjs.append(query.get_pulsar(flux_tup[0])['DECJ'][0])

    tab = QTable([np.arange(len(jnames))+1, jnames, p0s, dms, fluxs, flux_errs, rajs, decjs],
                 names=('Rank', 'JNAME', 'P0', 'DM', 'Flux', 'Flux_err', 'RAJ', 'DECJ'),
                 descriptions=(
                     'Rank',
                     'Pulsar name (J2000)',
                     'Pulsar spin period',
                     'Dispersion measure',
                     'Estimated flux density',
                     'Uncertainty in estimated flux density',
                     'Right ascension (J2000)',
                     'Declination (J2000)'
                 ))
    tab['P0'].info.format = '.3f'
    tab['DM'].info.format = '.3f'
    tab['Flux'].info.format = '.3f'
    tab['Flux_err'].info.format = '.3f'

    tab.write(args.fn_outfile, format='ascii.mrt', overwrite=True)

    if type(args.dmlim) is type(None):
        if args.lowfreqs:
            comment_title = f'Title: Pulsars with an estimated flux density of > {args.min_flux:.2f} mJy ' + \
            f'and flux density data available below 600 MHz.'
        else:
            comment_title = f'Title: Pulsars with an estimated flux density of > {args.min_flux:.2f} mJy.'
    else:
        if args.lowfreqs:
            comment_title = f'Title: Pulsars with an estimated flux density of > {args.min_flux:.2f} mJy, ' + \
            f'a DM > {args.dmlim:.2f} pc/cc, and flux density data available below 600 MHz.'
        else:
            comment_title = f'Title: Pulsars with an estimated flux density of > {args.min_flux:.2f} mJy ' + \
            f'and a DM > {args.dmlim:.2f} pc/cc.'

    wrapped_comment_title = textwrap.fill(comment_title, width=80)

    with open(args.fn_outfile, 'r') as f:
        lines = f.readlines()
    
    lines = lines[1:]
    lines.insert(0, wrapped_comment_title + '\n')
    with open(args.fn_outfile, 'w') as file:
        file.writelines(lines)


def main():
    args = parse_args()

    if type(args.dmlim) is type(None):
        cond = None
    else:
        cond = f'DM < {args.dmlim}'

    pulsars = list(np.loadtxt(args.fn_psrs, dtype=str, unpack=True, usecols=0))
    query = QueryATNF(psrs=pulsars, params=['JName', 'RAJD', 'DECJD'], condition=cond, checkupdate=True)
    print(f'Catalogue version {query.get_version}')

    flux_matrix = []
    cat_dict = collect_catalogue_fluxes()
    for pulsar in query.table['JNAME']:
        print(f'Estimating flux density for {pulsar}')
        fitted_flux, fitted_flux_err = get_flux(cat_dict, pulsar, args.obsfreq, args.lowfreqs, args.plot_best)
        if type(fitted_flux) is type(None):
            continue
        if fitted_flux > 30000:
            # Anything brighter than this probably wasn't modelled properly
            continue
        if fitted_flux > args.min_flux:
            flux_matrix.append((pulsar, fitted_flux, fitted_flux_err))

    # Sort pulsars by flux density
    flux_matrix.sort(key=lambda x: x[1], reverse=True)

    print_results_table(query, flux_matrix, args)


if __name__ == '__main__':
    main()
