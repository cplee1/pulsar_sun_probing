# pulsar_sun_probing
Code to assist with planning pulsar observations near the Sun.

## Scripts
This repository includes three Python scripts in `/scripts` which are designed to be
used in succession. The scripts require the following Python packages:

- [NumPy](https://numpy.org/)

- [Astropy](https://www.astropy.org/)

- [psrqpy](https://psrqpy.readthedocs.io/en/latest/) - for querying the ATNF pulsar catalogue

- [pulsar_spectra](https://pulsar-spectra.readthedocs.io/en/latest/) - for spectral modelling

- [tqdm](https://tqdm.github.io/) - for the progress bar

### find_pulsars_near_ecliptic.py
Compiles a list of all known pulsars lower than a given ecliptic latitude (default 15 deg).
Since converting from equitorial to ecliptic coordinates can take a while, the script will
store these coordinates in `data/ecliptic_coordinate_lookup.csv` so that they are not
re-computed each time the script runs. The list of pulsars near the ecliptic is then stored
in `data/pulsars_near_ecliptic.csv` which also includes the distance of each pulsar to
the ecliptic (|beta|) in degrees.

### find_brightest_pulsars.py
Finds the best-fit spectral model for each pulsar and interpolates the models to a
given frequency to estimate the flux density. Ranks pulsars based on estimated flux
densities and writes the results to `data/bright_pulsars_near_ecliptic.mrt` as
a [Machine Readable Table (MRT)](https://journals.aas.org/mrt-standards/).

To improve the reliability of the spectral fits, I have also included the `-l` option
to only include pulsars with flux density measurements available below 600 MHz.
Furthermore, the `-d DM` option can be used to restrict the list to pulsars below a
given DM limit (in pc/cc). The example data was created with the following command:

    find_brightest_pulsars.py -l -d 100


### plot_equatorial_skymap.py
Plots a full-sky map in equatorial (RA and DEC) coordinates, indicating the location
of the ecliptic, the latitude limit, and the pulsars selected by `find_brightest_pulsars.py`.
Saves the plot to `data/sunpath_skymap.png`.

![equatorial skymap](https://github.com/cplee1/pulsar_sun_probing/blob/786ab262d3f7f266b72f18c5704c67cf5530d72d/data/sunpath_skymap.png)