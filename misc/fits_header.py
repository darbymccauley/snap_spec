#################################################################
## This module places the collected data from the UC Berkeley 
## Leuschner Radio Observatory spectrometer into a fits file.

## ROUGH DRAFT

#################################################################

import astropy
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import time
from astropy.time import Time

# Leuschner Observatory coordinates
ALT, LAT, LON = 304.0, 37.9183, -122.1067

# Develop FITs file header as primary HDU
def fits_header(nspec, coords, coord_sys='ga'):
    """
    Creates the primary HDU (header) of the data collection FITS
    file. Writes in observation attributes such as time of 
    observation, number of spectra collected, and the coordinates
    of the observation target.

    Inputs:
    - nspec: Number of spectra to collect.
    - coords: Coordinate of the target.
        Format: (l/ra, b/dec)
    - coord_sys: Coordinate system used for ''coords''.
        Default is galactic coordinates. Can use galactic ('ga') 
        or equatorial ('eq') coordinate systems.
    Returns:
    - FITS file primary HDU information containing the attributes
    of the observation.
    """
    # Ensure that a proper coordinate system has been supplied
    if coord_sys != 'ga' and coord_sys != 'eq':
        raise ValueError('Invalid coordinate system supplied: ' + coord_sys)


    # Set times
    obs_start_unix_time = time.time() # unix time
    unix_time = Time(obs_start_unix_time, format='unix', location=(LON, LAT, ALT)) # unix time class
    obs_start_jd = unix_time.jd # convert unix time to julian date

    # Set the coordinates
    if coord_sys == 'ga':
        l, b = coords*u.degree
        c = SkyCoord(l=l, b=b, frame='galactic')
        equatorial = c.fk5 
        ra, dec = equatorial.ra, equatorial.dec
    elif coord_sys == 'eq':
        ra, dec = coords*u.degree
        c = SkyCoord(ra, dec)
        galactic = c.galactic
        l, b = galactic.l, galactic.b


    # Create header and write observation attributes into it
    header = fits.Header()
    # FPGA = 250MHz
    # ADC = 500MHz 
    # BW = 250MHz
    # N channels = 8192
    # Res = 250e6/8192 = 30.517kHz
    # FFT shift (rough)
    # acclen =? inttime (rough)
    header['NSPEC'] = (nspec, 'Number of spectra collected')
    header['L'] = (l.value, 'Galactic longitude [deg]')
    header['B'] = (b.value, 'Galactic latitude [deg]')
    header['RA'] = (ra.value, 'Right Ascension [deg]')
    header['DEC'] = (dec.value, 'Declination [deg]')
    header['JD'] = (obs_start_jd, 'Julian date of start time')
    header['UNIX'] = (obs_start_unix_time, 'Seconds since epoch')

    return fits.PrimaryHDU(header=header)
    
