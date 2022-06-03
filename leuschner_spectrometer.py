###########################################################################
## This module provides tools for interacting  with the spectrometer at 
## the UC Berkeley Leuschner Radio Observatory. 

## Darby McCauley

## ROUGH DRAFT
###########################################################################

import casperfpga
import numpy as np
import astropy
import time 
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from astropy.io import fits

# Leuschner Observatory coordinates
ALT, LAT, LON = 304.0, 37.9183, -122.1067

# Define the fpga
#fpga = casperfpga.CasperFpga('hostname')

# Create Spectrometer class
class Spectrometer(object):
    """
    Casperfpga interface to the SNAP spectrometer.
    """
    def __init__(self, hostname):
    """
    Create the interface to the SNAP.
    
    Inputs:
    - hostname: IP address of the fpga.
    """
        self.fpga = casperfpga.CasperFpga(hostname)
        
        self.hostname = hostname
        #self.mode = 
        #self.count =
        #self.scale =
        #self.boffile =
        #self.nchan = 
        #self.downsample = 
        #self.fft_shift = 
        # self.acc_len = 

    def check_if_connected():
    """
    Checks if the SNAP is connected and raises an IOError if the client
    cannot reach the SNAP.
    """
        
        if self.fpga.is_connected() == True:
            print('Connection to the SNAP established.')
        elif self.fpga.is_connected() == False:
            raise IOError('NOT connected to the SNAP.')

    # IS THIS FUNCTION NEEDED?
    #def check_if_running():
    """
    Checks to see if the fpga process for the spectrometer has been
    initialized on the SNAP.
    """
    

    def fits_header(self, nspec, coords, coord_sys='ga'):
    """
    Creats the primary HDU (header) of the data collection FITS file. 
    Writes in observation attributes such as time of observation, number of
    spectra collected, and the coordinates of the observation target.
    
    Inputs:
    - nspec: Number of spectra to collect.
    - coords: Coordinate(s) of the target.
        Format: (l/ra, b/dec)
    - coord_sys: Coordinate system used for ''coords''.
        Default is galactic coordinates. Takes in either galactic ('ga') or
        equatorial ('eq') coordinate systems.
    Returns:
    - FITS file primary HDU information containing the attributes of the 
    observation.
    """
        # Ensure that a proper coordinate system has been supplied
        if coord_sys != 'ga' and coord_sys != 'eq':
            raise ValueError('Invalid coordinate system supplied: ' + coord_sys)
        # Set times
        obs_start_unix = time.time() #unix time
        unix_class = Time(obs_start_unix, format='unix', location=(LON, LAT, ALT)) #unix time class
        obs_start_jd = unix_class.jd #convert unix time to julian date

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

    # CONFUSED ABOUT HOW TO WRITE THIS FUNCTION -- RACHEL'S init_spec()
    def initialize_spec(self, scale=False, force_restart=False):
    """
    Starts the bof process on the SNAP. First checks to see if the
    spectrometer is already running, then initializes it if it isn't.

    Inputs:
    - scale: Whether or not to scale down each integration by the total
    number of spectra per integration time.
    - force_restart: Restart the bof process even if it is already
    running.
    """
    
        # Check if the spectrometer is running
        if self.fpga.is_running() == True:
            print('Spectrometer is running.')
        elif self.fpga.is_running() == False:
            print('Starting the spectrometer...')

    
    

                








