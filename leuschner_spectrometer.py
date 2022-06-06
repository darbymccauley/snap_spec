###########################################################################
## This module provides tools for interacting  with the spectrometer at 
## the UC Berkeley Leuschner Radio Observatory. 

## Darby McCauley

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
        #self.fpgfile =
        #self.nchan = 
        #self.downsample = 
        #self.fft_shift = 
        #self.acc_len = 

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
        unix_object = Time(obs_start_unix, format='unix', location=(LON, LAT, ALT)) #unix time Time object
        obs_start_jd = unix_object.jd #convert unix time to julian date

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

    def make_fits_cols(name, data):
    """
    Create a FITS column of double-precision floating data.

    Inputs:
    - name: Name of the FITS column.
    - data: Array of data for the FITS column.
    """
        return fits.Column(name=name, 'D', array=data)


    # CONFUSED ABOUT THIS FUNCTION -- RACHEL'S init_spec()
    def initialize_spec(self, scale=False):
    """
    Starts the fpg process on the SNAP and initializes the 
    spectrometer.

    Inputs:
    - scale: Whether or not to scale down each integration by the 
    total number of spectra per integration time.
    - force_restart: Restart the fpg process even if it is already
    running.
    """
        print('Starting the spectrometer...')
        
        self.fpga.upload_to_ram_and_program(self.fpgfile)
        self.fpga.write_int('fft_shift', self.fft_shift)

        self.fpga.write_int('corr_0_acc_len') = _____
        self.fpga.write_int('corr_1_acc_len') = _____
        
        # Sync pulse lets the spectrometer know when to start.
        for i in (0,1,0):
            self.fpga.write_int('arm', i)

        self.count = self.fpga.read_int('corr_0_acc_cnt')
        self.count = self.fpga.read_int('corr_1_acc_cnt')

        print('Spectrometer is ready.')


    # PROBABLY NEEDS A LOT OF WORK
    def poll(self):
    """
    Waits until the integration count has been incrimented and
    returns the date of the integration in seconds (s) since 
    Jan 1, 1970 UTC.

    Returns:
    - obs_date: Unix time of the integration.
    """
        self.count = self.fpga.read_int('acc_num')
        while self.fpga.read_int('acc_num') == self.count:
            time.sleep(0.1)
        obs_date = time.time() - 0.5*self.int_time
        self.count = self.fpga.read_int('acc_num')
        return obs_date

        
    # DONT KNOW WHAT TO DO WITH THIS ATM
    #def read_bram(self, bram):


    def read_spec(self, filename, nspec, coords, coord_sys='ga', bandwidth=12e6):
    """
    Recieves data from the Leuschner spectrometer and saves it to a
    FITS file. The primary HDU contains information about the
    observation (coordinates, number of spectra collected, time,
    etc.) and spectrometer attributes used. Each set of spectra is
    stored in its own FITS table in the FITS file. The columns in
    each FITS table are ''auto0_real'', ''auto1_real'',
    ''cross_real'', and ''cross_imag''. All columns contained
    double-precision floating-point numbers.

    Inputs:
    - filename: Name of the output FITs file.
    - nspec: Number of spectra to collect.
    - coords: Coordinate(s) of the target.
        Format: (l/ra, b/dec)
    - coord_sys: Coordinate system used for ''coords''.
        Default is galactic coordinates. Takes in either galactic 
        ('ga') or equatorial ('eq') coordinate systems.
    Returns:
    - FITS file with collected spectrometer data.
    """
        # Ensure that a proper coordinate system has been supplied
        if coord_sys != 'ga' and coord_sys != 'eq':
            raise ValueError('Invalid coordinate system supplied: ' + coord_sys)

        # Make sure the spectrometer is actually running
        if self.fpga.is_running() == True:
            continue
        elif self.fpga.is_running() == False:
            self.initialize_spec()
        self.spec_props(bandwidth)

        # BRAM device names
        bram_names = ['auto0_real', 'auto1_real', 'cross_real', 'cross_imag']
        bram_devices = map(lambda name: 'spec_' + name, bram_names)

        # Create a primary FITS HDU for a table file
        hdulist = [self.fits_header(nspec, coords, coord_sys)]

        # Read some number of spectra to a FITS file
        print('Reading', nspec, 'spectra from the SNAP.')
        ninteg = 0
        while ninteg < nspec:
            # Update the counter and read the spectra from the SNAP
            spec_date = self.poll()
            try:
                spectra = map(self.read_bram(), bram_devices)
            except RuntimeError:
                print('WARNING: Cannot reach the SNAP. Skipping integration.')
                self.fpga.connect()
                continue

            # Create FITS columns with the data
            fcols = map(make_fits_cols, bram_names, spectra)
            hdulist.append(fits.BinTableHDU.from_columns(fcols))

            # Add the accumulation date in several formats to the header
            unix_spec = Time(spec_date, format='unix', location=(LON, LAT, ALT))
            julian_date_spec = unix_spec.jd
            utc_spec = time.asctime(time.gmtime(spec_date))
            hdulist[-1].header['JD'] = (julian_date_spec, 'Julian date of observation.')
            hdulist[-1].header['UTC'] = (utc_spec, 'UTC time of observation.')
            hdulist[-1].header['UNIX TIME'] = (spec_date, 'Seconds since epoch.')

            ninteg += 1
            integ_time = ninteg*self.int_time
            print('Integration count:', ninteg, '(' + str(integ_time), 's)')

        # Save the output file
        print('Saving spectra to output file:', filename)
        fits.HDUList(hdulist).writeto(filename, clobber=True)


    # DON'T THINK MODE REGISTER EXISTS SO IS NEEDED??
    def reconnect(self):
    """
    Runs if the spectrometer can't be reached in the middle of data
    collection. This should only be run if the fpg process for the
    spectrometer has already started.
    """
        while True:
            try:
                self.fpga.read_int('mode')
                break
            except:
                time.sleep(0.1)


    def set_fft_shift(self, fft_shift):
    """
    Allows the user to change the FFT shifting instructions on the
    SNAP.
    
    Inputs:
    - fft_shift: FFT shifting instructions for the SNAP.
    """
    self.fft_shift = int(fft_shift)
    self.fpga.write_int('fft_shift', self.fft_shift)
    for i in range(2):
        self.poll()


    def set_scale(self, scale):
    """
    Scales the spectra as desired by the number of spectra 
    integrated per accumulation.

    Inputs:
    - scale: Whether or not to downscale the spectra.
    """

    self.scale = int(scale)
    self.fpga.write_int('scale', self.scale)
    for i in range(2):
        self.poll()

    #def spec_props(self, bandwidth):
    #"""
    #Takes in a sample rate and stores the parameters computed from
    #it (bandwidth, resolution, etc...).
    #
    #Inputs:
    #"""








