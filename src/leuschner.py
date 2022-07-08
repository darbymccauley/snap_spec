##########################################################################
# This script allows the user to interact with the spectrometer at the
# the UC Berkeley Leuschner Radio Observatory. The spectrometer is
# configured to a SNAP board, which is controlled by a rPi4. The two 
# communicate via casperfpga through the tcpborphserver3 server,
# reconfigured for SNAPs.
# 
# Darby McCauley 2022
# darbymccauley@berkeley.edu
###########################################################################


import casperfpga
from hera_corr_f import SnapFengine
import ugradio

import numpy as np
import time 
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from astropy.io import fits


# Create Spectrometer class
class Spectrometer(object):
    """
    Casperfpga interface to the SNAP spectrometer.
    """

    def __init__(self):
        """
        Create the interface to the SNAP.
        """
        self.host = 'localhost'
        self.fpgfile = 'fpga/ugradio_corrspec_2022-02-22_0905.fpg'
        
        # self.scale = 0
        # # self.adc_rate = 500e6
        # self.downsample = 1<<3
        # self.bandwidth = 250e6
        # self.samp_rate = self.bandwidth*2
        # self.nchan = 1<<13
        # self.resolution = self.bandwidth/self.nchan
        # self.fft_shift = 1<<14
        # self.acc_len = 1<<27
        # self.clock_rate = self.downsample*self.samp_rate # or 10 MHz?
        # self.int_time = self.acc_len/self.clock_rate


        self.fpga = casperfpga.CasperFpga(self.host)
        self.adc = casperfpga.snapadc.SNAPADC(self.fpga) # 10 MHz reference signal default       
        self.s = SnapFengine(self.host, transport='default')
        

    def check_connection(self):
        """
        Checks if the SNAP is connected and raises an IOError if the 
        client cannot reach the SNAP.
        """
        if self.fpga.is_connected():
            print("Connection to the SNAP established.")
        elif not self.fpga.is_connected():
            raise IOError("NOT connected to the SNAP.")

    
    def check_running(self):
        """
        Checks if the fpga has been programmed and is running.
        Returns an IOError is the SNAP is not running and cannot
        program.
        """
        if self.fpga.is_running() and self.s.is_programmed():
            print("Fpga is programmed and running.")
        elif not self.fpga.is_running() or not self.s.is_programmed():
            print("WARNING: Fpga is not running. Programming...")
            self.fpga.upload_to_ram_and_program(self.fpgfile)
            self.s.fpga.upload_to_ram_and_program(self.fpgfile)
            if self.fpga.is_running() and self.s.is_programmed():
                print("Fpga is now programmed and running.")
            else:
                raise IOError("Cannot program fpga.")


    def fits_header(self, nspec, coords, coord_sys='ga'):
        """
        Creates the primary HDU (header) of the data collection FITS 
        file. Writes in observation attributes such as time of 
        observation, number of spectra collected, and the coordinates 
        of the observation target.
        
        Inputs:
        - nspec: Number of spectra to collect.
        - coords: Coordinate(s) of the target.
            Format: (l/ra, b/dec)
        - coord_sys: Coordinate system used for ''coords''.
            Default is galactic coordinates. Takes in either galactic 
            ('ga') or equatorial ('eq') coordinate systems.
        Returns:
        - FITS file primary HDU information containing the attributes 
        of the observation and spectrometer.
        """
        # Ensure that a proper coordinate system has been supplied
        if coord_sys != 'ga' and coord_sys != 'eq':
            raise ValueError("Invalid coordinate system supplied: " + coord_sys)
        # Set times
        obs_start_unix = time.time() #unix time
        unix_object = Time(obs_start_unix, format='unix', location=(ugradio.leo.lon, ugradio.leo.lat, ugradio.leo.alt)) #unix time Time object
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

        # Create header and write spectrometer and observation attributes into it
        header = fits.Header()

        header['NSPEC'] = (nspec, "Number of spectra collected")
        header['FPGFILE'] = (self.fpgfile, "FPGA FPG file")
        header['HOST'] = (self.s.fpga.host, "Host of the FPGA")
        # header['MODE'] = (self.mode, "Spectrometer mode")
        header['CLK'] = (self.s.fpga.estimate_fpga_clock(), "FPGA clock speed [MHz]")
        # header['ADC'] = (self.adc_rate, "ADC clock speed [Hz]")
        header['ADC_NAME'] = (self.s.adc.adc.name, "Name of ADC")
        # header['DOWNSAMPLE'] = (self.downsample, "ADC downsampling period")
        # header['SAMPRATE'] = (self.samp_rate, "Downsampled clock speed [Hz]")
        # header['BW'] = (self.bandwidth, "Bandwidth of spectra [Hz]")
        # header['NCHAN'] = (self.nchan, "Number of frequency channels")
        # header['RES'] = (self.resolution, "Frequency resolution [Hz]")
        # header['FFTSHIFT'] = (self.fft_shift, "FFT shifting instructions")
        # header['ACCLEN'] = (self.acc_len, "Number of clock cycles")
        # header['INTTIME'] = (self.int_time, "Integration time of spectra")
        # header['SCALE'] = (self.scale, "Average instead of sum on SNAP")
        # header['PYTHON'] = (PYTHON_VERSION, "Python version")
        # header['SRC'] = (SRC_CODE, "Source code")
        # header['CASPERFPGA'] = (CASPERFPGA_VERSION, "casperfpga code used")
        # header['HERA_CORR_F'] = (HERA_CORR_F_VERSION, "hera_corr_f code used")

        header['L'] = (l.value, "Galactic longitude [deg]")
        header['B'] = (b.value, "Galactic latitude [deg]")
        header['RA'] = (ra.value, "Right Ascension [deg]")
        header['DEC'] = (dec.value, "Declination [deg]")
        header['JD'] = (obs_start_jd, "Julian date of start time")
        header['UNIX'] = (obs_start_unix, "Seconds since epoch")

        return fits.PrimaryHDU(header=header)


    def make_fits_cols(self, name, data):
        """
        Create a FITS column of double-precision floating data.

        Inputs:
        - name: Name of the FITS column.
        - data: Array of data for the FITS column.
        """
        return fits.Column(name=name, format='D', array=data)


    def initialize(self):
        """
        Programs the fpga on the SNAP and initializes the spectrometer.
        """
        print('Starting the spectrometer...')
        
        # Program fpga
        self.check_connection()
        self.check_running()
        
        # Initialize and align ADCs
        print("Aligning and initializing ADCs...")
        try:
            self.s.adc.init()
            self.s.align_adc()
            print("ADCs aligned and initialized.")
        except:
            try: # try again (usually works after two attempts)
                self.s.adc.init()
                self.s.align_adc()
                print("ADCs aligned and initialized.")
            except:
                raise IOError("Could not align and initialize ADCs.")
        
        # Initialize other blocks and both correlators
        print("Initializing other blocks, including PFB and both correlators...")
        try:
            self.s.initialize()
        except:
            self.s.pfb.initialize()
            self.s.corr_0.initialize()
            self.s.corr_1.initialize()

        print('Spectrometer is ready.')

##############################################################

    # NEEDS A LOT OF WORK
    def poll(self):
        """
        Waits until the integration count has been incrimented and
        returns the date of the integration in seconds (s) since 
        Jan 1, 1970 UTC.

        Returns:
        - obs_date: Unix time of the integration.
        """
        self.count0 = self.fpga.read_int('corr_0_acc_cnt')
        while self.fpga.read_int('corr_0_acc_cnt') == self.count0:
            time.sleep(0.1)

        self.count1 = self.fpga.read_int('corr_1_acc_cnt')
        while self.fpga.read_int('corr_1_acc_cnt') == self.count1:
            time.sleep(0.1)

        obs_date = time.time() - 0.5*self.int_time
        self.count0 = self.fpga.read_int('corr_0_acc_cnt')
        self.count1 = self.fpga.read_int('corr_1_acc_cnt')
        return obs_date

        
    # NOT NEEDED ANYMORE???
    def read_bram(self, bram_name):
        """
        Reads out data from a SNAP BRAM. The data is stored in the SNAP
        as 32-bit fixed point numbers with the binary point at the 
        30th bit.

        Inputs:
        - bram: Name of the BRAM to read data from.
        Returns:
        - bram_fp: Array of floats of the SNAP BRAM values.
        """
        bram_size = 3*self.nchan
        bram_ints = np.fromstring(self.fpga.read(bram_name, bram_size), '>i4')

        # Remove DC offset
        bram_ints[0] = 0
        bram_ints[1] = 0
        bram_ints[-1] = 0
        bram_ints[-2] = 0

        bram_fp = bram_ints/float(1<<30)
        return bram_fp


    def read_spec(self, filename, nspec, coords, coord_sys='ga'):
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

        # Make sure the spectrometer is actually running. If not, initialize.
        if not self.s.fpga.is_running():
            self.initialize()

        # Make FITS file
        hdulist = fits.HDUList()

        # Name data columns
        data_names = ['auto0_real', 'auto1_real', 'cross_real', 'cross_imag']

        # Create a primary FITS HDU for a table file
        hdulist = [self.fits_header(nspec, coords, coord_sys)]

        # Read some number of spectra to a FITS file
        print('Reading', nspec, 'spectra from the SNAP.')
        ninteg = 0
        while ninteg < nspec:
            cols = []
            pols = [0,0], [1,1], [0,1] # Polarizations for auto0, auto1, cross
            try:
                auto0 = self.s.corr_0.get_new_corr(0,0)
                auto1 = self.s.corr_1.get_new_corr(1,1)
                auto0_real, auto1_real = auto0.real, auto1.real
                cross = self.s.corr_0.get_new_corr(0,1)
                cross_real, cross_imag = cross.real, cross.imag
                spectra = [auto0_real, auto1_real, cross_real, cross_imag]
                for i in range(len(data_names)):
                    cols.append(fits.Column(name=data_names[i], format='D', array=spectra[i]))

            except RuntimeError:
                print('WARNING: Cannot reach the SNAP. Skipping integration.')
                self.fpga.connect()
                continue

            # Add header and data to FITS
            hdulist.append(self.fits_header(nspec, coords, coord_sys))
            hdulist.append(fits.BinTableHDU.from_columns(cols, name='CORR_DATA'))
            
            # Add to counter
            ninteg += 1
            # integ_time = ninteg*self.int_time
            print('Integration count:', ninteg) #, '(' + str(integ_time), 's)')

        # Save the output file
        print('Saving spectra to output file:', filename)
        hdulist.writeto(filename+".fits", overwrite=True)
        hdulist.close()

    # DON'T KNOW ABOUT THESE NEXT THREE FUNCTIONS
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









