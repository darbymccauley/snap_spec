import casperfpga
from hera_corr_f import SnapFengine
import ugradio

import numpy as np
import time 
import logging
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from astropy.io import fits
import os, sys

DELAY_TIME = 0.1 # seconds
HOST = 'localhost'
FPGFILE = 'fpga/ugradio_corrspec_2022-02-22_0905.fpg'
STREAM_1 = 0
STREAM_2 = 1
LOGGER = 'spectrometer.log'
TRANSPORT = 'default'

# Create Spectrometer class
class Spectrometer(object):
    """
    Casperfpga interface to the SNAP spectrometer.
    """

    def __init__(self, host=HOST, fpgfile=FPGFILE, transport=TRANSPORT, stream_1=STREAM_1, stream_2=STREAM_2, logger=None):
        """
        Create the interface to the SNAP.

        Inputs:
        - host: IP address of SNAP.
        - fpgfile: design file used to program fpga.
        - transport: communication protocal.
        - stream_1, stream_2: SNAP ports used for correlation data aquisition.
        - logger: filename in which log is recorded. 
        """
        self.host = host
        self.fpgfile = fpgfile
        self.transport = transport

        if logger is None:
            self.logger = LOGGER
        elif logger is not None:
            self.logger = logger
        logging.basicConfig(filename=self.logger, 
                            format='%(asctime)s - %(levelname)s - %(message)s', 
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            level=logging.NOTSET)

        # Ports used for ADCs
        self.stream_1 = stream_1
        self.stream_2 = stream_2
        
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
        self.s = SnapFengine(self.host, transport=self.transport)

        
    def is_connected(self):
        """
        Check if the SNAP is connected.
        """
        if self.fpga.is_connected():
            return True
        else:
            logging.warning('SNAP is not connected')
            return False


    def is_running(self):
        """
        Check if the fpga has been programmed and is running.
        """
        if self.fpga.is_running() and self.s.is_programmed():
            return True
        else:
            logging.warning('SNAP is not programmed and running.')

  
    def program(self):
        """
        Program the fpga.
        """
        self.fpga.upload_to_ram_and_program(self.fpgfile)
        self.s.fpga.upload_to_ram_and_program(self.fpgfile)
        

    def initialize(self):
        """
        Programs the fpga on the SNAP and initializes the spectrometer.
        """
        logging.info('Starting the spectrometer.')
        
        # Program fpga
        self.program()
        
        # Initialize and align ADCs
        logging.info('Aligning and initializing ADCs...')
        try:
            self.s.adc.init()
            self.s.align_adc()
            logging.info('ADCs aligned and initialized.')
        except: ### No bare excepts
            try: # try again (usually works after two attempts)
                self.s.adc.init()
                self.s.align_adc()
                logging.info('ADCs aligned and initialized.')
            except:
                logging.error('Could not align and initialize ADCs.')
                raise IOError('Could not align and initialize ADCs.')

        # Initialize other blocks and both correlators
        logging.info('Initializing other blocks, including PFB and both correlators.')
        try:
            self.s.initialize()
        except:
            self.s.pfb.initialize()
            self.s.corr_0.initialize()
            self.s.corr_1.initialize()
        logging.info('Spectrometer is ready.')


    def make_PrimaryHDU(self, nspec, coords, coord_sys='ga'):
        """
        Make the PrimaryHDU of the FITS file. Serves as the header 
        containing metadata of the system as well as observation 
        attributes.
        
        Inputs:
        - nspec: Number of spectra to collect.
        - coords: Coordinate(s) of the target.
            Format: (l/ra, b/dec)
        - coord_sys: Coordinate system used for ''coords''.
            Default is galactic coordinates. Takes in either galactic 
            ('ga') or equatorial ('eq') coordinate systems.
        Returns:
        - PrimaryHDU information containing the attributes 
        of the observation and spectrometer.
        """
        # Ensure that a proper coordinate system has been supplied
        if coord_sys != 'ga' and coord_sys != 'eq':
            raise ValueError("Invalid coordinate system supplied: " + coord_sys)

        # Set times
        obs_start_unix = time.time() #unix time
        unix_object = Time(obs_start_unix, format='unix', 
                           location=(ugradio.leo.lon, ugradio.leo.lat, ugradio.leo.alt)) #unix time Time object
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
        
        # Make PrimaryHDU
        header = fits.Header()

        # Save metadata of the system and spectrometer
        header['NSPEC'] = (nspec, "Number of spectra collected")
        header['FPGFILE'] = (self.fpgfile, "FPGA FPG file")
        header['HOST'] = (self.s.fpga.host, "Host of the FPGA")
        header['TRANSPORTT'] = (self.transport, "Communication protocal")
        # header['CLK'] = (self.s.fpga.estimate_fpga_clock(), "FPGA clock speed [MHz]")
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
        header['PYTHON'] = (3.8, "Python version")
        header['SRC'] = ('https://github.com/darbymccauley/Leuschner_Spectrometer.git', "Source code")
        # header['CASPERFPGA'] = (CASPERFPGA_VERSION, "casperfpga code used")
        # header['HERA_CORR_F'] = (HERA_CORR_F_VERSION, "hera_corr_f code used")
        
        # Save observation attributes
        header['L'] = (l.value, "Galactic longitude [deg]")
        header['B'] = (b.value, "Galactic latitude [deg]")
        header['RA'] = (ra.value, "Right Ascension [deg]")
        header['DEC'] = (dec.value, "Declination [deg]")
        header['JD'] = (obs_start_jd, "Julian date of start time")
        header['UNIX'] = (obs_start_unix, "Seconds since epoch")

        primaryhdu = fits.PrimaryHDU(header=header)
        return primaryhdu


    def read_spec(self, filename, nspec, coords, coord_sys='ga'):
        """
        Recieves spectrometer data from the Leuschner spectrometer and 
        saves it to a FITS file. The primary HDU contains information about
        the observation (coordinates, number of spectra collected, time,
        etc.) and spectrometer attributes used. Each set of spectra is
        stored in its own FITS table in the FITS file. The columns in
        each FITS table are ''auto0_real'', ''auto1_real'',
        ''cross_real'', and ''cross_imag''. All columns contain
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
        - FITS file with autocorrelated spectrometer data.
        """
        # Make PrimaryHDU for FITS file
        primaryhdu = self.make_PrimaryHDU(nspec, coords, coord_sys)
        hdulist = fits.HDUList(hdus=[primaryhdu])

        # Define spectra to collect
        spectra = [('auto0_real', self.s.corr_0, (self.stream_1, self.stream_1)), # (0, 0)
                   ('auto1_real', self.s.corr_1, (self.stream_2, self.stream_2))] # (1, 1)
        data = {}
   
        cnt_0 = self.s.corr_0.read('acc_len', 4) # Read corr_0 count
        for ninteg in range(nspec):
            while self.s.corr_0.read('acc_cnt', 4) == cnt_0: # sleep until corr_0 count increased
                time.sleep(DELAY_TIME)
            for name, corr, (stream_1, stream_2) in spectra: # read the spectra from both corrs
                data[name] = corr.get_new_corr(stream_1, stream_2).real
            assert cnt_0 + 1 == self.s.corr_1.read('acc_cnt', 4) # assert corr_0's count increased by 1
            cnt_0 = cnt_0 + 1
            # Make BinTableHDU and append collected data
            data_list = [fits.Column(name=name, format='D', array=data[name]) for name, _, _ in spectra]
            bintablehdu = fits.BinTableHDU.from_columns(data_list, name='CORR_DATA')
            hdulist.append(bintablehdu)
           
        # Save the output file
        hdulist.writeto(filename, overwrite=True)
        hdulist.close()


    def read_corr(self, filename, nspec, coords, coord_sys='ga'):
        """
        Recieves correlation data from the Leuschner spectrometer and 
        saves it to a FITS file. The primary HDU contains information about
        the observation (coordinates, number of spectra collected, time,
        etc.) and spectrometer attributes used. Each set of spectra is
        stored in its own FITS table in the FITS file. The columns in
        each FITS table are ''auto0_real'', ''auto1_real'',
        ''cross_real'', and ''cross_imag''. All columns contain
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
        - FITS file with correlated spectrometer data.
        """
        primaryhdu = self.make_PrimaryHDU(nspec, coords, coord_sys)
        hdulist = fits.HDUList(hdus=[primaryhdu])

        # Read some number of spectra to a FITS file
        ninteg = 0
        while ninteg < nspec:
            spectra = [('cross', (self.stream_1, self.stream_2))] # (0, 1)
            data_list = []
            for name, (stream_1, stream_2) in spectra:
                cross = self.s.corr_0.get_new_corr(stream_1, stream_2)
                cross_real, cross_imag = cross.real, cross.imag
            
            data_list.append(fits.Column(name=name+'_real', format='D', array=cross_real))
            data_list.append(fits.Column(name=name+'_imag', format='D', array=cross_imag))

            bintablehdu = fits.BinTableHDU.from_columns(data_list, name='CORR_DATA')
            hdulist.append(bintablehdu)
            ninteg += 1
 
        # Save the output file
        hdulist.writeto(filename, overwrite=True)
        hdulist.close()
