print("Importing dependancies...")

import casperfpga
from hera_corr_f import SnapFengine
import ugradio

import numpy as np
import time 
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from astropy.io import fits


#print("\n\nBuilding spectrometer calss...")

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
        self.fpga = casperfpga.CasperFpga(self.host)
        self.fpga.upload_to_ram_and_program(self.fpgfile)
        #self.adc = casperfpga.snapadc.SNAPADC(self.fpga) # 10 MHz reference signal default
        self.s = SnapFengine(self.host, transport='default')
        self.s.fpga.upload_to_ram_and_program(self.fpgfile)

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
        Checks if the fpga of the spectrometer has been programed 
        and is running. Returns an IOError if the SNAP is not 
        running and cannot program.
        """
        if self.fpga.is_running() and self.s.is_programmed():
            print("Fpga is programmed and running.")
        elif not self.fpga.is_running() or not self.s.is_programmed():
            print("WARNING: Fpga is not running. Programming...")
            self.fpga.upload_to_ram_and_program(self.fpgfile)
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

        header['NSPEC'] = (nspec, 'Number of spectra collected')
        header['FPGFILE'] = (self.fpgfile, 'FPGA FPG file')
        header['HOST'] = (self.s.fpga.host, "Host of the FPGA")
        #header['MODE'] = (self.mode, 'Spectrometer mode')
        header['CLK'] = (self.s.fpga.estimate_fpga_clock(), 'FPGA clock speed [MHz]')
        header['ADC_NAME'] = (self.s.adc.adc.name, "Name of ADC")
        #header['ADC'] = (self.adc_rate, 'ADC clock speed [Hz]')
        #header['DOWNSAMPLE'] = (self.downsample, 'ADC downsampling period')
        #header['SAMPRATE'] = (self.samp_rate, 'Downsampled clock speed [Hz]')
        #header['BW'] = (self.bandwidth, 'Bandwidth of spectra [Hz]')
        #header['NCHAN'] = (self.nchan, 'Number of frequency channels')
        #header['RES'] = (self.resolution, 'Frequency resolution [Hz]')
        #header['FFTSHIFT'] = (self.fft_shift, 'FFT shifting instructions')
        #header['ACCLEN'] = (self.acc_len, 'Number of clock cycles')
        #header['INTTIME'] = (self.int_time, 'Integration time of spectra')
        #header['SCALE'] = (self.scale, 'Average instead of sum on SNAP')

        header['L'] = (l.value, 'Galactic longitude [deg]')
        header['B'] = (b.value, 'Galactic latitude [deg]')
        header['RA'] = (ra.value, 'Right Ascension [deg]')
        header['DEC'] = (dec.value, 'Declination [deg]')
        header['JD'] = (obs_start_jd, 'Julian date of start time')
        header['UNIX'] = (obs_start_unix, 'Seconds since epoch')

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
        Starts the fpg process on the SNAP and initializes the 
        spectrometer.
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

        # Initialize other blocks and correlators
        print("Initializing other blocks, including PFB and both correlators...")
        try:
            self.s.initialize()
        except:
            self.s.pfb.initialize()
            self.s.corr_0.initialize()
            self.s.corr_1.initialize()

        print('Spectrometer is ready.')

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
        if not self.fpga.is_running():
            self.initialize()

        # Make FITS file
        hdulist = fits.HDUList()

        # Name data columns
        data_names = ['auto0_real', 'auto1_real', 'cross_real', 'cross_imag']

        # Read some number of spectra to a FITS file
        print('Reading', nspec, 'spectra from the SNAP.')
        ninteg = 0
        while ninteg < nspec:
            cols = []
            pols = [0,0], [1,1], [0,1] # Polarizations for auto0, auto1, cross
            try:
                auto0 = self.s.corr_0.get_new_corr(*pols[0])
                auto0_real = auto0.real
                
                auto1 = self.s.corr_1.get_new_corr(*pols[1])
                auto1_real = auto0.real

                cross = self.s.corr_0.get_new_corr(*pols[2])
                cross_real, cross_imag = cross.real, cross_imag

                spectra = [auto0_real, auto1_real, cross_real, cross_imag]
                for i in range(len(data_names)):
                    cols.append(fits.Column(name=data_names[i],format='D', array=spectra[i]))

            except RuntimeError:
                print("WARNING: Cannot reach the SNAP. Skipping integration.")
                self.fpga.connect()
                continue

            # Add header and data to FITS
            hdulist.append(self.fits_header(nspec, coords, coord_sys))
            hdulist.append(fits.BinTableHDU.from_columns(cols, name='CORR_DATA'))

            # Add to counter
            ninteg += 1
            #integ_time = ninteg*self.int_time
            print("Integration count:", ninteg) #, "(" + str(integ_time), "s)")

        # Save the output file:
        print("Saving spectra to output file: ", filename+".fits")
        hdulist.writeto(filename+".fits", overwrite=True)
        hdulist.close()


# print("Spectrometer has been built.")

# Open IPython to continue session
import IPython
Q = input("\n\nDo you wish to continue in IPython? [y/n] ")
if Q == "y":
    print("\nOpening IPython...:")
    IPython.embed()
elif Q == "n":
    print("Goodbye.")

