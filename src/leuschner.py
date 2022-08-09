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
import tqdm
import random


DELAY_TIME = 0.1 # seconds
HOST = 'localhost'
FPGFILE = 'fpga/ugradio_corrspec_2022-02-22_0905.fpg'
STREAM_1 = 0
STREAM_2 = 1
LOGGER = 'spectrometer.log'
TRANSPORT = 'default'
ACC_LEN = 38150
SPEC_PER_ACC = 8

# Create Spectrometer class
class Spectrometer:
    """
    Casperfpga interface to the SNAP spectrometer.
    """

    def __init__(self, host=HOST, fpgfile=FPGFILE, transport=TRANSPORT, stream_1=STREAM_1, stream_2=STREAM_2, logger=None, acc_len=ACC_LEN, spec_per_acc=SPEC_PER_ACC):
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
                            level=logging.WARNING)

        # Ports used for ADCs
        self.stream_1 = stream_1
        self.stream_2 = stream_2
        
        self.acc_len = acc_len
        self.spec_per_acc = spec_per_acc

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
            return False

  
    def program(self):
        """
        Program the fpga.
        """
        self.fpga.upload_to_ram_and_program(self.fpgfile)
        self.s.fpga.upload_to_ram_and_program(self.fpgfile)
        logging.info('FPGA programmed.')

 
    def alignFrameClock_discover_snap(self, chipsel=None, chip_lanes=None, retry=True):
        """
        Align frame clock with data frame. (Sourced from hera_corr_f.)
        
        Inputs: 
        - chipsel (list): which ADC chip to align. Default is all chips. 
        - chips_lanes: which lanes of the ADC(s) to align. Default is all lanes.
        - retry: whether to retry alignment if initial attempt fails.

        Returns:
        - Which chips failed and the corresponding lanes that errored.
        """
        if chip_lanes is None:
                if chipsel is None:
                    chip_lanes = {chip:self.s.adc.laneList for chip in self.s.adc.adcList}
                elif chipsel is not None:
                    chip_lanes = {chip:self.s.adc.laneList for chip in chipsel}
        logging.debug('Aligning frame clock on ADCs/lanes: %s' % \
                        str(chip_lanes))

        failed_chips = {}
        self.s.adc.setDemux(numChannel=1)
        for chip, lanes in chip_lanes.items():
            self.s.adc.selectADC(chip)
            self.s.adc.adc.test('dual_custom_pat', self.s.adc.p1, self.s.adc.p2)
            ans1 = self.s.adc._signed(self.s.adc.p1, self.s.adc.RESOLUTION)
            ans2 = self.s.adc._signed(self.s.adc.p2, self.s.adc.RESOLUTION)
            failed_lanes = []
            for cnt in range(2*self.s.adc.RESOLUTION):
                slipped = False
                self.s.adc.snapshot() # make bitslip "take" (?!) XXX
                d = self.s.adc.readRAM(chip).reshape(-1, self.s.adc.RESOLUTION)
                # sanity check: these failures mean line clock errors
                failed_lanes += [L for L in lanes
                        if np.any(d[0::2,L] != d[0,L]) or \
                           np.any(d[1::2,L] != d[1,L])]
                lanes = [L for L in lanes if L not in failed_lanes]
                for lane in lanes:
                    if not d[0, lane] in [ans1, ans2]:
                        if cnt == 2*self.s.adc.RESOLUTION - 1:
                            # Failed on last try
                            failed_lanes += [lane]
                        self.s.adc.bitslip(chip, lane)
                        slipped = True
                if not slipped:
                    break
            self.s.adc.adc.test('off')
            if len(failed_lanes) > 0:
                failed_chips[chip] = failed_lanes
        self.s.adc.setDemux(numChannel=self.s.adc.num_chans)
        if len(failed_chips) > 0 and retry:
            if self.s.adc._retry_cnt < self.s.adc._retry:
                self.s.adc._retry_cnt += 1
                logging.info('retry=%d/%d redo Line on ADCs/lanes: %s' % \
                            (self.s.adc._retry_cnt, self.s.adc._retry, failed_chips))
                self.alignLineClock_discover_snap(chipsel=failed_chips) # retry using my functions again
                return self.alignFrameClock_discover_snap(chipsel=failed_chips) # retry using my functions again
        return failed_chips


    def alignLineClock_discover_snap(self, chipsel=None, chip_lanes=None, ker_size=5):
        """
        Find a tap for the line clock that produces reliable bit
        capture from ADC.

        Inputs: 
        - chipsel (list): which ADC chip to align. Default is all chips. 
        - chips_lanes: which lanes of the ADC(s) to align. Default is all lanes.
        - retry: whether to retry alignment if initial attempt fails.

        Returns:
        - Which chips failed and the corresponding lanes that errored.        
        """
        if chip_lanes is None:
                if chipsel is None:
                    chip_lanes = {chip:self.s.adc.laneList for chip in self.s.adc.adcList}
                elif chipsel is not None:
                    chip_lanes = {chip:self.s.adc.laneList for chip in chipsel}
        logging.info('Aligning line clock on ADCs/lanes: %s' % \
                        str(chip_lanes))
        try:
            self.s.adc._find_working_taps(ker_size=ker_size)
        except(RuntimeError):
            logging.info('Failed to find working taps.')
            return chip_lanes # total failure
        self.s.adc.setDemux(numChannel=1)
        for chip, lanes in chip_lanes.items():
            self.s.adc.selectADC(chip)
            taps = self.s.adc.working_taps[chip]
            tap = random.choice(taps)
            for L in self.s.adc.laneList: # redo all lanes to be the same
                self.s.adc.delay(tap, chip, L)
            # Remove from future consideration if tap doesn't work out
            self.s.adc.working_taps[chip] = taps[np.abs(taps - tap) >= ker_size//2]
            logging.info('Setting ADC=%d tap=%s' % (chip, tap))
        self.s.adc.setDemux(numChannel=self.s.adc.num_chans)
        return {} # success


    def rampTest_discover_snap(self, chipsel=None, nchecks=300, retry=False):
        """
        (Sourced from hera_corr_f.)
        
        Inputs:
        - chipsel (list): which ADC chip to align. Default is all chips. 
        - nchecks: number of times to check ramp test passing.
        - retry: whether to retry if initial attempt fails.

        Returns:
        - Which chips failed and the corresponding lanes that errored.
        """
        if chipsel is None:
            chips = self.s.adc.adcList
        elif chipsel is not None:
            chips = chipsel
        logging.debug('Ramp test on ADCs: %s' % str(chips))
        failed_chips = {}
        self.s.adc.setDemux(numChannel=1)
        predicted = np.arange(128).reshape(-1,1)
        self.s.adc.selectADC(chips) # specify chips
        self.s.adc.adc.test("en_ramp")
        for cnt in range(nchecks):
            self.s.adc.snapshot()
            for chip, d in self.s.adc.readRAM(ram=chips, signed=False).items():
                ans = (predicted + d[0,0]) % 256
                failed_lanes = np.sum(d != ans, axis=0)
                if np.any(failed_lanes) > 0:
                    failed_chips[chip] = np.where(failed_lanes)[0]
            if (retry is False) and len(failed_chips) > 0:
                # can bail out if we aren't retrying b/c we don't need list of failures.
                break
        self.s.adc.selectADC(chips) # specify chips
        self.s.adc.adc.test('off')
        self.s.adc.setDemux(numChannel=self.s.adc.num_chans)
        if len(failed_chips) > 0 and retry:
            if self.s.adc._retry_cnt < self.s.adc._retry:
                self.s.adc._retry_cnt += 1  
                logging.info('retry=%d/%d redo Line/Frame on ADCs/lanes: %s' % \
                                (self.s.adc._retry_cnt, self.s.adc._retry, failed_chips))
                self.alignLineClock_discover_snap(failed_chips) # retry using my functions again
                self.alignFrameClock_discover_snap(failed_chips) # retry using my functions again
                return self.rampTest_discover_snap(chipsel=chipsel, nchecks=nchecks, retry=retry) # retry using my functions again
        return failed_chips


    def align_adc_discover_snap(self, chipsel=None, chip_lanes=None, ker_size=5, retry=True, nchecks=300, force=False, verify=True):
        """
        Align clock and data lanes of ADC. (Sourced from hera_corr_f.)
        """
        if force:
            self.s._set_adc_status(0)
        if self.s.adc_is_configured():
            return
        fails = self.alignLineClock_discover_snap(chipsel=chipsel, chip_lanes=chip_lanes, ker_size=ker_size)
        if len(fails) > 0:
            logging.warning("alignLineClock failed on: " + str(fails))
        fails = self.alignFrameClock_discover_snap(chipsel=chipsel, chip_lanes=chip_lanes, retry=retry)
        if len(fails) > 0:
            logging.warning("alignFrameClock failed on: " + str(fails))
        fails = self.rampTest_discover_snap(chipsel=chipsel, nchecks=nchecks, retry=False)
        if len(fails) > 0:
            logging.warning("rampTest failed on: " + str(fails))
        else:
            self.s._set_adc_status(1)  # record status
        if verify:
            assert(self.s.adc_is_configured())  # errors if anything failed
        # Otherwise, finish up here.
        self.s.adc.selectADC(chipsel)
        self.s.adc.adc.selectInput([1, 1, 3, 3])
        self.s.adc.set_gain(4)
        

    def initialize_discover_SNAP(self, chipsel=None, chip_lanes=None, ker_size=5, retry=True, nchecks=300, force=False, verify=True):
        """
        Program the fpga on the SNAP and initialize the 1st ADC on the
        Discover SNAP board.
        """
        if chipsel is None:
            self.initialize()
        
        logging.info('Initializing SNAP...')

        # Program fpga
        self.program()

        self.s.corr_0.set_acc_len(self.acc_len)
        self.s.corr_1.set_acc_len(self.acc_len)

        # Initialize and align the Nth ADC
        while self.s.adc_is_configured() == 0:
            self.s.adc.init()
            self.align_adc_discover_snap(chipsel=chipsel,
                                        chip_lanes=chip_lanes, 
                                        ker_size=ker_size, 
                                        retry=retry, 
                                        nchecks=nchecks, 
                                        force=force, 
                                        verify=verify)
        
        # Initialize other blocks and both correlators
        try:
            self.s.initialize()
        except UnicodeDecodeError: # XXX address this issue later with Aaron
            self.s.pfb.initialize()
            self.s.corr_0.initialize()
            self.s.corr_1.initialize()
        logging.info('Spectrometer initialized.')


    def initialize(self):
        """
        Programs the fpga on the SNAP and initializes the spectrometer.
        """
        logging.info('Initializing the spectrometer...')
        
        # Program fpga
        self.program()

        self.s.corr_0.set_acc_len(self.acc_len)
        self.s.corr_1.set_acc_len(self.acc_len)
        
        # Initialize and align ADCs
        while self.s.adc_is_configured() == 0:
            self.s.adc.init()
            self.s.align_adc()        

        # Initialize other blocks and both correlators
        try:
            self.s.initialize()
        except UnicodeDecodeError: # XXX address this issue later with Aaron
            self.s.pfb.initialize()
            self.s.corr_0.initialize()
            self.s.corr_1.initialize()
        logging.info('Spectrometer initialized.')


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
        header['NSPEC'] = (nspec, 'Number of spectra collected')
        header['FPGFILE'] = (self.fpgfile, 'FPGA FPG file')
        header['HOST'] = (self.s.fpga.host, 'Host of the FPGA')
        header['TRANSP'] = (self.transport, 'Communication protocal (transport)')
        header['ACCLEN'] = (self.acc_len, 'Number of clock cycles')
        header['SPEC/ACC'] = (self.spec_per_acc, 'Spectra per accumulation')
        header['STREAM_1'] = (self.stream_1, 'First ADC port used')
        header['STREAM_2'] = (self.stream_2, 'Second ADC port used')
        #header['LOGGER'] = (self.logger, 'Logger file')

        header['PYTHON'] = (3.8, 'Python version')
        header['SRC'] = ('https://github.com/darbymccauley/Leuschner_Spectrometer.git', 'Source code')
        # header['CASPERFPGA'] = (CASPERFPGA_VERSION, 'casperfpga code used')
        # header['HERA_CORR_F'] = (HERA_CORR_F_VERSION, 'hera_corr_f code used')
        
        # Save observation attributes
        header['L'] = (l.value, 'Galactic longitude [deg]')
        header['B'] = (b.value, 'Galactic latitude [deg]')
        header['RA'] = (ra.value, 'Right Ascension [deg]')
        header['DEC'] = (dec.value, 'Declination [deg]')
        header['JD'] = (obs_start_jd, 'Julian date of start time')
        header['UNIX'] = (obs_start_unix, 'Seconds since epoch')

        primaryhdu = fits.PrimaryHDU(header=header)
        return primaryhdu


    def wait_for_cnt(self):
        """
        Waits for corr_0 acc_cnt to increase by 1. 
        Returns the count read from the register.
        (Sourced with modifications from hera_corr_f.)
        """
        cnt_0 = self.s.corr_0.read_uint('acc_cnt')
        while self.s.corr_0.read_uint('acc_cnt') < (cnt_0+1):
            time.sleep(0.1)
        return cnt_0


    def get_new_corr(self, corr, pol1, pol2):
        """
        Reads and returns the spectra collected given a set of polarizations.

        Inputs:
        - corr: which correlator to use
        - pol1: first polarization
        - pol2: second polarization

        Returns: correlated data, either auto or cross depending on choice of pol1 and pol2.

        (Sourced with modifications from hera_corr_f.)
        """
        corr.set_input(pol1, pol2)
        spec = corr.read_bram(flush_vacc=False)/float(self.acc_len*self.spec_per_acc)
        if pol1 == pol2:
            return spec.real + 1j*np.zeros(len(spec))
        else:
            return spec


    def read_spec(self, filename, nspec, coords, coord_sys='ga', progress=False):
        """
        Recieves spectrometer data from the Leuschner spectrometer and 
        saves it to a FITS file. The primary HDU contains information about
        the observation (coordinates, number of spectra collected, time,
        etc.) and spectrometer attributes used. Each set of spectra is
        stored in its own FITS table in the FITS file. The columns in
        each FITS table are ''auto0_real'' and ''auto1_real'',
        for each polarization's auto-correlation. All columns contain
        double-precision floating-point numbers.

        Inputs:
        - filename: Name of the output FITs file.
        - nspec: Number of spectra to collect.
        - coords: Coordinate(s) of the target.
            Format: (l/ra, b/dec)
        - coord_sys: Coordinate system used for ''coords''.
            Default is galactic coordinates. Takes in either galactic 
            ('ga') or equatorial ('eq') coordinate systems.
        - progress: shows progress bar for data collection.
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
        
        # Collect spectra
        if progress is True:
            logging.info('Reading %s spectra from the SNAP' % str(nspec))
            for ninteg in tqdm.tqdm(range(nspec), desc='Progress'):
                cnt_0 = self.wait_for_cnt()        
                for name, corr, (stream_1, stream_2) in spectra: # read the spectra from both corrs
                    data[name] = self.get_new_corr(corr, stream_1, stream_2).real
                cnt_1 = self.s.corr_1.read_uint('acc_cnt')
                assert cnt_0 + 1 == cnt_1 # assert corr_0's count increased and matches corr_1's count

                # Make BinTableHDU and append collected data
                data_list = [fits.Column(name=name, format='D', array=data[name]) for name, _, _ in spectra]
                bintablehdu = fits.BinTableHDU.from_columns(data_list, name='CORR_DATA')
                hdulist.append(bintablehdu)

        if progress is False:
            logging.info('Reading %s spectra from the SNAP' % str(nspec))
            for ninteg in range(nspec):
                cnt_0 = self.wait_for_cnt()
                for name, corr, (stream_1, stream_2) in spectra: # read the spectra from both corrs
                    data[name] = self.get_new_corr(corr, stream_1, stream_2).real
                cnt_1 = self.s.corr_1.read_uint('acc_cnt')
                assert cnt_0 + 1 == cnt_1 # assert corr_0's count increased and matches corr_1's count

                # Make BinTableHDU and append collected data
                data_list = [fits.Column(name=name, format='D', array=data[name]) for name, _, _ in spectra]
                bintablehdu = fits.BinTableHDU.from_columns(data_list, name='CORR_DATA')
                hdulist.append(bintablehdu)
           
        # Save the output file
        hdulist.writeto(filename, overwrite=True)
        hdulist.close()


    def read_corr(self, filename, nspec, coords, coord_sys='ga', progress=False):
        """
        Recieves correlation data from the Leuschner spectrometer and 
        saves it to a FITS file. The primary HDU contains information about
        the observation (coordinates, number of spectra collected, time,
        etc.) and spectrometer attributes used. Each set of spectra is
        stored in its own FITS table in the FITS file. The columns in
        each FITS table are ''cross_real'' and ''cross_imag'', for the real
        and imaginary components of the cross-correlated data. All columns 
        contain double-precision floating-point numbers.

        Inputs:
        - filename: Name of the output FITs file.
        - nspec: Number of spectra to collect.
        - coords: Coordinate(s) of the target.
            Format: (l/ra, b/dec)
        - coord_sys: Coordinate system used for ''coords''.
            Default is galactic coordinates. Takes in either galactic 
            ('ga') or equatorial ('eq') coordinate systems.
        - progress: shows progress bar for data collection.
        
        Returns:
        - FITS file with correlated spectrometer data.
        """
        primaryhdu = self.make_PrimaryHDU(nspec, coords, coord_sys)
        hdulist = fits.HDUList(hdus=[primaryhdu])

        if progress is True:
            logging.info('Reading %s spectra from the SNAP' % str(nspec))
            pbar = tqdm.tqdm(total=nspec, desc='Progress')
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
                pbar.update()
            pbar.close()

        if progress is False:
            logging.info('Reading %s spectra from the SNAP' % str(nspec))
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


# Things to address:
#### UnicodeDecodeError in s.initialize()
#### Other things to be added to PrimaryHDU
#### logging stuff
