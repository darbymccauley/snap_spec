import casperfpga
from hera_corr_f import SnapFengine
from hera_corr_f import blocks as snap_blocks
#import ugradio

import numpy as np
import time 
import logging
#from astropy.coordinates import SkyCoord
#import astropy.units as u
#from astropy.time import Time
#from astropy.io import fits
#import tqdm
import os


DELAY_TIME = 0.1 # seconds
HOST = 'localhost'
FPGFILE = os.path.join(os.path.dirname(__file__), 'fpga', 'ugradio_corrspec_2022-02-22_0905.fpg')
STREAM_1 = 0
STREAM_2 = 1
LOGGER = 'spectrometer.log'
TRANSPORT = 'default'
ACC_LEN = 38150
SPEC_PER_ACC = 8
CHIP_SEL_DISCOVER_SNAP = [0] # for the Discover SNAP spectrometer


class UGRadioSnap(SnapFengine):
    """
    Interface to the SNAP.

    Arguments:
    - host: IP address of SNAP.
    - fpgfile: design file used to program fpga.
    - stream_1, stream_2: SNAP ports used for correlation data.
    - is_discover: bool describing if the discover snap is being used.
    - acc_len: accumulation length.
    - spec_per_acc: number of spectra collected per accumulation
    """
    def __init__(self,
                 host=HOST,
                 fpgfile=FPGFILE,
                 is_discover=True,
                 stream_1=STREAM_1, 
                 stream_2=STREAM_2, 
                 acc_len=ACC_LEN, 
                 spec_per_acc=SPEC_PER_ACC):
        try:
            super().__init__(host, ant_indices=None, logger=None, 
                             transport=TRANSPORT, redishost=None)
        except casperfpga.transport_katcp.KatcpRequestFail:
            fpga = casperfpga.CasperFpga(host)
            fpga.upload_to_ram_and_program(fpgfile)
            super().__init__(host, ant_indices=None, logger=None, 
                             transport=TRANSPORT, redishost=None)
            

        self.host = host
        self.fpgfile = fpgfile
        self.transport = TRANSPORT

        # Ports used for ADCs
        self.stream_1 = stream_1
        self.stream_2 = stream_2
        self.acc_len = acc_len
        self.spec_per_acc = spec_per_acc

        self.corr_0 = self.corr
        self.corr_1 = snap_blocks.Corr(self.fpga, 'corr_1')

        # blocks initialized in this (significant) order
        self.blocks = [
            self.synth,
            self.adc,
            self.sync,
            self.noise,
            self.input,
            self.delay,
            self.pfb,
            self.eq,
            # self.eq_tvg, # temporarily removed // not needed for spec
            #self.reorder, # not needed for spec
            #self.packetizer, # not needed for spec
            #self.eth,
            self.corr_0,
            self.corr_1,
            # self.phase_switch XXX katcp error -- address with Aaron
        ]

        if is_discover:
            self.chips = CHIP_SEL_DISCOVER_SNAP
        else:
            self.chips = None

    def _add_i2c(self):
        '''Bypasses initialization of HERA FEM/PAM which are not
        in this system.'''
        self.i2c_initialized = True


    def initialize(self, mode='corr', force=False, sample_rate=500., verify=False):
        """
        Do all programming, initialization, and synchronization.
        """
        logging.info('Initializing the spectrometer...')

        if not self.is_programmed() or force:
            self.program(self.fpgfile, force=force)
            self.initialize_adc(sample_rate=sample_rate)
            self.align_adc(force=force)

        super().initialize(force=force, verify=verify)
        self.corr_0.set_acc_len(self.acc_len)
        self.corr_1.set_acc_len(self.acc_len)
        
        # Set FFT shift
        self.pfb.set_fft_shift(0xfff)

        # Set Eq. coeffs
        value = 160 # debug session with Aaron -- may need revision
        cos = np.ones(self.eq.ncoeffs) * value
        for stream in [self.stream_1, self.stream_2]:
            self.eq.set_coeffs(stream=stream, coeffs=cos)
        
        self.synchronize()

        assert mode in ('spec', 'corr')
        self.mode = mode
        if mode == 'spec':
            # Define spectra to collect
            self.corr_0.set_input(self.stream_1, self.stream_1)  # (0, 0)
            self.corr_1.set_input(self.stream_2, self.stream_2)  # (1, 1)
        else:  # corr mode
            self.corr_0.set_input(self.stream_1, self.stream_2)  # (0, 1)
            # corr_1 unused
        logging.info('Spectrometer initialized.')


    def synchronize(self, manual=True):
        """
        Synchronize DSP logic.
        Called by self.startup()
        """
        self.sync.set_delay(0)
        if manual:
            self.sync.arm_sync()
            for i in range(3):
                self.sync.sw_sync()
        else:
            self.sync.wait_for_sync()
            self.sync.arm_sync()


    def align_adc(self, chips=None, force=False, verify=True):
        """
        Align clock and data lanes of ADC.
        Called by self.startup()
        """
        # repeating this code to bypass non-functioning adcs on discover snap
        if chips is None:
            chips = self.chips
        if force:
            self._set_adc_status(0)
        if self.adc_is_configured():
            return
        chips_lanes = {chip:self.adc.laneList for chip in chips}
        fails = self.adc.alignLineClock(chips_lanes=chips_lanes)
        if len(fails) > 0:
            self.logger.warning("alignLineClock failed on: " + str(fails))
        fails = self.adc.alignFrameClock(chips_lanes=chips_lanes)
        if len(fails) > 0:
            self.logger.warning("alignFrameClock failed on: " + str(fails))
        fails = self.adc.rampTest(chips=list(chips_lanes.keys()))
        if len(fails) > 0:
            self.logger.warning("rampTest failed on: " + str(fails))
        else:
            self._set_adc_status(1)  # record status
        if verify:
            assert(self.adc_is_configured())  # errors if anything failed
        # Otherwise, finish up here.
        self.adc.selectADC()
        self.adc.adc.selectInput([1, 1, 3, 3]) 
        self.adc.set_gain(4)


    def read_data(self, prev_cnt=None):
        """
        Read auto/cross spectra, based on mode set in set_spec_corr_mode.
        Arguments:
            prev_cnt [int]: previous acc_cnt. If provided, wait for current
                            acc_cnt to advance past prev_cnt before proceeding.
        Returns:
            rv [dict]: dict of spectra (auto0/1 or corr01), along with acc_cnt
                       time of data acquisition.
        """
        acc_cnt = self.corr_0.read_uint('acc_cnt')
        if prev_cnt != None:
            while acc_cnt == prev_cnt:
                time.sleep(DELAY_TIME)
                acc_cnt = self.corr_0.read_uint('acc_cnt')
            assert acc_cnt == prev_cnt + 1  # errors if missed an integration
        t = time.time()
        spec0 = self.corr_0.read_bram(flush_vacc=False) / float(self.acc_len * self.spec_per_acc)
        if self.mode == 'spec':
            spec1 = self.corr_1.read_bram(flush_vacc=False) / float(self.acc_len * self.spec_per_acc)
            assert self.corr_1.read_uint('acc_cnt') == acc_cnt  # error if integration advances mid read
            rv = {'auto0': spec0.real, 'auto1': spec1.real}
        else:  # corr mode
            assert self.corr_0.read_uint('acc_cnt') == acc_cnt  # error if integration advances mid read
            rv = {'corr01': spec0}
        rv.update({'acc_cnt': acc_cnt, 'time': t})
        return rv
        

#LEUSCH_LAT, LEUSCH_LON, LEUSCH_ALT = ugradio.leo.lat, ugradio.leo.lon, ugradio.leo.alt
#NCH_LAT, NCH_LON, NCH_ALT = ugradio.nch.lat, ugradio.nch.lon, ugradio.nch.alt
#
#
#class Spectrometer(LeuschFengine):
#    """
#    Casperfpga interface to the SNAP spectrometer.
#    Create the interface to the SNAP.
#
#    Inputs:
#    - host: IP address of SNAP.
#    - fpgfile: design file used to program fpga.
#    - transport: communication protocal.
#    - stream_1, stream_2: SNAP ports used for correlation data aquisition.
#    - logger: filename in which log is recorded.
#    - is_discover: boolean describing if the discover snap is being used.
#    - location: location of observation. Either 'Leuschner' or 'NCH' (New Campbell Hall).
#        (Default is 'Leuschner'.)
#    - acc_len: accumulation length.
#    - spec_per_acc: number of spectra collected per accumulation
#    """
#    def __init__(self,
#                 host=HOST,
#                 fpgfile=FPGFILE,
#                 is_discover=False,
#                 stream_1=STREAM_1, 
#                 stream_2=STREAM_2, 
#                 acc_len=ACC_LEN, 
#                 spec_per_acc=SPEC_PER_ACC,
#                 location='NCH'):
#        super().__init__(host, fpgfile=fpgfile, is_discover=is_discover,
#                         stream_1=stream_1, stream_2=stream_2,
#                         acc_len=acc_len, spec_per_acc=spec_per_acc)
#
#        self.location = location
#        if self.location == 'Leuschner':
#            self.LAT, self.LON, self.ALT = LEUSCH_LAT, LEUSCH_LON, LEUSCH_ALT
#        elif self.location == 'NCH':
#            self.LAT, self.LON, self.ALT = NCH_LAT, NCH_LON, NCH_ALT
#
#        if is_discover:
#            self.snap = 'Discover'
#        elif not is_discover:
#            self.snap = 'NAN'
#        
#        self.fpga.upload_to_ram_and_program(fpgfile)        
#
#    def make_PrimaryHDU(self, nspec, coords, coord_sys='ga'):
#        """
#        Make the PrimaryHDU of the FITS file. Serves as the header 
#        containing metadata of the system as well as observation 
#        attributes.
#        
#        Inputs:
#        - nspec: Number of spectra to collect.
#        - coords: Coordinate(s) of the target.
#            Format: (l/ra, b/dec)
#        - coord_sys: Coordinate system used for ''coords''.
#            Default is galactic coordinates. Takes in either galactic 
#            ('ga') or equatorial ('eq') coordinate systems.
#        Returns:
#        - PrimaryHDU information containing the attributes 
#        of the observation and spectrometer.
#        """
#        # Ensure that a proper coordinate system has been supplied
#        if coord_sys != 'ga' and coord_sys != 'eq':
#            raise ValueError("Invalid coordinate system supplied: " + coord_sys)
#
#        # Set times
#        obs_start_unix = time.time() #unix time
#        unix_object = Time(obs_start_unix, format='unix', 
#                           location=(self.LON, self.LAT, self.ALT)) #unix time Time object
#        obs_start_jd = unix_object.jd #convert unix time to julian date
#
#        # Set the coordinates
#        if coord_sys == 'ga':
#            l, b = coords*u.degree
#            c = SkyCoord(l=l, b=b, frame='galactic')
#            equatorial = c.fk5
#            ra, dec = equatorial.ra, equatorial.dec
#        elif coord_sys == 'eq':
#            ra, dec = coords*u.degree
#            c = SkyCoord(ra, dec)
#            galactic = c.galactic
#            l, b = galactic.l, galactic.b
#        
#        # Make PrimaryHDU
#        header = fits.Header()
#
#        # Save metadata of the system and spectrometer
#        header['NSPEC'] = (nspec, 'Number of spectra collected')
#        header['FPGFILE'] = (self.fpgfile, 'FPGA FPG file')
#        header['SNAP'] = (self.snap, 'SNAP board')
#        header['HOST'] = (self.s.fpga.host, 'Host of the FPGA')
#        header['ACCLEN'] = (self.acc_len, 'Number of clock cycles')
#        header['SPEC/ACC'] = (self.spec_per_acc, 'Spectra per accumulation')
#        header['STREAM_1'] = (self.stream_1, 'First ADC port used')
#        header['STREAM_2'] = (self.stream_2, 'Second ADC port used')
#
#        header['PYTHON'] = (3.8, 'Python version')
#        header['SRC'] = ('https://github.com/darbymccauley/Leuschner_Spectrometer.git', 'Source code')
#        # header['CASPERFPGA'] = (CASPERFPGA_VERSION, 'casperfpga code used')
#        # header['HERA_CORR_F'] = (HERA_CORR_F_VERSION, 'hera_corr_f code used')
#        
#        # Save observation attributes
#        header['L'] = (l.value, 'Galactic longitude [deg]')
#        header['B'] = (b.value, 'Galactic latitude [deg]')
#        header['RA'] = (ra.value, 'Right Ascension [deg]')
#        header['DEC'] = (dec.value, 'Declination [deg]')
#        header['JD'] = (obs_start_jd, 'Julian date of start time')
#        header['UNIX'] = (obs_start_unix, 'Seconds since epoch')
#
#        primaryhdu = fits.PrimaryHDU(header=header)
#        return primaryhdu
#
#
#    def _wait_for_cnt(self):
#        """
#        Waits for corr_0 acc_cnt to increase by 1. 
#        Returns the count read from the register.
#        (Sourced with modifications from hera_corr_f.)
#        """
#        cnt_0 = self.corr_0.read_uint('acc_cnt')
#        while self.corr_0.read_uint('acc_cnt') < cnt_0 + 1:
#            time.sleep(DELAY_TIME)
#        return cnt_0
#
#
#    def _get_new_corr(self, corr, pol1, pol2):
#        """
#        Reads and returns the spectra collected given a set of polarizations.
#
#        Inputs:
#        - corr: which correlator to use
#        - pol1: first polarization
#        - pol2: second polarization
#
#        Returns: correlated data, either auto or cross depending on choice of pol1 and pol2.
#
#        (Sourced with modifications from hera_corr_f.)
#        """
#        corr.set_input(pol1, pol2)
#        spec = corr.read_bram(flush_vacc=False) / float(self.acc_len * self.spec_per_acc)
#        if pol1 == pol2:
#            return spec.real + 1j * np.zeros(len(spec))
#        else:
#            return spec
#
#    
#    def _progress_bar(self, nspec, progress=False):
#        """
#        Add a progress bar to keep track of spectra collection.
#        
#        Inputs:
#        - nspec: Number of spectra to collect.
#        - progress: Add progress bar.
#        """
#        if progress == False:
#            return range(nspec)
#        elif progress == True:
#            return tqdm.tqdm(range(nspec), desc='Progress')
#
#
#    def read_spec(self, filename, nspec, coords, coord_sys='ga', progress=False):
#        """
#        Recieves spectrometer data from the Leuschner spectrometer and 
#        saves it to a FITS file. The primary HDU contains information about
#        the observation (coordinates, number of spectra collected, time,
#        etc.) and spectrometer attributes used. Each set of spectra is
#        stored in its own FITS table in the FITS file. The columns in
#        each FITS table are ''auto0_real'' and ''auto1_real'',
#        for each polarization's auto-correlation. All columns contain
#        double-precision floating-point numbers.
#
#        Inputs:
#        - filename: Name of the output FITs file.
#        - nspec: Number of spectra to collect.
#        - coords: Coordinate(s) of the target.
#            Format: (l/ra, b/dec)
#        - coord_sys: Coordinate system used for ''coords''.
#            Default is galactic coordinates. Takes in either galactic 
#            ('ga') or equatorial ('eq') coordinate systems.
#        - progress: shows progress bar for data collection.
#        Returns:
#        - FITS file with autocorrelated spectrometer data.
#        """
#        # Make PrimaryHDU for FITS file
#        primaryhdu = self.make_PrimaryHDU(nspec, coords, coord_sys)
#        hdulist = fits.HDUList(hdus=[primaryhdu])
#
#        # Define spectra to collect
#        spectra = [('auto0_real', self.corr_0, (self.stream_1, self.stream_1)), # (0, 0)
#                   ('auto1_real', self.corr_1, (self.stream_2, self.stream_2))] # (1, 1)
#        data = {}
#        
#        # Collect spectra
#        logging.info('Reading %s spectra from the SNAP.' % str(nspec))
#        for ninteg in self._progress_bar(nspec, progress):
#            cnt_0 = self._wait_for_cnt()
#            for name, corr, (stream_1, stream_2) in spectra: # read the spectra from both corrs
#                data[name] = self._get_new_corr(corr, stream_1, stream_2).real
#            cnt_1 = self.corr_1.read_uint('acc_cnt')
#            assert cnt_0 + 1 == cnt_1 # assert corr_0's count increased and matches corr_1's count
#
#            # Make BinTableHDU and append collected data
#            data_list = [fits.Column(name=name, format='D', array=data[name]) for name, _, _ in spectra]
#            bintablehdu = fits.BinTableHDU.from_columns(data_list, name='CORR_DATA')
#            hdulist.append(bintablehdu)
#
#        # Save the output file
#        hdulist.writeto(filename, overwrite=True)
#        hdulist.close()
#
#
#    def read_corr(self, filename, nspec, coords, coord_sys='ga', progress=False):
#        """
#        Recieves correlation data from the Leuschner spectrometer and 
#        saves it to a FITS file. The primary HDU contains information about
#        the observation (coordinates, number of spectra collected, time,
#        etc.) and spectrometer attributes used. Each set of spectra is
#        stored in its own FITS table in the FITS file. The columns in
#        each FITS table are ''cross_real'' and ''cross_imag'', for the real
#        and imaginary components of the cross-correlated data. All columns 
#        contain double-precision floating-point numbers.
#
#        Inputs:
#        - filename: Name of the output FITs file.
#        - nspec: Number of spectra to collect.
#        - coords: Coordinate(s) of the target.
#            Format: (l/ra, b/dec)
#        - coord_sys: Coordinate system used for ''coords''.
#            Default is galactic coordinates. Takes in either galactic 
#            ('ga') or equatorial ('eq') coordinate systems.
#        - progress: shows progress bar for data collection.
#        
#        Returns:
#        - FITS file with correlated spectrometer data.
#        """
#        primaryhdu = self.make_PrimaryHDU(nspec, coords, coord_sys)
#        hdulist = fits.HDUList(hdus=[primaryhdu])
#
#        logging.info('Reading %s spectra from the SNAP' % str(nspec))
#        spectra = [('cross', (self.stream_1, self.stream_2))] # (0, 1)
#
#        # Collect spectra
#        for ninteg in self._progress_bar(nspec, progress):
#            data_list = []
#            for name, (stream_1, stream_2) in spectra:
#                cross = self.corr_0.get_new_corr(stream_1, stream_2)
#                cross_real, cross_imag = cross.real, cross.imag
#            data_list.append(fits.Column(name=name+'_real', format='D', array=cross_real))
#            data_list.append(fits.Column(name=name+'_imag', format='D', array=cross_imag))
#
#            bintablehdu = fits.BinTableHDU.from_columns(data_list, name='CORR_DATA')
#            hdulist.append(bintablehdu)
# 
#        # Save the output file
#        hdulist.writeto(filename, overwrite=True)
#        hdulist.close()
