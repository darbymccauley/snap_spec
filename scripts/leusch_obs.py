from leuschner import Spectrometer
import ugradio
import numpy as np
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astropy import units as u
from astropy.time import Time

import argparse

parser = argparse.ArgumentParser(description='Leuschner observing scipt.')

parser.add_argument('coord_file', type=str, help='Numpy file containing observation coordinate pairs in list format [l/ra, b/dec].')
parser.add_argument('coord_sys', type=str, help='Coordinate system of coord_file (either eq or ga).')
parser.add_argument('nspec', type=int, help='Number of spectra to collect per coordinate observation.')
parser.add_argument('filename', type=str, help='Filename in which observation data will be saved to as a FITS file. DO NOT include .fits extension.')

args = parser.parse_args()
args.coord_file = COORD_FILE
args.coord_sys = COORD_SYS
args.nspec = NSPEC
args.filename = FILENAME

filename_main, filename_noise_on, filename_noise_off = FILENAME+'_main', FILENAME+'_noise_on', FILENAME+'_noise_off'


LAT, LON, ALT = ugradio.leo.lat, ugradio.leo.lon, ugradio.leo.alt # Leuschner coords
ALT_MIN, ALT_MAX = ugradio.leusch.ALT_MIN, ugradio.leusch.ALT_MAX # min and max altitude of Leuschner
AZ_MIN, AZ_MAX = ugradio.leusch.AZ_MIN, ugradio.leusch.AZ_MAX # min and max azimuth of Leuschner

MAX_POINTING_ERRORS = 10

class LeuschObs():
    def __init__(self):
        # instantiate useful modules
        self.spec = Spectrometer(is_discover=True, location='Leuschner')
        self.telescope = ugradio.leusch.LeuschTelescope()
        self.noise = ugradio.leusch.LeuschNoise()
        self.synth = ugradio.agilent.SynthDirect()

        self.spec.initialize() # initialize spectrometer for use

        self.coordinates = np.load(COORD_FILE)


    def calc_pos_altaz(self, coords, coord_sys):
        """
        Convert input coordinates to [alt, az].

        Inputs:
            - coords (list): Input coordinates
            - coord_sys (str): Coordinate system of coords ('eq' or 'ga')
        Returns:
            - alt, az [degrees]: Altitude and azimuth
        """
        loc = EarthLocation(lat=LAT*u.deg, lon=LON*u.deg, height=ALT*u.m)
        time = Time(ugradio.timing.utc(fmt='%Y-%m-%d %X'))
        if coord_sys == 'ga':
            l, b = coords
            c = SkyCoord(l=l*u.degree, b=b*u.degree, frame='galactic')
        elif coord_sys == 'eq':
            ra, dec = coords
            c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='fk5')
        AltAz = c.transform_to(AltAz(obstime=time, location=loc))
        alt, az = AltAz.alt.degree, AltAz.az.degree
        return alt, az
        

    def observe(self):
        nfails = 0
        end_obs = False
        for coord in self.coordinates:
            alt, az = self.calc_pos_altaz(*coord) # calculate alt, az
            # check that coordinate is not outside Leuschner's pointing ranges
            if alt >= ALT_MAX or alt <= ALT_MIN:
                continue
            elif az >= AZ_MAX or az <= AZ_MIN:
                continue
            else: # if within range, point and collect data
                print('Moving to position... \n')
                try:
                    self.telescope.point(alt, az) # point telescope to alt, az
                    altaz_string = '{0:0.3f}, {1:0.3f}'.format(alt, az)
                    print('Current alt, az (degrees): {0:0.4f}, {1:0.4f}'.format(alt, az))
                    self.synth.set_frequency(635, 'MHz')
                    self.spec.read_spec(filename=filename_main+'_'+altaz_string+'.fits', nspec=NSPEC, coords=coord, coord_sys=COORD_SYS) # collect main data
                    self.synth.set_frequency(670, 'MHz')
                    self.noise.on() # turn on noise diode
                    self.spec.read_spec(filename=filename_noise_on+'_'+altaz_string+'.fits', nspec=NSPEC, coords=coord, coord_sys=COORD_SYS) # collect noise-on data
                    self.noise.off() # turn off noise diode
                    self.spec.read_spec(filename=filename_noise_off+'_'+altaz_string+'.fits', nspec=NSPEC, coords=coord, coord_sys=COORD_SYS) # collect noise-off data
                except: # if pointing fails
                    print('Failed to point to input coordinate [{0:0.3f}, {1:0.3f}] (alt, az = [{2:0.3f}, {3:0.3f}])'.format(coord[0], coord[1], alt, az), '\nPointing to next coordinate...')
                    nfails += 1 
                    if nfails == MAX_POINTING_ERRORS:
                        end_obs = True
                        break
                    continue
            if end_obs:
                print('Too many pointing errors. Ending observation.')
                break
        print('Stowing the telescope.')
        self.telescope.stow() # stow when finished


if __name__ == "__main__":
    obs = LeuschObs()
    obs.observe()