#################################################################
## A very simple receiver script for the Leuschner Radio 
## Observatory at UC Berkeley.
##
## Darby Mccauley
##
##################################################################


from astropy.io import fits
import sys
import leuschner_spec

def unpack_attr(obs_file):
"""
Unpack the observation atrributes from a FITS file.

Inputs:
- obs_file: Observation FITS file.
Returns:
- nspec: Number of spectra collected.
- coords: Observed coordinates.
- coords_sys: Coordinate system used by ''coords''.
"""
    obs_hdu = fits.open(name=obs_file)
    header = obs_hdu[0].header
    nspec = header['NSPEC']
    coords_sys = header['COORDSYS']
    
    # Extract coordinates
    if coords_sys == 'ga':
        l = header['L']
        b = header['b']
        coords = (l, b)
    elif coords_sys == 'eq':
        ra = header['RA']
        dec = header['DEC']
        coords = (ra, dec)
    else:
        raise ValueError ('Invalid coordinate system supplied: ' + coords_sys)
        sys.exit(1)
    obs_hdu.close()

    return nspec, coords, coords_sys


if __name__ == '__main__':
    usage = 'Usage: Python3 leuschner_spec_rx.py obs_attr.fits out_file.fits'
    assert len(sys.argv) == 3, usage
    obs_file = sys.argv[1]
    out_file = sys.argv[2]
    nspec, coords, system = unpack_attr(obs_file)


    spec = leuschner_spec.Spectrometer(hostname='')
    spec.check_if_connected()
    spec.initialize_spec()
    spec.read_spec(out_file, nspec, coords, coords_sys)
    sys.exit(0)

