"""
UC Berkeley UGRadio SNAP spectrometer control software.
"""

from distutils.core import setup
import glob

if __name__ == '__main__':
    setup(name = 'snap_spec',
        description = 'Control software for the UC Berkeley UGRadio lab SNAP spectrometer and correlator.',
        long_description = __doc__,
        author = 'Darby McCauley',
        author_email = 'darbymccauley@berkeley.edu',
        url = 'https://github.com/darbymccauley/snap_spec',
        packages = ['snap_spec'],
        package_dir = {'snap_spec':'src'},
        package_data = {'snap_spec': ['fpga/ugradio_corrspec_2022-02-22_0905.fpg']},
        scripts = glob.glob('scipts/*.py')
    )  
