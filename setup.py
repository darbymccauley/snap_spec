"""
UC Berkeley Leuschner radio SNAP spectrometer control software.
"""

from distutils.core import setup
import glob

if __name__ == '__main__':
    setup(name = 'leuschner_spec',
          description = 'Leuschner radio spectrometer control software.',
          long_description = __doc__,
          author = 'Darby McCauley',
          author_email = 'darbymccauley@berkeley.edu',
          url = 'https://github.com/darbymccauley/Leuschner_Spectrometer',
          package_dir = {'':'src'},
          py_modules = ['leuschner_spec'],
          scripts = glob.glob('scipts/*.py')
          )  
