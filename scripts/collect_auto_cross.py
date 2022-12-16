from leuschner import Spectrometer

# build and initialize
spec = Spectrometer(is_discover=True, location='NCH')
spec.initialize()


# Data files 
auto_file = 'interf_test_auto.fits'
cross_file = 'interf_test_cross.fits'

# Collect auto-corr
spec.read_spec(filename=auto_file, nspec=10, coords=[60,12], progress=True)

# Collect cross-corr
spec.read_corr(filename=cross_file, nspec=10, coords=[60,12], progress=True)




