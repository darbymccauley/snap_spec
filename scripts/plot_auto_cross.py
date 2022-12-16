import numpy as np
import matplotlib.pyplot as plt 
from astropy.io import fits
from average_spec import average

# Set constants
SAMPLE_RATE = 500 # [MHz]
NCHANS = 1024

# Freqs
freqs = np.fft.fftshift(np.fft.fftfreq(2*NCHANS, 1/SAMPLE_RATE))[NCHANS:]

# Data files 
auto_file = 'interf_test_auto.fits'
cross_file = 'interf_test_cross.fits'

# Read in data
auto_data = fits.open(auto_file)
cross_data = fits.open(cross_file)

# Average spectra
auto_avg0 = average(auto_data, 'auto0_real')
auto_avg1 = average(auto_data, 'auto1_real')
cross_real = average(cross_data, 'cross_real')
cross_imag = average(cross_data, 'cross_imag')

# # Plot
# fig, [(ax0, ax1), (ax2, ax3)] = plt.subplots(2, 2)
# # auto0
# ax0.semilogy(freqs, auto_avg0)
# ax0.set_title('auto0')
# # auto1
# ax1.semilogy(freqs, auto_avg1)
# ax1.set_title('auto1')
# # cross_real
# ax2.semilogy(freqs, cross_real)
# ax2.set_title('cross_real')
# # cross_imag
# ax3.semilogy(freqs, cross_imag)
# ax3.set_title('cross_imag')

plt.figure()
plt.plot(freqs, auto_avg0, label='auto0')
plt.plot(freqs, auto_avg1, label='auto1')
plt.plot(freqs, np.sqrt(np.abs(cross_real)**2 + np.abs(cross_imag)**2), label='abs') 
plt.plot(freqs, cross_real, label='cross_real')
plt.plot(freqs, cross_imag, label='cross_imag')
plt.legend()
plt.show()

