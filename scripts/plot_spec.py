import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.io import fits
from leuschner import Spectrometer
from average_spec import average

parser = argparse.ArgumentParser(description='Plot spectra')
parser.add_argument('filename', type=str, help='fits file name')

args = parser.parse_args()
FILENAME = args.filename

# freqs = np.linspace(0, 1024, 1024)
freqs = np.fft.fftshift(np.fft.fftfreq(2*1024, 1/500))[1024:]

data = fits.open(FILENAME)

avg0 = average(data, 'auto0_real')
avg1 = average(data, 'auto1_real')

# plt.figure()
fig, [ax0, ax1] = plt.subplots(1,2, figsize=(10,4), constrained_layout=True)
plt.suptitle('Averaged auto-correlation spectra', fontsize=15)
ax0.semilogy(freqs, avg1, 'k')
ax1.semilogy(freqs, avg0, 'k')
ax0.set_title('0th polarization')
ax1.set_title('1st polarization')
ax0.set_xlabel('Frequency [MHz]')
ax1.set_xlabel('Frequency [MHz]')
ax0.set_ylabel('Power [arb.]')
plt.show()

