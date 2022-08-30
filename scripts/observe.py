import casperfpga
from hera_corr_f import SnapFengine
import argparse
import numpy as np
import matplotlib.pyplot as plt
from leuschner import Spectrometer

parser = argparse.ArgumentParser(description='Run spectrometer.')
parser.add_argument('fpga', type=str, help='fpg file to program with')
parser.add_argument('is_discover', type=bool, help='is the discover SNAP being used?')
# parser.add_argument('correlation', type=str, help='auto-correlation (auto) or cross-correlation (cross)?')

args = parser.parse_args()
FPGA = args.fpga
DISCOVER = args.is_discover
# CORRELATION = args.correlation 

fpga = casperfpga.CasperFpga('localhost')
fpga.upload_to_ram_and_program(FPGA)
s = SnapFengine('localhost', transport='default')

# Build spectrometer and initialize
spec = Spectrometer(is_discover=DISCOVER)
spec.initialize()

# Observe and plot
x = np.linspace(0, 250, 1024, endpoint=False)

for i in range(2):
    data = spec.s.corr_0.get_new_corr(i,i)
    plt.plot(x, data)
plt.show()

import IPython; IPython.embed()
