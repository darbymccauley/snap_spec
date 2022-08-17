#! /usr/bin/env python3

import casperfpga
from hera_corr_f import SnapFengine
import argparse

from leuschner import Discover_Spectrometer, Spectrometer

parser = argparse.ArgumentParser(description='Run spectrometer.')
parser.add_argument('snap', type=str, help='which snap to configure to')
parser.add_argument('fpga', type=str, help='fpg file')
args = parser.parse_args()
SNAP = args.snap
FPGA = args.fpga

fpga = casperfpga.CasperFpga('localhost')
fpga.upload_to_ram_and_program(FPGA)
s = SnapFengine('localhost', transport='default')


if SNAP == 'Default':
    spec = Spectrometer(fpgfile=FPGA)
    spec.initialize()
elif SNAP == 'Discover':
    spec = Discover_Spectrometer(fpgfile=FPGA)
    spec.initialize()
else:
    raise IOError('Snap argument not recognized. Please specify SNAP type: Default or Discover.')

import IPython; IPython.embed()
