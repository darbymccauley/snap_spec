#! /usr/bin/env python3

import casperfpga
from hera_corr_f import SnapFengine
import argparse

from leuschner import Spectrometer

parser = argparse.ArgumentParser(description='Run spectrometer.')
parser.add_argument('snap', types=str, help='which snap to configure to')
parser.add_argument('fpga', type=str, help='fpg file')
args = parser.parse_args()
SNAP = args.snap
FPGA = args.fpga

fpga = casperfpga.CasperFpga('localhost')
fpga.upload_to_ram_and_program(FPGA)
s = SnapFengine('localhost', transport='default')

spec = Spectrometer(fpgfile=FPGA)

if SNAP is None:
    spec.initialize()
elif SNAP is 'discover':
    spec.initialize_discover_SNAP(chipsel=[0])

import IPython; IPython.embed()
