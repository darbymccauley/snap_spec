#! /usr/bin/env python3

import casperfpga
from hera_corr_f import SnapFengine
import argparse

from leuschner import Spectrometer

parser = argparse.ArgumentParser(description='Run spectrometer.')
parser.add_argument('fpga', type=str, help='fpg file')
args = parser.parse_args()
FPGA = args.fpga

fpga = casperfpga.CasperFpga('localhost')
fpga.upload_to_ram_and_program(FPGA)
s = SnapFengine('localhost', transport='default')

spec = Spectrometer(fpgfile=FPGA)
spec.initialize()

import IPython; IPython.embed()
