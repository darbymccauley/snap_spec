#! /usr/bin/env python3

import casperfpga
from hera_corr_f import SnapFengine
import argparse

from leuschner import Spectrometer

parser = argparse.ArgumentParser(description='Run spectrometer.')
parser.add_argument('fpga', type=str, help='fpg file to program')
parser.add_argument('is_discover', type=bool, help='Is the discover SNAP being used?')

args = parser.parse_args()
FPGA = args.fpga
DISCOVER = args.is_discover

fpga = casperfpga.CasperFpga('localhost')
fpga.upload_to_ram_and_program(FPGA)
s = SnapFengine('localhost', transport='default')

spec = Spectrometer(is_discover=DISCOVER)
spec.initialize()

import IPython; IPython.embed()
