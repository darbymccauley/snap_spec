import casperfpga
from hera_corr_f import SnapFengine

import argparse


#--------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Collect and plot SNAP corrlator spectra.")

parser.add_argument('-pol1', action='store', dest='pol1', help='Specify a first polarization to correlate.')

parser.add_argument('-pol2', action='store', dest='pol2', help='Specify a second polarization to correlate.')

args = parser.parse_args()
print('pol1 = {0}'.format(int(args.pol1)))
print('pol2 = {0}'.format(int(args.pol2)))
#--------------------------------------------------------------------------

def take_corr_data(pol1, pol2):

    # Make fpga object and program the fpga
    fpga = casperfpga.CasperFpga('localhost')
    fpga.upload_to_ram_and_program('ugradio_corrspec_2022-02-22_0905.fpg')

    # Make SnapFengine object, program fpga, initialize adcs, and align adcs
    s = SnapFengine('localhost', transport='default')
    s.fpga.upload_to_ram_and_program('ugradio_corrspec_2022-02-22_0905.fpg')
    if s.is_programmed():
        print("\nFpga programed. \nInitializing ADCs...")
    else:
        raise IOError("\nCannot program fpga.")
    try:
        s.adc.init()
        s.align_adc()
        print("\nADCs aligned and initialized.")
    except:
        print("\nCould not properly initialize the ADCs first time around. Trying again...")
        try:
            s.adc.init()
            s.align_adc()
            print("\nADCs aligned and initialized.")
        except:
            raise IOError("\nSecond attempt at ADC initialization failed.")
            
    # Initialize blocks
    print("\nInitializing other blocks, including PFB and correlator blocks.")
    try:
        s.initialize()
    except:
        s.pfb.initialize()
        s.corr.initialize()
    print("\nBlocks initialized. \nNow collecting and plotting spectra...")

    # Collect some spectra and display it
    s.corr.plot_corr(pol1, pol2, show=True)

if __name__ == "__main__":
    take_corr_data(int(args.pol1), int(args.pol2))
