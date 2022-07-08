import casperfpga
from hera_corr_f import SnapFengine

# Make fpga object and program fpga
fpga = casperfpga.CasperFpga('localhost')
fpga.upload_to_ram_and_program('ugradio_corrspec_2022-02-22_0905.fpg')

# Make SnapFengine object, program fpga, and initialize and align ADCs
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
    print("\nCount not properly initialize the ADCs first time around. Trying again...")
    try:
        s.adc.init()
        s.align_adc()
        print("\nADCs aligned and initialized.")
    except:
        raise IOError("\nSecond attempt at ADC initialization failed.")

# Initilize other blocks
print("\nInitializing other blocks, including PFB and (both) correlator blocks...")
try:
    s.initialize()
except:
    s.pfb.initialize()
    s.corr_0.initialize()
    s.corr_1.initialize()
print("\nBlocks initialized.")

# Open IPython to continue session
import IPython
Q = input("\n\nDo you wish to continue in IPython? [y/n] ")
if Q == "y":
    print("\nOpening IPython...:")
    IPython.embed()
elif Q == "n":
    print("Goodbye.")
