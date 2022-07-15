import numpy as np
import sys 
import os
cwd = os.getcwd()
sys.path.append(os.path.abspath(cwd+'/src'))
from leuschner import Spectrometer

# Open IPython to continue session as desired
import IPython
Q = input("\n\nDo you wish to continue in IPython? [y/n] ")
if Q == "y":
    print("\nOpening IPython...:")
    IPython.embed()
elif Q == "n":
    print("Goodbye.")
