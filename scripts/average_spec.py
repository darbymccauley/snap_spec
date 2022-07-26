#########################################
# Simple function that takes the average
# of the collected correlator data.
#########################

import numpy as np

def average(corr_data, name, axis=0):
    """
    Computes the average of the collected data.
    
    Inputs:
       - corr_data: correlator/spectrometer data
       - name = name of the column desired from fits file
            (i.e. 'auto0_real', 'auto1_real', 'cross_real', 'cross_imag')
       - axis: axis along which to compute the average
    Returns: averaged data
    """
    stack = np.array([corr_data[i].data[name] for i in range(1, len(corr_data)-1)])
    avg = np.mean(stack, axis=axis)
    return avg

