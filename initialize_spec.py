#################################################################
## This module initializes the SNAP board

## ROUGH DRAFT

#################################################################

def check_connected(timeout=10):
    """
    Checks if the SNAP is connected. Raises an IOError if the
    client cannot reach the SNAP to connect.
    
    Inputs:
    - timeout: The amount of time (in seconds) to wait before 
    raising IOError.
    """
    if ________:
        raise IOError('Cannot connect to the SNAP.')


    print('Connection to SNAP established.')



def initialize_spec(scale=False, force_restart=False):
    """
    This function starts the fpga process on the SNAP. First it 
    checks to see if the spectrometer is already running, then 
    initializes it if it isn't.

    Inputs:
    - scale: Whether or not to scale down each integration by
    the total number of spectra per integration time.
    - force_restart: Restart the fpga process even if it is 
    already running.
    """

