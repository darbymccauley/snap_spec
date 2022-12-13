import numpy as np
import matplotlib.pyplot as plt
import time
from leuschner import Spectrometer

spec = Spectrometer(is_discover=True)
spec.initialize()

plt.ion()
fig = plt.subplot((1,2, constrained_layout=True))
line, _ = plt.plot()
time = time.time()
while True: # infinite loop
    if time.time - time > 1: # update every 1s
        time = time.time()
        # read data
        if col_name == 'auto':
            spectra = [('auto0_real', self.s.corr_0, (self.stream_1, self.stream_1)), # (0, 0)
                   ('auto1_real', self.s.corr_1, (self.stream_2, self.stream_2))] # (1, 1)
            data = {}

        # data = #XXX
        line.set_ydata(data)
        fig.canvas.draw()
        fig.canvas.flush_events()

