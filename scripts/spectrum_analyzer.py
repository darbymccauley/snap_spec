import numpy as np
import matplotlib.pyplot as plt
import time as Time
from leuschner import Spectrometer
import argparse

parser = argparse.ArgumentParser(description='Real-time plotter of collected spectra.')
parser.add_argument('stream0', type=int, help='stream 0 (either 0 or 1)')
parser.add_argument('stream1', type=int, help='stream 1 (either 0 or 1)')

args = parser.parse_args()
STREAM0 = args.stream0
STREAM1 = args.stream1

# spec = Spectrometer(is_discover=True)
# spec.initialize()
# corr = spec.s.corr_0


class SpecAnalyzer(Spectrometer):
    def show(self):
        self.initialize()
        corr = self.s.corr_0
        freqs = np.fft.fftfreq(2*1024, 1/500)[:1024] 
        data = self._get_new_corr(corr, STREAM0, STREAM1)
        fig, ax = plt.subplots(figsize=(8,4))
        line, = plt.semilogy(freqs, data)
        plt.ylim(1e-1, 1e2)
        # plt.grid()
        plt.show()
        while True:
            data = self._get_new_corr(corr, STREAM0, STREAM1)
            line.set_ydata(data)
            fig.canvas.draw()
            fig.canvas.flush_events()
            Time.sleep(0.1)

plt.ion()
sa = SpecAnalyzer(is_discover=True)
sa.show()
