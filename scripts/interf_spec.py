import leuschner
import sys
import time
import numpy as np
from leuschner import Spectrometer

spec = Spectrometer(is_discover=True)
spec.startup()
#spec.initialize()
#try:
#    spec.initialize_adc()
#except(AttributeError):
#    pass
#spec.s.synth.initialize(verify=False)
#spec.s.adc.init()
#spec.s.initialize()
##spec.s.initialize_adc()
##spec.s.sync.wait_for_sync()
#spec.s.sync.set_delay(0)
#spec.s.sync.arm_sync()
#spec.initialize()

print(spec.acc_len)

filename = sys.argv[-1]
t0 = time.time()
data = []
print('Starting at', t0)
while True:
    d = spec.corr_0.get_new_corr(0,1)
    print('Read spec', time.time() - t0)
    data.append(d)
    if len(data) % 10 == 5:
        print('Saving', len(data))
        np.savez(sys.argv[-1], data=data)
