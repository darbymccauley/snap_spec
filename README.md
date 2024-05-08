# UGRadio SNAP Spectrometer

Control code for the spectrometer used in the UC Berkeley undergraduate radio labratory (Ay121). This is a dual-polarization spectrometer with correlator functionality. This spectrometer is deployed on two of the class's three telescopes: the UC Berkeley Leuschner Radio Observatory's 3.5m dish (lab 4), and the interferometer (lab 3) located at the top of Campbell Hall. 

## Overview:
- Designed and compiled using Simulink
- 250MHz bandwidth
- 1024 frequency channels (frequency resolution of 244kHz)
- Runs on a [SNAP](https://casper.astro.berkeley.edu/wiki/SNAP) board
- Using the GPIO pins, the SNAP is connected to an RPi4
- User interaction via Python running on RPi4

Completed under the instruction and advising of Professor Aaron Parsons of UC Berkeley.

2022 Darby McCauley: darbymccauley@berkeley.edu

---

## Hardware:
The spectrometer design is programmed onto a Xilink Kintex 7 FPGA located on the SNAP, which can be reconfigured on the fly. The SNAP is also equipped with a TI LMX2581 frequency synthesizer, three HMCAD1511 ADCs, and two 10Gb Ethernet ports. We power the board using a 12V power supply. The SNAP is supplied with a 1 PPS and reference clock. The ADC is clocked at 500MHz for a total resulting bandwidth of 250MHz. 

<p align="center">
    <img align="center" src="/images/SNAP.png" >
    <img align="center" src="/images/hardware_block_diagram.png" width="50%">
</p>

(***Left:*** SNAP board. Sourced from Figure 2b [Hickish et al. 2016](https://ui.adsabs.harvard.edu/abs/2016JAI.....541001H/abstract). ***Right:*** SNAP board connections.)

---

## Spectrometer Design
<img src="/images/snap_spec_simulink.png" alt="snap_spec_simulink"/>

<img align="right" width=200 height=500 src="/images/snap_spec_design.png" alt="snap_spec_design"/>
Both polarization streams first pass through the ADC, an 8-bit sampler clocked at 500MHz. Next, the streams are sent through a  4-tap polyphase filter bank (PFB), resulting in unsigned 18b real and imaginary components. We then descramble the FFT and perform a bit selection to reduce the data to unsigned 4b real and imaginary. Lastly, the data is sent through a correlator block, which will output the desired spectra depending on the mode of the spectrometer. The output is the correlated data, with 1024 frequency channels per spectra.

---

## Software
The user can communicate with the SNAP board via a RPi4. <b>snap_spec</b> is a Python package created to interface the two and allow the client to collect data from the SNAP. <b>snap_spec</b> inherits relavant parts from the control software used for the Hydrogen Epoch of Re-ionization Array's (HERA) SNAP [F-engine](https://github.com/HERA-Team/hera_corr_f/tree/master), <b>hera_corr_f</b>. Using <b>hera_corr_f</b>, we initialize the other parts of the SNAP board, such as the ADC and the synthesizer, and set the FFT shift, equalization coefficients, and number of spectra that are integrated together.

Each ADC contains two internal correlator blocks, one for each ADC input. The user will have to specify which mode to take: 'spec' or 'corr'. 'spec' will return the auto-correlations of each stream, while 'corr' returns the cross-correlation of the two streams.

A very simple start-up example (in 'corr' mode):
```
import snap_spec

# Instantiate and initialize the spectrometer
snap = snap_spec.UGRadioSnap(host='localhost', stream_1=0, stream_2=1)
snap.initialize(mode='corr', sample_rate=500)

# Collect one cross-correlated spectra
data = snap.read_data()
```
The output is a dictionary containing the correlated spectra, fpga count, and time of the data's acquisition.
