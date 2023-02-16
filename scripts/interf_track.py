#! /usr/bin/env python
'''
Track a source with the interferometer. Can plot data on the fly.
'''

import numpy as np
import optparse
import time, sys
import ugradio
from ugradio.hp_multi import HP_Multimeter
from ugradio.interf import Interferometer
from ugradio.interf import AZ_MIN, AZ_MAX
#from ugradio.interf_delay import DelayClient
from astropy.utils.iers import conf
conf.auto_max_age = None

import astropy.coordinates, astropy.time

# Parse command-line arguments
o = optparse.OptionParser()
o.set_usage('interf_track.py [options]')
o.set_description(__doc__)
#o.add_option('-d', '--delay', dest='delay', action='store_true', default=False,
#    help='Use the delay lines to remove geometric delay (EXPERIMENTAL)')
o.add_option('-s', '--src', dest='src', default='Sun',
    help='Named source to track.')
o.add_option('--no_rec', dest='no_rec', action='store_true', default=False,
    help='Turn off data recording.')
o.add_option('--plot', dest='plot', action='store_true', default=False,
    help='Plot data in updating chart.')

opts, args = o.parse_args(sys.argv[1:])

t = astropy.time.Time(time.time(), format='unix')

# Get object to track
if opts.src == 'Sun':
    src = astropy.coordinates.get_sun(t)
elif opts.src == 'crab':
    src = astropy.coordinates.SkyCoord('05h34m32.0s', '+22d00m48s', equinox='J2000', frame='icrs') # from simbad
else:
    raise ValueError('Unrecognized source: %s' % src)

# Observing from New Campbell Hall
obs = astropy.coordinates.EarthLocation(
            lon=ugradio.nch.lon,
            lat=ugradio.nch.lat,
            height=ugradio.nch.alt,
)

#DELAY_LINE = opts.delay
TRACKING_DT = 10 # s
MAX_PLOT = 1024
PLOT = opts.plot
#FQ = 10.5 # GHz

def flip_point(alt, az):
    '''For interferometer, use over-the-top pointing to extend az bounds.'''
    if az < AZ_MIN:
        az = az + 180
        alt = 180 - alt
    elif az > AZ_MAX:
        az = az - 180
        alt = 180 - alt
    return alt, az


# XXX need replacement functionality from astropy
#freqs = np.array([FQ])
#aa = aipy.cal.get_aa(CALFILE, freqs)
#src = aipy.cal.get_catalog(CALFILE, [SRC])[SRC]
#bl = aa.get_baseline(0, 1, src='z')

print('LOCATION:', obs.lat, obs.lon)
print('Current JD:', t.jd)
print(opts.src, '(RA, DEC)', src.ra, src.dec)

interf = Interferometer()

if opts.no_rec:
    hpm = None
else:
    hpm = HP_Multimeter()

if PLOT:
    import matplotlib.pyplot as plt
    plt.ion()
    data = np.zeros(MAX_PLOT)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot = plt.plot(data)[0]

#dc = DelayClient()
cnt = 0

try:
    while True:
        t = astropy.time.Time(time.time(), format='unix')
        altaz = astropy.coordinates.AltAz(obstime=t, location=obs)
        pnt = src.transform_to(altaz)
        print('Track Iter: %d at JD %f' % (cnt, t.jd))

        # XXX need replacement functionality from astropy
        #src_pos = src.get_crds('top', ncrd=3)
        #bl_proj = np.dot(src_pos, bl) # units of ns
        #print('Baseline:', aa.get_baseline(0, 1, src=src))
        #print('Source delay [ns]:', bl_proj)
#if DELAY_LINE:
        #    dc.delay_ns(bl_proj)

        print('   ', opts.src, '(AZ, ALT)', pnt.az, pnt.alt)
        alt, az = flip_point(pnt.alt.deg, pnt.az.deg)
        if cnt == 0:
            print('    Slewing...')
            interf.point(alt, az, wait=True)
            print('    Pointing locked at', interf.get_pointing())
            if hpm:
                hpm.start_recording(dt=1)
        else:
            print('    Pointing (AZ, ALT)', az, alt)
            interf.point(alt, az, wait=False)
            if hpm:
                print('    REC STATUS:', hpm.get_recording_status())
                print('    Voltage:', hpm._volts[-1], 'at', hpm._times[-1])
                if PLOT:
                    volts, _ = hpm.get_recording_data()
                    volts = volts[-MAX_PLOT:]
                    plot.set_ydata(volts)
                    plot.set_xdata(np.arange(len(volts)))
                    plt.xlim(0, len(volts))
                    plt.ylim(min(volts), max(volts))
                    fig.canvas.draw()
                    fig.canvas.flush_events()
                if cnt % 10 == 0:
                    volts, times = hpm.get_recording_data()
                    print('    Backing up to backup.npz')
                    np.savez('backup.npz', volts=volts, times=times)
        cnt += 1
        time.sleep(TRACKING_DT)

except(KeyboardInterrupt):
    if hpm:
        hpm.end_recording()
    print('Stowing')
    #interf.stow()

import IPython; IPython.embed()
