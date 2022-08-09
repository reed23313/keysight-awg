# This program has been adapted from Keysight's SD1 docs
# It allows for synchronous capture of a device's response given an arbitrary stimulus
# ----------

import sys
import matplotlib.pyplot as plt
import numpy as np
import pxi_modules

sys.path.append('C:\Program Files (x86)\Keysight\SD1\Libraries\Python')
from keysightSD1 import AIN_Impedance as imp
from keysightSD1 import AIN_Coupling as cpl
from keysightSD1 import SD_TriggerModes as trg
from keysightSD1 import SD_AIN_TriggerMode as atrg
from keysightSD1 import SD_TriggerBehaviors as trgb
from keysightSD1 import SD_TriggerExternalSources as trgext

CHASSIS = 1

# AWG constants
AWG_CHANNELS = [1, 2, 3, 4]
AWG_DELAYS = [0, 2, 7, 12] # delay in ns
AWG_AMPLITUDE = [0.7, 0.7, 0.7, 0.7] # full scale in V
AWG_BUFFER_LENGTH = 500 # number of samples to hold in buffer

AWG_PULSES = [
    [0.5, 1, 0.7, 0.2], # best (low ripple pulse) we have so far: 2.3ns FWHM
    [0.75, 1, 0.3, 0.1], # this is pretty good, moderate ripple, 1.9ns FWHM
    [0.6, 1, 0.5, 0.1], # this is pretty good, low ripple, 2.1ns FWHM
    [0.7, 1, 0.4, 0.0] # also pretty good, moderate ripple
]

# DAQ constants
TSAMP = 2e-9
DAQ_CHANNELS = [1, 2, 3, 4]
DAQ_POINTS_PER_CYCLE = 500 # number of samples in a single acquisition cycle/frame
DAQ_CYCLES = 2 # number of acquisition cycles/frames
DAQ_TRIG_DELAY = 0 # used for analog trigger to set the delay from when the trigger condition is met
DAQ_FULL_SCALE = [1, 1, 1, 1] # full scale in V
DAQ_IMPEDANCE = [imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50]
DAQ_COUPLING = [cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC]

awg_data = [
    ((AWG_PULSES[0] + [0]*16)*4 + [0]*20)*5,
    ((AWG_PULSES[0] + [0]*16)*4 + [0]*20)*5,
    ((AWG_PULSES[0] + [0]*16)*4 + [0]*20)*5,
    ((AWG_PULSES[0] + [0]*16)*4 + [0]*20)*5
    #([1]*98 + [0.5, 0.1, -0.1, -0.5] + [-1]*98)*5
]

# add AWG and DAQ channels, set AWG buffer contents
awg = pxi_modules.AWG("M3202A", 1, 5, AWG_BUFFER_LENGTH)
awg.add_channels(AWG_CHANNELS)
for n,c in enumerate(AWG_CHANNELS):
    awg.set_channel_delay(c, AWG_DELAYS[n])
    awg.set_buffer_contents(c, awg_data[n], AWG_AMPLITUDE[n])

daq = pxi_modules.DAQ("M3102A", 1, 7, DAQ_POINTS_PER_CYCLE, DAQ_CYCLES)
daq.add_channels(DAQ_CHANNELS, DAQ_FULL_SCALE, DAQ_IMPEDANCE, DAQ_COUPLING)

# set up triggering modes
# by default, use EXTTRIG (or EXTTRIG_CYCLE in the case of the AWG) to allow
# for synchronization of the AWG and DAQ so no samples are lost/missed
awg.set_trigger_mode(trg.EXTTRIG_CYCLE)
daq.set_trigger_mode(trg.EXTTRIG)
# AWG autotriggering will cause the AWG to start sending samples immediately
# after it starts. this is useful if you just want to send a repeating signal
# and probe the response with an oscilloscope
# be sure to also set the number of cycles in the awg.launch_channels() command
# awg.set_trigger_mode(trg.AUTOTRIG)
# example analog trigger (like an oscilloscope trigger)
# daq.set_trigger_mode(trg.ANALOGTRIG, trigger_delay = DAQ_TRIG_DELAY, analog_trig_chan = 1, analog_trig_type = atrg.RISING_EDGE, threshold = 0.2)
# software / HVI trigger, the PC will send a command to the DAQ after each cycle to capture another one
# decent for debugging if you don't care about any synchronization of data
# (i.e. capture cycle to capture cycle or relative to the AWG)
# daq.set_trigger_mode(trg.SWHVITRIG)

# use external triggers (for both AWG and DAQ) with PXI backplane rising edge trigger
# then use digitizer to send a PXI trigger using PXItriggerWrite(pxislot, trigger_n)
# there's a roughly 300ns delay between when the DAQ starts capturing samples to when the AWG starts sending them,
# so the number of capture samples per DAQ cycle should be about 200 plus whatever number of samples is required for the experiment
# it's important to set up PXI trigger connections in the Keysight Connection Expert application
# (there doesn't seem to be any way to do this with the SD1 API)
# the setting should be persistent across boots
# since the digitizer is in slot 7 (which is part of PXI bus 2) and the AWG is in slot 5 (PXI bus 1),
# the route from PXI bus 2 to bus 1 must also be enabled

awg.launch_channels(sum(2**(c-1) for c in AWG_CHANNELS), DAQ_CYCLES)
daq_data = daq.capture(sum(2**(c-1) for c in DAQ_CHANNELS))

# very important to close AWG, otherwise the buffers will not get properly flushed
print("closing AWG and DAQ")
awg.stop()
daq.stop()

print("plotting results")
plt.figure()
for n,channel in enumerate(DAQ_CHANNELS):
    plt.plot(np.linspace(0,(DAQ_POINTS_PER_CYCLE*DAQ_CYCLES-1)*TSAMP*1e9,DAQ_POINTS_PER_CYCLE*DAQ_CYCLES), (daq_data[n]/2**15)*DAQ_FULL_SCALE[n], label = f'ch{channel}')
plt.legend()
plt.xlabel("t [ns]")
plt.ylabel("V [V]")
plt.show()