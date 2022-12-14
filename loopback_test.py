# This program has been adapted from Keysight's SD1 docs
# It allows for synchronous capture of a device's response given an arbitrary stimulus
# channels 1 and 2 emulate ntron/SNSPD pulses
# pulse height of SNSPD waveforms is inconsistent
# this is mostly due to the low sample rate of the 500MS/s digitizer
# ----------

import sys
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pxi_modules
import scipy.io as sio
import scipy.signal as sigproc

from pulse_lib import *

sys.path.append('C:\Program Files (x86)\Keysight\SD1\Libraries\Python')
from keysightSD1 import AIN_Impedance as imp
from keysightSD1 import AIN_Coupling as cpl
from keysightSD1 import SD_TriggerModes as trg

CHASSIS = 1

# AWG constants
AWG_TSAMP = 1e-9 # 1GS/s
AWG_CHANNELS = [1, 2, 3, 4]
#AWG_DELAYS = [0, 2, 7, 12] # delay in ns
AWG_DELAYS = [0, 0, 0, 0] # delay in ns
AWG_AMPLITUDE = [0.7, 0.7, 0.7, 0.7] # full scale in V
AWG_BUFFER_LENGTH = 1000 # number of samples to hold in buffer

# DAQ constants
DAQ_TSAMP = 2e-9 # 500MS/s
DAQ_CHANNELS = [1, 2]
DAQ_POINTS_PER_CYCLE = 500 # number of samples in a single acquisition cycle/frame
DAQ_CYCLES = 1 # number of acquisition cycles/frames
DAQ_TRIG_DELAY = 180 # set the delay until capturing samples from when the trigger condition is met
DAQ_FULL_SCALE = [1, 1, 1, 1] # full scale in V
DAQ_IMPEDANCE = [imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50]
DAQ_COUPLING = [cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC]

# load .mat data to simulate SNSPD/ntron pulses
mat = sio.loadmat('csvs\SPG717_sr2loop_noinput_clk1_shunt123_ntron_shiftreg_2loop_no_input_2022-05-17 17-33-06.mat')
c3_full_bw = mat['C3y'][1,2:]/0.35 # scale to roughly unit amplitude
c4_full_bw = mat['C4y'][1,2:]/0.35 # scale to roughly unit amplitude
c3 = c3_full_bw[::10] # downsample 10x to go from 10GS/s to 1GS/s
c4 = c4_full_bw[::10] # downsample 10x to go from 10GS/s to 1GS/s

awg_data = [
    #c3,
    #c4,
    ((AWG_LONG_PULSE + [0]*29)*4 + [0]*20)*4 + [0]*200,
    ((AWG_LONG_PULSE + [0]*29)*4 + [0]*20)*4 + [0]*200,
    ((AWG_LONG_PULSE + [0]*29)*4 + [0]*20)*4 + [0]*200,
    ((AWG_LONG_PULSE + [0]*29)*4 + [0]*20)*4 + [0]*200
    #([1]*98 + [0.5, 0.1, -0.1, -0.5] + [-1]*98)*5
]

# add AWG and DAQ channels, set AWG buffer contents
awg = pxi_modules.AWG("M3202A", 1, 5, AWG_BUFFER_LENGTH)
awg.add_channels(AWG_CHANNELS)
for n,c in enumerate(AWG_CHANNELS):
    awg.set_channel_delay(c, AWG_DELAYS[n])
    awg.set_channel_amplitude(c, AWG_AMPLITUDE[n])
    awg.set_channel_offset(c, 0)
    awg.set_buffer_contents(c, awg_data[n])

daq = pxi_modules.DAQ("M3102A", 1, 7, DAQ_POINTS_PER_CYCLE, DAQ_CYCLES)
daq.add_channels(DAQ_CHANNELS, DAQ_FULL_SCALE, DAQ_IMPEDANCE, DAQ_COUPLING)

# set up triggering modes
# by default, use EXTTRIG (or EXTTRIG_CYCLE in the case of the AWG) to allow
# for synchronization of the AWG and DAQ so no samples are lost/missed
awg.set_trigger_mode(trg.EXTTRIG_CYCLE)
daq.set_trigger_mode(trg.EXTTRIG, trigger_delay = DAQ_TRIG_DELAY)
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
try:
    awg.launch_channels(sum(2**(c-1) for c in AWG_CHANNELS), 0)
    daq_data = daq.capture(sum(2**(c-1) for c in DAQ_CHANNELS))
except Exception as e:
    print(e)
    awg.stop()
    daq.stop()

input('press any key to stop awg')

# very important to close AWG, otherwise the buffers will not get properly flushed
print("closing AWG and DAQ")
awg.stop()
daq.stop()


print("plotting results")
fig, axs = plt.subplots(2,1)
tvec = np.linspace(0,(DAQ_POINTS_PER_CYCLE*DAQ_CYCLES-1)*DAQ_TSAMP*1e9,DAQ_POINTS_PER_CYCLE*DAQ_CYCLES)
colors = list(mcolors.TABLEAU_COLORS.keys())
for n,channel in enumerate(DAQ_CHANNELS):
    # do peak finding on SNSPD data
    peaks, _ = sigproc.find_peaks(daq_data[n], height=3*np.std(daq_data[n]), distance=4)
    axs[0].plot(tvec, (daq_data[n]/2**15)*DAQ_FULL_SCALE[n], label = f'ch{channel}', color=colors[n])
    axs[0].plot(tvec[peaks], (daq_data[n][peaks]/2**15)*DAQ_FULL_SCALE[n], "x", color=colors[n])
axs[0].legend()
axs[0].set_xlabel("t [ns]")
axs[0].set_ylabel("V [V]")
axs[1].plot(np.linspace(0,1e9*10*AWG_BUFFER_LENGTH*0.1*AWG_TSAMP,10*AWG_BUFFER_LENGTH), c3_full_bw, label = 'saved scope data c3')
axs[1].plot(np.linspace(0,1e9*10*AWG_BUFFER_LENGTH*0.1*AWG_TSAMP,10*AWG_BUFFER_LENGTH), c4_full_bw, label = 'saved scope data c4')
axs[1].legend()
axs[1].set_xlabel("t [ns]")
axs[1].set_ylabel("V [V]")
fig.tight_layout()
plt.show()