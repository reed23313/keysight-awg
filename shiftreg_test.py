# This program has been adapted from Keysight's SD1 docs
# It allows for synchronous capture of a device's response given an arbitrary stimulus
# It is designed to test a superconducting shift register made with nanocryotrons

# test pattern: 256 bytes (0x00, 0x01, ... 0xff)
# send pattern 489 times to 1Mbit of data (BER down to 1e-6) ---> this should be reasonable
# send pattern 488282 times to get 1Gbit of data (BER down to 1e-9) ---> only do this if lower BER is really needed

# sweep Vinput, Vclk_in, Vclk_shift, Vclk_readout

import sys
import random
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
PLOT_RESULTS = True
RAND_SEED = None # if None, use current system time as random seed
random.seed(RAND_SEED)

# AWG constants
AWG_FSAMP = int(1e9) # 1GS/s
AWG_CHANNELS = [1, 2, 3, 4]
AWG_DELAYS = [0, 0, 3, 0] # delay in ns
AWG_AMPLITUDE = [0.7, 0.7, 0.7, 0.7] # full scale in V

# DAQ constants
DAQ_FSAMP = int(500e6) # 500MS/s
DAQ_CHANNELS = [1, 2, 3, 4]
DAQ_TRIG_DELAY = 184 # set the delay until capturing samples from when the trigger condition is met
DAQ_FULL_SCALE = [1, 1, 1, 1] # full scale in V
DAQ_IMPEDANCE = [imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50]
DAQ_COUPLING = [cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC]

def make_word(value, bits, pulse, freq, fs):
    period_samples = int(fs/freq)
    if len(pulse) > period_samples:
        raise ValueError(f"frequency {freq} is too high for the selected pulse length of {len(pulse)} at sample rate {fs}")
    word = []
    for i in range(bits):
        if value & 1:
            # binary 1
            word += pulse + [0]*(period_samples - len(pulse))
        else:
            # binary 0
            word += [0]*period_samples
        value >>= 1
    return word

TEST_CYCLES = 1
WORD_SIZE = 8
NUM_WORDS = 2**WORD_SIZE
BIT_RATE = int(100e6) # sym/s (maybe up to 150MHz; 100MHz looks good however)
AWG_BUFFER_SIZE = NUM_WORDS*WORD_SIZE*(AWG_FSAMP//BIT_RATE)
DAQ_BUFFER_SIZE = NUM_WORDS*WORD_SIZE*(DAQ_FSAMP//BIT_RATE)
input_signal = []
for i in range(NUM_WORDS):
    input_signal += make_word(random.randint(0, NUM_WORDS-1), WORD_SIZE, AWG_SHORT_PULSES[0], BIT_RATE, AWG_FSAMP)
clock_word = make_word(NUM_WORDS-1, WORD_SIZE, AWG_SHORT_PULSES[0], BIT_RATE, AWG_FSAMP)
clock_signal = clock_word*NUM_WORDS

print(AWG_BUFFER_SIZE)
print(len(input_signal))
print(len(clock_signal))
print(len(clock_word))

assert len(input_signal) == AWG_BUFFER_SIZE
assert len(clock_signal) == AWG_BUFFER_SIZE

# add AWG and DAQ channels, set AWG buffer contents
awg = pxi_modules.AWG("M3202A", 1, 5, AWG_BUFFER_SIZE)
awg.add_channels(AWG_CHANNELS)
for n,c in enumerate(AWG_CHANNELS):
    awg.set_channel_delay(c, AWG_DELAYS[n])
    awg.set_channel_amplitude(c, AWG_AMPLITUDE[n])
    if c == 1:
        awg.set_buffer_contents(c, input_signal)
    else:
        awg.set_buffer_contents(c, clock_signal)

daq = pxi_modules.DAQ("M3102A", 1, 7, DAQ_BUFFER_SIZE, TEST_CYCLES)
daq.add_channels(DAQ_CHANNELS, DAQ_FULL_SCALE, DAQ_IMPEDANCE, DAQ_COUPLING)

# set up triggering modes
# by default, use EXTTRIG (or EXTTRIG_CYCLE in the case of the AWG) to allow
# for synchronization of the AWG and DAQ so no samples are lost/missed
awg.set_trigger_mode(trg.EXTTRIG_CYCLE)
daq.set_trigger_mode(trg.EXTTRIG, trigger_delay = DAQ_TRIG_DELAY)

try:
    awg.launch_channels(sum(2**(c-1) for c in AWG_CHANNELS), TEST_CYCLES)
    daq_data = daq.capture(sum(2**(c-1) for c in DAQ_CHANNELS))
except Exception as e:
    print(e)
    awg.stop()
    daq.stop()

# very important to close AWG, otherwise the buffers will not get properly flushed
print("closing AWG and DAQ")
awg.stop()
daq.stop()

if PLOT_RESULTS:
    print("plotting results")
    plt.figure()
    tvec = np.linspace(0,(DAQ_BUFFER_SIZE*TEST_CYCLES-1)/DAQ_FSAMP*1e9,DAQ_BUFFER_SIZE*TEST_CYCLES)
    colors = list(mcolors.TABLEAU_COLORS.keys())
    for n,channel in enumerate(DAQ_CHANNELS):
        # do peak finding on channel data and calculate bit error rate
        print(f'standard deviation: {np.std(daq_data[n])}')
        print(f'rms_est (assuming sinusoid): {np.quantile(daq_data[n])}')
        threshold = 3*np.std(daq_data[n])
        peaks, _ = sigproc.find_peaks(daq_data[n], height=3*np.std(daq_data[n]), distance=4)
        # if std-dev is comparable (i.e. within a factor of 2) to (Q(0.95)-Q(0.05))/sqrt(2)
        # in case pattern is very dense
        plt.plot(tvec, n + (daq_data[n]/2**15)*DAQ_FULL_SCALE[n], label = f'ch{channel}', color=colors[n])
        plt.plot(tvec[peaks], n + (daq_data[n][peaks]/2**15)*DAQ_FULL_SCALE[n], "x", color=colors[n])
    plt.legend()
    plt.xlabel("t [ns]")
    plt.ylabel("V [V]")
    plt.show()