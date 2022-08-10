# This program has been adapted from Keysight's SD1 docs
# It allows for synchronous capture of a device's response given an arbitrary stimulus
# It is designed to test a superconducting shift register made with nanocryotrons

# test pattern: 500 random bytes
# send pattern 250 times to 1Mbit of data (BER down to 1e-6) ---> this should be reasonable
# send pattern 250000 times to get 1Gbit of data (BER down to 1e-9) ---> only do this if lower BER is really needed

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
AWG_DELAYS = [0, 0, 5, 10] # delay in ns
AWG_AMPLITUDE = [0.7, 0.7, 0.7, 0.7] # full scale in V

# DAQ constants
DAQ_FSAMP = int(500e6) # 500MS/s
DAQ_CHANNELS = [1, 2, 3, 4]
# set the delay until capturing samples from when the trigger condition is met
# 184 samples seems to be an intrinsic delay between when the DAQ and AWG start given the PXI bus trigger from the DAQ
# -1000 samples is used so that the peak finding algorithm can get an estimate of the noise level of each channel
DAQ_NOISE_CAL_BUFFER_LENGTH = 1000
DAQ_TRIG_DELAY = 184 - DAQ_NOISE_CAL_BUFFER_LENGTH
DAQ_FULL_SCALE = [1, 1, 1, 1] # full scale in V
DAQ_IMPEDANCE = [imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50]
DAQ_COUPLING = [cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC]

def make_word(value, bits, pulse, freq, fs):
    period_samples = int(fs)//int(freq)
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
NUM_WORDS = 500
BIT_RATE = int(100e6) # sym/s (maybe up to 150MHz; 100MHz looks good however)
AWG_BUFFER_SIZE = NUM_WORDS*WORD_SIZE*(AWG_FSAMP//BIT_RATE)
DAQ_BUFFER_SIZE = NUM_WORDS*WORD_SIZE*(DAQ_FSAMP//BIT_RATE)
input_signal = []
for i in range(NUM_WORDS):
    input_signal += make_word(random.randint(0, NUM_WORDS-1), WORD_SIZE, AWG_SHORT_PULSES[0], BIT_RATE, AWG_FSAMP)
# clocks just get 1111....
clock_word = make_word(2**WORD_SIZE-1, WORD_SIZE, AWG_SHORT_PULSES[0], BIT_RATE, AWG_FSAMP)
clock_signal = clock_word*NUM_WORDS

# make sure we've created the buffers properly
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
    # unit conversion to V
    for n in range(len(DAQ_CHANNELS)):
        daq_data[n] = (daq_data[n]/2**15)*DAQ_FULL_SCALE[n]
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
    t_units = 1e-6
    siprefix = {
        1e-24:'y', 1e-21:'z', 1e-18:'a', 1e-15:'f', 1e-12:'p', 1e-9:'n',
        1e-6:'u', 1e-3:'m', 1e-2:'c', 1e-1:'d', 1e3:'k', 1e6:'M', 1e9:'G',
        1e12:'T', 1e15:'P', 1e18:'E', 1e21:'Z', 1e24:'Y'
    }
    tvec = np.linspace(0,(DAQ_BUFFER_SIZE*TEST_CYCLES-1)/DAQ_FSAMP/t_units,DAQ_BUFFER_SIZE*TEST_CYCLES)
    tvec = tvec - DAQ_NOISE_CAL_BUFFER_LENGTH/DAQ_FSAMP/t_units
    colors = list(mcolors.TABLEAU_COLORS.keys())
    for n,channel in enumerate(DAQ_CHANNELS):
        # do peak finding on channel data and calculate bit error rate
        # use first 1000 samples, since they are going to be all zero based on the delay
        # zscore reflects the SNR of the readout circuit --- 500 is suitable for loopback,
        # but most likely a lower threshold (e.g. 5sigma) is needed for a real test
        zscore = 50
        threshold = zscore*np.std(daq_data[n,:1000])
        distance = int(0.8*(DAQ_FSAMP/BIT_RATE))
        peaks, _ = sigproc.find_peaks(daq_data[n], height=threshold, distance=distance)
        # calculate bit error rate based on peaks
        plt.plot(tvec, n + daq_data[n], label = f'ch{channel}', color=colors[n])
        plt.plot(tvec[peaks], n + daq_data[n][peaks], "x", color=colors[n])
    plt.legend()
    plt.xlabel(f"t [{siprefix[t_units]}s]")
    plt.ylabel("V [V]")
    plt.show()
