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

# use AWG to mimic SNSPD pulse
snspd_pulse = sio.loadmat('csvs/SPG717_sr2loop_noinput_clk1_shunt123_ntron_shiftreg_2loop_no_input_2022-05-17 17-33-06.mat')['C4y'][1,3500:3700]

CHASSIS = 1
PLOT_RESULTS = False
DEBUG_BER_CHECKER = True
DEBUG_BER_PLOT = True
DEBUG_BER = 0.01
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
NUM_WORDS = 1250
BIT_RATE = int(100e6) # sym/s (maybe up to 150MHz; 100MHz looks good however)
AWG_BUFFER_SIZE = NUM_WORDS*WORD_SIZE*(AWG_FSAMP//BIT_RATE)
DAQ_BUFFER_SIZE = NUM_WORDS*WORD_SIZE*(DAQ_FSAMP//BIT_RATE) + DAQ_NOISE_CAL_BUFFER_LENGTH
input_signal = []
error_signal = [] # generate copy of input signal with some errors based on DEBUG_BER
bits = np.zeros(NUM_WORDS*WORD_SIZE)
for i in range(NUM_WORDS):
    val = random.randint(0, NUM_WORDS-1)
    input_signal += make_word(val, WORD_SIZE, AWG_SHORT_PULSES[0], BIT_RATE, AWG_FSAMP)
    new_val = []
    for b in range(WORD_SIZE):
        if random.random() < DEBUG_BER:
            # flip each bit w.p. DEBUG_BER
            new_val.append(1 - (val & 1))
        else:
            new_val.append(val & 1)
        bits[i*WORD_SIZE + b] = val & 1
        val >>= 1
    error_signal += make_word(sum(i*2**n for n,i in enumerate(new_val)), WORD_SIZE, AWG_SHORT_PULSES[0], BIT_RATE, AWG_FSAMP)
# clocks just get 1111....
clock_word = make_word(2**WORD_SIZE-1, WORD_SIZE, AWG_SHORT_PULSES[0], BIT_RATE, AWG_FSAMP)
clock_signal = clock_word*NUM_WORDS

# make sure we've created the buffers properly
assert len(input_signal) == AWG_BUFFER_SIZE
assert len(error_signal) == AWG_BUFFER_SIZE
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
        if DEBUG_BER_CHECKER and c == 2:
            awg.set_buffer_contents(c, error_signal)
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

# calculate bit error rate
one_to_zero = np.zeros(len(DAQ_CHANNELS))
zero_to_one = np.zeros(len(DAQ_CHANNELS))
ber = np.zeros(len(DAQ_CHANNELS))
for n,channel in enumerate(DAQ_CHANNELS):
    if DEBUG_BER_CHECKER and channel != 2:
        continue
    if not DEBUG_BER_CHECKER and channel != 3:
        continue
    # to determine threshold, use first 1000 samples: they're all zero based on the delay
    # zscore reflects the SNR of the readout circuit --- 500 is suitable for loopback,
    # but most likely a lower threshold (e.g. 5sigma) is needed for a real test
    zscore = 50
    threshold = zscore*np.std(daq_data[n,:1000])
    # get cross correlation between daq data and 2x decimated awg signal to help with alignment for BER calc
    corr = sigproc.correlate(daq_data[n], input_signal[::2])
    corr /= np.max(corr)
    lags = sigproc.correlation_lags(len(daq_data[n]), len(input_signal[::2]))
    best_lag = lags[np.argmax(corr)]
    daq_data_delayed = np.zeros(DAQ_BUFFER_SIZE)
    if best_lag > 0:
        daq_data_delayed[:len(daq_data[n])-best_lag] = daq_data[n,best_lag:]
    else:
        print("WARNING: got unphysical optimal lag, BER rate estimate should be set to >0.5")
        best_lag = -best_lag
        daq_data_delayed[best_lag:] = daq_data[n,:len(daq_data[n])-best_lag]
    for bit in range(WORD_SIZE*NUM_WORDS):
        # each symbol is DAQ_FSAMP//BIT_RATE samples
        search_min = max(0, bit*(DAQ_FSAMP//BIT_RATE) - DAQ_FSAMP//(2*BIT_RATE))
        search_max = min(DAQ_BUFFER_SIZE, bit*(DAQ_FSAMP//BIT_RATE) + DAQ_FSAMP//(2*BIT_RATE))
        if bits[bit]:
            # make sure that there wasn't a 1 -> 0 error
            # check that there's a peak in the vicinity
            if max(daq_data_delayed[search_min:search_max]) < threshold:
                #print("1->0 error at t = ", tvec_daq[search_min + np.argmax(daq_data_delayed[search_min:search_max])])
                one_to_zero[n] += 1
        else:
            # make sure there wasn't a 0 -> 1 error
            # check that there aren't any peaks in the vicinity
            if max(daq_data_delayed[search_min:search_max]) > threshold:
                #print("0->1 error at t = ", tvec_daq[search_min + np.argmax(daq_data_delayed[search_min:search_max])])
                zero_to_one[n] += 1
    ber[n] = (one_to_zero[n] + zero_to_one[n])/(WORD_SIZE*NUM_WORDS)

print("num (1->0 errors) = ", one_to_zero)
print("num (0->1 errors) = ", zero_to_one)
print("BER = ", ber)

t_units = 1e-9
siprefix = {
    1e-24:'y', 1e-21:'z', 1e-18:'a', 1e-15:'f', 1e-12:'p', 1e-9:'n',
    1e-6:'u', 1e-3:'m', 1e-2:'c', 1e-1:'d', 1e3:'k', 1e6:'M', 1e9:'G',
    1e12:'T', 1e15:'P', 1e18:'E', 1e21:'Z', 1e24:'Y'
}

if DEBUG_BER_PLOT:
    print("plotting BER results")
    tvec_daq = np.linspace(0,(DAQ_BUFFER_SIZE*TEST_CYCLES-1)/DAQ_FSAMP/t_units,DAQ_BUFFER_SIZE*TEST_CYCLES)
    tvec_awg = np.linspace(0,(AWG_BUFFER_SIZE*TEST_CYCLES-1)/AWG_FSAMP/t_units,AWG_BUFFER_SIZE*TEST_CYCLES)
    colors = list(mcolors.TABLEAU_COLORS.keys())
    fig, axs = plt.subplots(3,1)
    axs[0].plot(tvec_daq, daq_data[1], label = 'daq data', color=colors[0])
    axs[0].plot(tvec_awg, input_signal, label = 'awg output', color=colors[1])
    axs[0].legend()
    axs[0].set_xlabel(f"t [{siprefix[t_units]}s]")
    axs[0].set_ylabel("V [V]")
    axs[1].plot(lags/DAQ_FSAMP/t_units, corr, label = "corr(daq_data, awg_data)", color=colors[0])
    axs[1].legend()
    axs[1].set_xlabel(f"t [{siprefix[t_units]}s]")
    axs[2].plot(tvec_daq, daq_data_delayed, label = 'daq data', color=colors[0])
    axs[2].plot(tvec_awg[::2], input_signal[::2], label = 'awg output', color=colors[1])
    axs[2].plot(np.linspace(0,(WORD_SIZE*NUM_WORDS-1)/BIT_RATE/t_units,WORD_SIZE*NUM_WORDS), 0.6*bits, '.', label = 'bits', color=colors[2])
    axs[2].legend()
    axs[2].set_xlabel(f"t [{siprefix[t_units]}s]")
    axs[2].set_ylabel("V [V]")
    plt.show()

if PLOT_RESULTS:
    print("plotting results")
    plt.figure()
    tvec = np.linspace(0,(DAQ_BUFFER_SIZE*TEST_CYCLES-1)/DAQ_FSAMP/t_units,DAQ_BUFFER_SIZE*TEST_CYCLES)
    tvec = tvec - DAQ_NOISE_CAL_BUFFER_LENGTH/DAQ_FSAMP/t_units
    colors = list(mcolors.TABLEAU_COLORS.keys())
    for n,channel in enumerate(DAQ_CHANNELS):
        plt.plot(tvec, n + daq_data[n], label = f'ch{channel}', color=colors[n])
    plt.legend()
    plt.xlabel(f"t [{siprefix[t_units]}s]")
    plt.ylabel("V [V]")
    plt.show()
