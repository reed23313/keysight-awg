# test pattern: 500 random bytes
# send pattern 250 times to 1Mbit of data (BER down to 1e-6) ---> this should be reasonable

import random
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import scipy.io as sio
import scipy.signal as sigproc
from pulse_lib import *

RAND_SEED = None # if None, use current system time as random seed
random.seed(RAND_SEED)

AWG_FSAMP = int(1e9) # 1GS/s
DAQ_FSAMP = int(500e6) # 500MS/s

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
NUM_WORDS = 25
BIT_RATE = int(100e6) # sym/s (maybe up to 150MHz; 100MHz looks good however)
AWG_BUFFER_SIZE = NUM_WORDS*WORD_SIZE*(AWG_FSAMP//BIT_RATE)
DAQ_BUFFER_SIZE = NUM_WORDS*WORD_SIZE*(DAQ_FSAMP//BIT_RATE)
input_signal = []
# daq data is artificially generated based on oscilloscope data of readout ntron pulses
# plus some jitter and added noise
daq_data_offset = 20 # artificially create some delay here
daq_data = np.zeros(DAQ_BUFFER_SIZE + daq_data_offset)

snspd_pulse = sio.loadmat('csvs/SPG717_sr2loop_noinput_clk1_shunt123_ntron_shiftreg_2loop_no_input_2022-05-17 17-33-06.mat')['C4y'][1,3500:3700]

BER = 0.01

bits = np.zeros(NUM_WORDS*WORD_SIZE)

for i in range(NUM_WORDS):
    word = random.randint(0, NUM_WORDS-1)
    input_signal += make_word(word, WORD_SIZE, AWG_SHORT_PULSES[0], BIT_RATE, AWG_FSAMP)
    # add some jitter to the timing of each SNSPD pulse
    jitter = random.randint(0, 40) # each sample is 100ps, which will correspond to 0.05 samples at DAQ
    for j in range(WORD_SIZE):
        val = word & 1
        bits[i*WORD_SIZE + j] = val
        word >>= 1
        # at a frequency of BER, flip the bit
        if random.random() < BER:
            val = 1 - val
        if not val:
            daq_data_offset += DAQ_FSAMP//BIT_RATE
            continue
        downsampled_snspd_pulse = snspd_pulse[jitter:jitter+100:20]
        daq_data[daq_data_offset:daq_data_offset+5] = downsampled_snspd_pulse
        daq_data_offset += DAQ_FSAMP//BIT_RATE
    added_noise = np.array([random.gauss(0, 0.02) for j in range(WORD_SIZE*(DAQ_FSAMP//BIT_RATE))])
    daq_data[daq_data_offset-WORD_SIZE*(DAQ_FSAMP//BIT_RATE):daq_data_offset] += added_noise

daq_data = daq_data[:DAQ_BUFFER_SIZE]

# make sure we've created the buffers properly
assert len(input_signal) == AWG_BUFFER_SIZE

print("plotting results")
t_units = 1e-9
siprefix = {
    1e-24:'y', 1e-21:'z', 1e-18:'a', 1e-15:'f', 1e-12:'p', 1e-9:'n',
    1e-6:'u', 1e-3:'m', 1e-2:'c', 1e-1:'d', 1e3:'k', 1e6:'M', 1e9:'G',
    1e12:'T', 1e15:'P', 1e18:'E', 1e21:'Z', 1e24:'Y'
}
tvec_daq = np.linspace(0,(DAQ_BUFFER_SIZE*TEST_CYCLES-1)/DAQ_FSAMP/t_units,DAQ_BUFFER_SIZE*TEST_CYCLES)
tvec_awg = np.linspace(0,(AWG_BUFFER_SIZE*TEST_CYCLES-1)/AWG_FSAMP/t_units,AWG_BUFFER_SIZE*TEST_CYCLES)
colors = list(mcolors.TABLEAU_COLORS.keys())

# find peaks in signal
threshold = 0.3*max(daq_data)
distance = int(0.8*(DAQ_FSAMP/BIT_RATE))
peaks, _ = sigproc.find_peaks(daq_data, height=threshold, distance=distance)

# get cross correlation between daq data and 2x decimated awg signal to help with alignment for BER calc
corr = sigproc.correlate(daq_data, input_signal[::2])
corr /= np.max(corr)
lags = sigproc.correlation_lags(len(daq_data), len(input_signal[::2]))
best_lag = lags[np.argmax(corr)]
daq_data_delayed = np.zeros(DAQ_BUFFER_SIZE)
if best_lag > 0:
    daq_data_delayed[:len(daq_data)-best_lag] = daq_data[best_lag:]
else:
    daq_data_delayed[best_lag:] = daq_data[:len(daq_data)-best_lag]

one_to_zero = 0
zero_to_one = 0
for bit in range(WORD_SIZE*NUM_WORDS):
    # each symbol is DAQ_FSAMP//BIT_RATE samples
    search_min = max(0, bit*(DAQ_FSAMP//BIT_RATE) - DAQ_FSAMP//(2*BIT_RATE))
    search_max = min(DAQ_BUFFER_SIZE, bit*(DAQ_FSAMP//BIT_RATE) + DAQ_FSAMP//(2*BIT_RATE))
    if bits[bit]:
        # make sure that there wasn't a 1 -> 0 error
        # check that there's a peak in the vicinity
        if max(daq_data_delayed[search_min:search_max]) < threshold:
            print("1->0 error at t = ", tvec_daq[search_min + np.argmax(daq_data_delayed[search_min:search_max])])
            one_to_zero += 1
    else:
        # make sure there wasn't a 0 -> 1 error
        # check that there aren't any peaks in the vicinity
        if max(daq_data_delayed[search_min:search_max]) > threshold:
            print("0->1 error at t = ", tvec_daq[search_min + np.argmax(daq_data_delayed[search_min:search_max])])
            zero_to_one += 1

fig, axs = plt.subplots(3,1)
axs[0].plot(tvec_daq, daq_data, label = 'daq data', color=colors[0])
axs[0].plot(tvec_daq[peaks], daq_data[peaks], "x", color=colors[0])
axs[0].plot(tvec_awg, input_signal, label = 'awg output', color=colors[1])
axs[0].legend()
axs[0].set_xlabel(f"t [{siprefix[t_units]}s]")
axs[0].set_ylabel("V [V]")
axs[1].plot(lags/DAQ_FSAMP/t_units, corr)
axs[1].set_xlabel(f"t [{siprefix[t_units]}s]")
axs[2].plot(tvec_daq, daq_data_delayed, label = 'daq data', color=colors[0])
axs[2].plot(tvec_daq, input_signal[::2], label = 'awg output', color=colors[1])
axs[2].plot(np.linspace(0,(WORD_SIZE*NUM_WORDS-1)/BIT_RATE/t_units,WORD_SIZE*NUM_WORDS), 0.6*bits, '.', label = 'bits', color=colors[2])
axs[2].legend()
axs[2].set_xlabel(f"t [{siprefix[t_units]}s]")
axs[2].set_ylabel("V [V]")
print(one_to_zero, " 1->0 errors")
print(zero_to_one, " 0->1 errors")
print("BER = ", (one_to_zero + zero_to_one)/(WORD_SIZE*NUM_WORDS))
plt.show()
