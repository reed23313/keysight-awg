# This program has been adapted from Keysight's SD1 docs
# It allows for synchronous capture of a device's response given an arbitrary stimulus
# It is designed to test a superconducting shift register made with nanocryotrons

# test pattern: 500 random bytes
# send pattern 250 times to 1Mbit of data (BER down to 1e-6) ---> this should be reasonable
# send pattern 250000 times to get 1Gbit of data (BER down to 1e-9) ---> only do this if lower BER is really needed

# sweep Vinput, Vclk_in, Vclk_shift, Vclk_readout

import sys
import traceback
import time
import datetime
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
# snspd_pulse = sio.loadmat('csvs/SPG717_sr2loop_noinput_clk1_shunt123_ntron_shiftreg_2loop_no_input_2022-05-17 17-33-06.mat')['C4y'][1,3500:3700]
# for savefiles
timestamp = ''.join('_'.join(str(datetime.datetime.now()).split('.')[0].split(' ')).split(':'))
print(timestamp)

# make plots nicer
t_units = 1e-9
siprefix = {
    1e-24:'y', 1e-21:'z', 1e-18:'a', 1e-15:'f', 1e-12:'p', 1e-9:'n',
    1e-6:'u', 1e-3:'m', 1e-2:'c', 1e-1:'d', 1e3:'k', 1e6:'M', 1e9:'G',
    1e12:'T', 1e15:'P', 1e18:'E', 1e21:'Z', 1e24:'Y'
}

# custom exception class for lag calculation
class LagException(Exception):
    pass

CHASSIS = 1
SAVE = 1

# debug options
PLOT_RESULTS = False
DEBUG_OPTIMAL_AMPLITUDES = False
SAVE_TRANSIENT = False
DEBUG_DELAY_CALC = True
DEBUG_BER_PLOT = False
DEBUG_BER = 0.01
RAND_SEED = 1000*((ord('S')*256 + ord('P'))*256 + ord('G')) + 717 # if None, use current system time as random seed
random.seed(RAND_SEED)

# determine these by setting DEBUG_OPTIMAL_AMPLITUDES = True
# Vinput, Vclkin, Vclksh, Vclkro
#AWG_OPTIMAL_AMPLITUDES = [0.28, 0.25, 0.19, 0.58] # amplitude to reset shift reg state
#AWG_OPTIMAL_AMPLITUDES = [0.28, 0.22, 0.19, 0.58] # optimal for short pulses at 100MHz
#AWG_OPTIMAL_AMPLITUDES = [0.19, 0.16, 0.35, 0.75] # optimal for short pulses at 75MHz
#AWG_OPTIMAL_AMPLITUDES = [0.20, 0.19, 0.35, 0.64] # optimal for short pulses at 50MHz
# AWG_OPTIMAL_AMPLITUDES = [0.20, 0.20, 0.35, 0.64] # optimal for short pulses at 50MHz (editing)
AWG_OPTIMAL_AMPLITUDES = [0.20, 0.20, 0.20]
AWG_OPTIMAL_BIAS = [0.30, 0.30]
#AWG_OPTIMAL_AMPLITUDES = [0.25, 0.21, 0.35, 0.64] # optimal for short pulses at 25MHz
#AWG_OPTIMAL_AMPLITUDES = [0.27, 0.22, 0.32, 0.75] # optimal for short pulses at 10MHz
# AWG_OPTIMAL_AMPLITUDES = [0.15, 0.11, 0.16, 0.65] # too big for long pulses (long pulse ~ Keysight 33600A)
# "current source" impedance:
# 1.96k, 1.96k, 1.96k, 11k
# clock amplitude sweeps
N_VINPUT = 1
N_VBIAS = 1
#N_VCLKSH = 16
#N_VCLKRO = 16
##################
# N_VCLKSH = 51
# N_VCLKRO = 51
##################

V_IN_READ_MIN = 0.26
V_IN_READ_MAX = 0.26
V_BIAS_MIN = 0.26
V_BIAS_MAX = 0.26
V_INPUT_SWEEP = np.linspace(V_IN_READ_MIN,V_IN_READ_MAX,N_VINPUT)
V_BIAS_SWEEP = np.linspace(V_BIAS_MIN,V_BIAS_MAX,N_VBIAS)

########################
AWG_TSAMP = 1e-9 # 1GS/s
AWG_CHANNELS_IN = [1, 2, 3]
AWG_CHANNELS_BIAS = [1, 2]
#AWG_DELAYS = [0, 2, 7, 12] # delay in ns

# AWG_BUFFER_LENGTH = 1000 # number of samples to hold in buffer# AWG constants

# DAQ constants
DAQ_TSAMP = 2e-9 # 500MS/s
DAQ_CHANNELS = [1]
DAQ_CYCLES = 1 # number of acquisition cycles/frames
DAQ_TRIG_DELAY = 180 # set the delay until capturing samples from when the trigger condition is met
DAQ_FULL_SCALE = [1] # full scale in V
DAQ_IMPEDANCE = [imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50]
DAQ_COUPLING = [cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC]


TEST_CYCLES = 1 # 1Msamp
NUM_WORDS = 100 # should be suitably large enough so that we get good word entropy
SYMS_PER_WORD = 100 # word size in symbols
AWG_FSAMP = int(1e9) # 1GS/s
DAQ_FSAMP = int(500e6) # 500MS/s
BIT_RATE = int(50e6) # sym/s (maybe up to 150MHz; 100MHz looks good however)
AWG_WORD_SIZE = SYMS_PER_WORD*(AWG_FSAMP//BIT_RATE) # word size in samples
DAQ_WORD_SIZE = SYMS_PER_WORD*(DAQ_FSAMP//BIT_RATE) # word size in samples
AWG_BUFFER_SIZE = NUM_WORDS*AWG_WORD_SIZE # buffer size in samples
DAQ_BUFFER_SIZE = NUM_WORDS*DAQ_WORD_SIZE
PULSE_SHAPE = AWG_SHORT_PULSES[0]
# PULSE_SHAPE = AWG_LONG_PULSE
AWG_DELAYS = [int(1e9/(2*BIT_RATE)), 0, int(1e9/(2*BIT_RATE))] # delay in ns
########################

# V_INPUT_SWEEP = np.linspace(AWG_OPTIMAL_AMPLITUDES[0], AWG_OPTIMAL_AMPLITUDES[0], N_VINPUT)
# V_BIAS_SWEEP = np.linspace(AWG_OPTIMAL_AMPLITUDES[1], AWG_OPTIMAL_AMPLITUDES[1], N_VBIAS)
#V_CLKSH_SWEEP = np.linspace(AWG_OPTIMAL_AMPLITUDES[2], AWG_OPTIMAL_AMPLITUDES[2], N_VCLKSH)
#V_CLKRO_SWEEP = np.linspace(AWG_OPTIMAL_AMPLITUDES[3], AWG_OPTIMAL_AMPLITUDES[3], N_VCLKRO)
#V_CLKSH_SWEEP = np.linspace(0.25, 0.31, N_VCLKSH) # sweeps for 100MHz
#V_CLKRO_SWEEP = np.linspace(0.72, 0.78, N_VCLKRO)
#V_CLKSH_SWEEP = np.linspace(0.20, 0.44, N_VCLKSH) # sweeps for 10MHz (narrow)
#V_CLKRO_SWEEP = np.linspace(0.60, 0.90, N_VCLKRO)
#V_CLKSH_SWEEP = np.linspace(0.10, 0.60, N_VCLKSH) # sweeps for 10MHz (wide)
#V_CLKRO_SWEEP = np.linspace(0.50, 1.20, N_VCLKRO)
##################
# V_CLKSH_SWEEP = np.linspace(0.25, 0.50, N_VCLKSH) # sweeps for 75MHz
# V_CLKRO_SWEEP = np.linspace(0.50, 0.90, N_VCLKRO)
##################
#V_CLKSH_SWEEP = np.linspace(0.25, 0.50, N_VCLKSH) # sweeps for 50MHz
#V_CLKRO_SWEEP = np.linspace(0.50, 0.90, N_VCLKRO)
#V_CLKSH_SWEEP = np.linspace(0.10, 0.60, N_VCLKSH) # sweeps for 25MHz
#V_CLKRO_SWEEP = np.linspace(0.50, 1.20, N_VCLKRO)

# generate configs
# N_CONFIGS = N_VINPUT*N_VBIAS
# TEST_CONFIGURATIONS = np.zeros((N_CONFIGS,4))
# i = 0
# for v_input in V_INPUT_SWEEP:
#     for v_clk_in in V_BIAS_SWEEP:
#         for v_clk_sh in V_CLKSH_SWEEP:
#             for v_clk_ro in V_CLKRO_SWEEP:
#                 TEST_CONFIGURATIONS[i,:] = np.array([v_input, v_clk_in, v_clk_sh, v_clk_ro])
#                 i = i + 1

N_CONFIGS = N_VINPUT*N_VBIAS
peaks1 = np.zeros((N_CONFIGS, DAQ_BUFFER_SIZE*TEST_CYCLES),dtype=np.int32)
peaks2 = np.zeros((N_CONFIGS, DAQ_BUFFER_SIZE*TEST_CYCLES),dtype=np.int32)
TEST_CONFIGURATIONS = np.zeros((N_CONFIGS,2))
i = 0
for v_in in V_INPUT_SWEEP:
    for v_bias in V_BIAS_SWEEP:
        TEST_CONFIGURATIONS[i,:] = np.array([v_in, v_bias])
        i = i + 1

BIT_RATE = int(50e6) # sym/s (maybe up to 150MHz; 100MHz looks good however)

# # AWG constants
# AWG_FSAMP = int(1e9) # 1GS/s
# AWG_CHANNELS = [1, 2, 3, 4]
# AWG_DELAYS = [0, 0, int(1e9/(2*BIT_RATE)), int(1e9/BIT_RATE)] # delay in ns

# # DAQ constants
# DAQ_FSAMP = int(500e6) # 500MS/s
# DAQ_CHANNELS = [3]
# #DAQ_CHANNELS = [1, 2, 3]
# # set the delay until capturing samples from when the trigger condition is met
# # 184 samples seems to be an intrinsic delay between when the DAQ and AWG start given the PXI bus trigger from the DAQ
# DAQ_TRIG_DELAY = 184
# #DAQ_FULL_SCALE = [0.125, 0.125, 0.125] # full scale in V
# DAQ_FULL_SCALE = [1, 1, 1] # full scale in V
# DAQ_IMPEDANCE = [imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50]
# DAQ_COUPLING = [cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC]

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

# filter parameters for bit error rate calculations
FILT_B, FILT_A = sigproc.butter(1, DAQ_FSAMP/5, fs=DAQ_FSAMP, btype='low', analog=False)

# test vector parameters
# BER calculation takes roughly 120ms per channel per 1Msamp
if DEBUG_OPTIMAL_AMPLITUDES:
    TEST_CYCLES = 1
    NUM_WORDS = 10
else:
    TEST_CYCLES = 1 # 1Msamp
    NUM_WORDS = 1250 # should be suitably large enough so that we get good word entropy
SYMS_PER_WORD = 8 # word size in symbols
AWG_WORD_SIZE = SYMS_PER_WORD*(AWG_FSAMP//BIT_RATE) # word size in samples
DAQ_WORD_SIZE = SYMS_PER_WORD*(DAQ_FSAMP//BIT_RATE) # word size in samples
AWG_BUFFER_SIZE = NUM_WORDS*AWG_WORD_SIZE # buffer size in samples
DAQ_BUFFER_SIZE = NUM_WORDS*DAQ_WORD_SIZE
PULSE_SHAPE = AWG_SHORT_PULSES[0]
#PULSE_SHAPE = AWG_LONG_PULSE
input_signal = np.zeros(AWG_BUFFER_SIZE)
clock1_signal = np.zeros(AWG_BUFFER_SIZE)
clock2_signal = np.zeros(AWG_BUFFER_SIZE)
bits = np.zeros(NUM_WORDS*SYMS_PER_WORD)
clock_word = make_word(2**SYMS_PER_WORD-1, SYMS_PER_WORD, PULSE_SHAPE, BIT_RATE, AWG_FSAMP)
# generate test vectors
for i in range(NUM_WORDS):
    val = random.randint(0, 2**(SYMS_PER_WORD-1))
    if i == NUM_WORDS - 1:
        # send 0 for the last word so next cycle doesn't have an accidental bit flip at the beginning due to a held over circulating current
        # i.e., make sure that there are an equal number of pulses for each channel
        val = 0
    input_signal[i*AWG_WORD_SIZE:(i+1)*AWG_WORD_SIZE] = make_word(val, SYMS_PER_WORD, PULSE_SHAPE, BIT_RATE, AWG_FSAMP)
    for b in range(SYMS_PER_WORD):
        bits[i*SYMS_PER_WORD + b] = val & 1
        val >>= 1
    # clocks just get 1111....
    clock1_signal[i*AWG_WORD_SIZE:(i+1)*AWG_WORD_SIZE] = clock_word
    clock2_signal[i*AWG_WORD_SIZE:(i+1)*AWG_WORD_SIZE] = clock_word

# add AWG and DAQ channels
awg_in = pxi_modules.AWG("M3202A", 1, 5, AWG_BUFFER_SIZE)
awg_bias = pxi_modules.AWG("M3201A", 1, 6, AWG_BUFFER_SIZE)
daq = pxi_modules.DAQ("M3102A", 1, 7, DAQ_BUFFER_SIZE, TEST_CYCLES)
awg_in.add_channels(AWG_CHANNELS_IN)
awg_bias.add_channels(AWG_CHANNELS_BIAS)
daq.add_channels(DAQ_CHANNELS, DAQ_FULL_SCALE, DAQ_IMPEDANCE, DAQ_COUPLING)

# set up triggering modes
# by default, use EXTTRIG (or EXTTRIG_CYCLE in the case of the AWG) to allow
# for synchronization of the AWG and DAQ so no samples are lost/missed
awg_in.set_trigger_mode(trg.EXTTRIG_CYCLE)
daq.set_trigger_mode(trg.EXTTRIG, trigger_delay = DAQ_TRIG_DELAY)

noise_stddev = np.zeros(len(DAQ_CHANNELS))
noise_maxabs = np.zeros(len(DAQ_CHANNELS))
energy_stddev = np.zeros(len(DAQ_CHANNELS))
lags = np.zeros(len(DAQ_CHANNELS), dtype=np.int32)

# get noise level of DAQ inputs and measure cable delay
try:
    daq.set_capture_cycles(1)
    # first send a bunch of zeros to check noise level on DAQ inputs
    for n,c in enumerate(AWG_CHANNELS_IN):
        # set capture sizes and cycle counts so we don't gather too much data
        # weird bugs happen when changing the DAQ capture points,
        # so just decrease the number of capture cycles
        # set delays and amplitudes to zero; we don't want to send any data
        awg_in.set_channel_delay(c, AWG_DELAYS[n])
        if True:#c == 1:
            awg_in.set_channel_amplitude(c, 0)
            awg_in.set_buffer_contents(c, np.zeros(AWG_BUFFER_SIZE))
        else:
            awg_in.set_channel_amplitude(c, AWG_OPTIMAL_AMPLITUDES[n])
            awg_in.set_buffer_contents(c, clock1_signal)
    # actually capture noise data
    awg_in.launch_channels(sum(2**(c-1) for c in AWG_CHANNELS_IN), 1)
    awg_bias.launch_channels(sum(2**(c-1) for c in AWG_CHANNELS_BIAS), 1)
    daq_data = daq.capture(sum(2**(c-1) for c in DAQ_CHANNELS))
    # unit conversion to V and get noise properties
    for n in range(len(DAQ_CHANNELS)):
        daq_data[n] = (daq_data[n]/2**15)*DAQ_FULL_SCALE[n]
        daq_lpf = sigproc.lfilter(FILT_B, FILT_A, daq_data[n])
        energy = np.trapz(np.square(np.reshape(daq_lpf, (NUM_WORDS*SYMS_PER_WORD,DAQ_FSAMP//BIT_RATE))), axis=1)
        noise_stddev[n] = np.std(daq_data[n])
        noise_maxabs[n] = np.max(np.abs(daq_data[n]))
        energy_stddev[n] = np.std(energy)
    # then send a frame of random samples to determine the lag
    # for each channel (due to lightspeed delay of coax)
    bias_zerobits = np.zeros(AWG_BUFFER_SIZE)

    for n,c in enumerate(AWG_CHANNELS_IN):
        awg_in.set_channel_delay(c, AWG_DELAYS[n])
        awg_in.set_channel_amplitude(c, AWG_OPTIMAL_AMPLITUDES[n])
        if c == 1:
            awg_in.set_buffer_contents(c, input_signal)
        if c == 2:
            awg_in.set_buffer_contents(c, clock1_signal)
        if c == 3:
            awg_in.set_buffer_contents(c, clock2_signal)

    for n,c in enumerate(AWG_CHANNELS_BIAS):
        awg_bias.set_channel_offset(c, AWG_OPTIMAL_BIAS[n])
        awg_bias.set_channel_amplitude(c, 0)
        awg_bias.set_buffer_contents(c, bias_zerobits)

    # capture lag data
    awg_in.launch_channels(sum(2**(c-1) for c in AWG_CHANNELS_IN), 1)
    daq_data = daq.capture(sum(2**(c-1) for c in DAQ_CHANNELS))
    for n in range(len(DAQ_CHANNELS)):
        # unit conversion to V
        daq_data[n] = (daq_data[n]/2**15)*DAQ_FULL_SCALE[n]
        # get cross correlation between daq data and 2x decimated awg signal to help with alignment for BER calc
        corr = sigproc.correlate(daq_data[n], input_signal[::2])
        corr_lags = sigproc.correlation_lags(len(daq_data[n]), len(input_signal[::2]))
        lags[n] = corr_lags[np.argmax(corr)]
        if lags[n] < 0:
            if not DEBUG_OPTIMAL_AMPLITUDES:
                raise LagException("WARNING: got unphysical optimal lag, BER rate estimate has been set to 1 (indicating > 0.5 BER)")
            else:
                print(f"WARNING: got unphysical optimal lag on channel {DAQ_CHANNELS[n]}, BER rate estimate has been set to 1 (indicating > 0.5 BER)")
except pxi_modules.KeysightException as e:
    print("Encountered fatal exception when commanding Keysight equipment, exiting now")
    print(e)
    awg_in.stop()
    daq.stop()
    exit()
except LagException as e:
    print("Encountered fatal exception when calculating lag")
    print(e)
    awg_in.stop()
    daq.stop()
    exit()

if SAVE_TRANSIENT:
    with open(f'csvs/shiftreg_transient_{int(BIT_RATE/1e6)}MHz_{timestamp}.npy', 'wb') as f:
        np.savez(
            f,
            timestamp=timestamp,
            bitrate=BIT_RATE,
            awg_fsamp=AWG_FSAMP,
            daq_fsamp=DAQ_FSAMP,
            v_input=AWG_OPTIMAL_AMPLITUDES[0],
            v_clkin=AWG_OPTIMAL_AMPLITUDES[1],
            v_clksh=AWG_OPTIMAL_AMPLITUDES[2],
            v_clkro=AWG_OPTIMAL_AMPLITUDES[3],
            awg_waveforms=awg_in.waveforms,
            daq_data=daq_data
            )

# plot AWG and DAQ channel data and cross correlation between DAQ channels and input_signal
if DEBUG_OPTIMAL_AMPLITUDES or DEBUG_DELAY_CALC:
    tvec_daq = np.linspace(0,(DAQ_BUFFER_SIZE-1)/DAQ_FSAMP,DAQ_BUFFER_SIZE)
    tvec_awg = np.linspace(0,(AWG_BUFFER_SIZE-1)/AWG_FSAMP,AWG_BUFFER_SIZE)
    tvec_energy = np.linspace(0,(SYMS_PER_WORD*NUM_WORDS-1)/BIT_RATE,NUM_WORDS*SYMS_PER_WORD)
    corr_fig, corr_ax = plt.subplots()
    time_fig, time_ax = plt.subplots(4,1,sharex=True)
    for n in range(len(AWG_CHANNELS_IN)):
        time_ax[0].plot(tvec_awg/t_units, (len(AWG_CHANNELS_IN) - n - 1) + awg_in.waveforms[n]*0.9, label = f"AWG {AWG_CHANNELS_IN[n]}")
    daq_offset = 0
    daq_lpf_offset = 0
    daq_lpf_energy_offset = 0
    for n in range(len(DAQ_CHANNELS)):
        daq_lpf = sigproc.lfilter(FILT_B, FILT_A, daq_data[n])
        energy = np.trapz(np.square(np.reshape(daq_lpf, (NUM_WORDS*SYMS_PER_WORD,DAQ_FSAMP//BIT_RATE))), axis=1)
        time_ax[1].plot(tvec_daq/t_units, daq_offset - max(daq_data[n]) + daq_data[n], label = f"DAQ {DAQ_CHANNELS[n]}")
        time_ax[2].plot(tvec_daq/t_units, daq_lpf_offset - max(daq_lpf) + daq_lpf, label = f"lpf(DAQ {DAQ_CHANNELS[n]})")
        time_ax[3].step(tvec_energy/t_units, daq_lpf_energy_offset - max(energy) + energy, label = f"<lpf(DAQ {DAQ_CHANNELS[n]})^2>")
        daq_offset -= 1.2*(max(daq_data[n]) - min(daq_data[n]))
        daq_lpf_offset -= 1.2*(max(daq_lpf) - min(daq_lpf))
        daq_lpf_energy_offset -= 1.2*(max(energy) - min(energy))
        corr_ax.plot(corr_lags, n + corr / np.max(corr), label = f"corr(DAQ {n+1}, AWG {DAQ_CHANNELS[n]})")
    corr_ax.legend()
    time_ax[0].legend()
    time_ax[1].legend()
    time_ax[2].legend()
    time_ax[3].legend()
    time_ax[0].set_xlabel(f"t [{siprefix[t_units]}s]")
    time_ax[1].set_xlabel(f"t [{siprefix[t_units]}s]")
    corr_ax.set_xlabel(f"samples")
    time_ax[0].set_ylabel("V [V]")
    time_ax[1].set_ylabel("V [V]")
    corr_ax.set_ylabel("corr")
    plt.show()
print("noise levels (stddev) = ", noise_stddev)
print("noise levels (maxabs) = ", noise_maxabs)
print("energy levels (stddev) = ", energy_stddev)

if DEBUG_OPTIMAL_AMPLITUDES:
    awg_in.stop()
    daq.stop()
    exit()
# reset capture cycles to original value
daq.set_capture_cycles(TEST_CYCLES)

one_to_zero = np.zeros((N_CONFIGS,len(DAQ_CHANNELS)))
zero_to_one = np.zeros((N_CONFIGS,len(DAQ_CHANNELS)))
ber = np.zeros((N_CONFIGS,len(DAQ_CHANNELS)))
# set up channel delays and waveforms
for n,c in enumerate(AWG_CHANNELS_IN):
    awg_in.set_channel_delay(c, AWG_DELAYS[n])
    if c == 1:
        awg_in.set_buffer_contents(c, input_signal)
    if c == 2:
        awg_in.set_buffer_contents(c, clock1_signal)
    if c == 3:
        awg_in.set_buffer_contents(c, clock2_signal)
try:
    t0 = time.time()
    for cfg in range(N_CONFIGS):
        amplitudes = TEST_CONFIGURATIONS[cfg]
        # set up amplitudes
        for n,c in enumerate(AWG_CHANNELS_IN):
            awg_in.set_channel_amplitude(c, amplitudes[0])
            awg_in.set_channel_offset(c, 0)
 
        for n,c in enumerate(AWG_CHANNELS_BIAS):
            awg_in.set_channel_amplitude(c, 0)
            awg_in.set_channel_offset(c, amplitudes[1])

        # capture data
        #awg.launch_channels(sum(2**(c-1) for c in AWG_CHANNELS), TEST_CYCLES)
        if cfg == 0:
            awg_in.launch_channels(sum(2**(c-1) for c in AWG_CHANNELS_IN), TEST_CYCLES)
            awg_bias.launch_channels(sum(2**(c-1) for c in AWG_CHANNELS_BIAS), TEST_CYCLES)
        else:
            awg_in.launch_channels_start_only(sum(2**(c-1) for c in AWG_CHANNELS_IN))
            awg_bias.launch_channels_start_only(sum(2**(c-1) for c in AWG_CHANNELS_BIAS))

        daq_data_delayed = np.zeros((len(DAQ_CHANNELS),DAQ_BUFFER_SIZE*TEST_CYCLES))
        # unit conversion to V
        for n in range(len(DAQ_CHANNELS)):
            daq_data[n] = (daq_data[n]/2**15)*DAQ_FULL_SCALE[n]
            daq_data_delayed[n,:DAQ_BUFFER_SIZE*TEST_CYCLES-lags[n]] = daq_data[n,lags[n]:]
        # calculate bit error rate
        #t0 = time.process_time()
        for n,channel in enumerate(DAQ_CHANNELS):
            if channel != 3:
                continue
            # zscore reflects the SNR of the readout circuit --- 200-500 is suitable for loopback,
            # but most likely a lower threshold (e.g. 20 sigma) will be needed for a real test
            zscore = 20
            # if threshold is too high, then there will be incorrectly many 1->0 errors
            # if threshold is too low, then there will be incorrectly many 0->1 errors
            daq_symbols = np.zeros(SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES)
            #threshold = zscore*noise_stddev[n]
            threshold = 0.3*np.max(daq_data[n])
            #threshold = 20*energy_stddev[n]
            #daq_lpf = sigproc.lfilter(FILT_B, FILT_A, daq_data_delayed[n])
            #energy = np.trapz(np.square(np.reshape(daq_lpf, (NUM_WORDS*SYMS_PER_WORD*TEST_CYCLES,DAQ_FSAMP//BIT_RATE))), axis=1)
            #daq_symbols = energy > threshold
            # get peaks and then bin each peak into a symbol location/time
            # allow 3/4 of a clock/symbol period of separation between peaks
            peaks, _ = sigproc.find_peaks(daq_data[n], height=threshold, distance=(7*DAQ_FSAMP)//(8*BIT_RATE))
            daq_symbols[(peaks - lags[n]) // (DAQ_FSAMP//BIT_RATE)] = 1
            long_bits = np.tile(bits, TEST_CYCLES)
            zero_to_one[cfg,n] = np.sum(long_bits < daq_symbols)
            one_to_zero[cfg,n] = np.sum(long_bits > daq_symbols)
            ber[cfg,n] = (one_to_zero[cfg,n] + zero_to_one[cfg,n])/(SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES)
        #t1 = time.process_time()
        #ber_calc_time = t1 - t0
        #print("it took ", ber_calc_time, " to calculate the bit error rate")
        if (cfg % max(1,(N_CONFIGS//100))) == 0:
            print(cfg*100/N_CONFIGS, "% done")
    t1 = time.time()
    ber_calc_time = t1 - t0
    print("it took ", ber_calc_time, "s to calculate the bit error rate")
except pxi_modules.KeysightException as e:
    print("Encountered fatal exception when commanding Keysight equipment, exiting now")
    print(e)
    awg_in.stop()
    daq.stop()
    exit()
except (Exception, KeyboardInterrupt) as e:
    print("Caught generic exception, closing AWG/DAQ now")
    print(e)
    traceback.print_exc()
    awg_in.stop()
    daq.stop()
    exit()

print("num (1->0 errors) = ", one_to_zero)
print("num (0->1 errors) = ", zero_to_one)
print("BER = ", ber)
with open(f'csvs/shiftreg_ber_{int(BIT_RATE/1e6)}MHz_{timestamp}.npy', 'wb') as f:
    np.savez(
        f,
        timestamp=timestamp,
        bitrate=BIT_RATE,
        awg_fsamp=AWG_FSAMP,
        daq_fsamp=DAQ_FSAMP,
        v_input=V_INPUT_SWEEP,
        v_bias=V_BIAS_SWEEP,
        test_bits=SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES,
        zscore=zscore,
        ber=ber,
        one_to_zero=one_to_zero,
        zero_to_one=zero_to_one
        )

if SAVE == 1:
    print("saving results")
    # save_data = {"ch1": daq_data[0], "ch2": daq_data[1], "ch3": daq_data[2], "ch4": daq_data[3], "tvec": tvec, "peaks1": peaks1, "peaks2": peaks2, "peaks3": peaks3, "peaks4": peaks4}
    save_data = {"timestamp": timestamp, "bitrate": BIT_RATE, "awg_fsamp": AWG_FSAMP, "daq_fsamp": DAQ_FSAMP, "v_input": V_INPUT_SWEEP, "v_bias": V_BIAS_SWEEP, "test_bits": SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES, "zscore": zscore, "ber": ber, "one_to_zero": one_to_zero, "zero_to_one": zero_to_one}
    mat = sio.savemat(f"memory_cell_{timestamp}.mat", save_data)

# SEI ARRIVATO QUI

fig, ax = plt.subplots()
# reshape ber into a 2D plot based on Vclksh and Vclkro
ber_shro = np.zeros((N_VINPUT, N_VBIAS))
for cfg in range(N_CONFIGS):
    ampl = TEST_CONFIGURATIONS[cfg]
    ber_shro[np.where(V_INPUT_SWEEP == ampl[3])[0][0], np.where(V_BIAS_SWEEP == ampl[2])[0][0]] = ber[cfg,DAQ_CHANNELS.index(3)]
ber_min = 1/(SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES)
im = ax.imshow(np.clip(ber_shro,ber_min,1), origin='lower', norm=mcolors.LogNorm(vmin=ber_min, vmax=1))
t = np.logspace(np.ceil(np.log10(ber_min)), 0, 5)
cb = fig.colorbar(im, ax=ax)
cb.set_label('BER')
ax.set_title(f'Bit error rate {round(SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES/1e6,3)}MS at {int(BIT_RATE/1e6)}MS/s ({timestamp}) threshold = max')
ax.set_xticks((N_VINPUT//5)*np.arange(6), np.round(V_INPUT_SWEEP[::(N_VINPUT//5)]*1e3).astype(np.int32))
ax.set_yticks((N_VBIAS//5)*np.arange(6), np.round(V_BIAS_SWEEP[::(N_VBIAS//5)]*1e3).astype(np.int32))
ax.set_xlabel("V_input [mV]")
ax.set_ylabel("V_bias [mV]")
plt.savefig(f'csvs/shiftreg_ber_{int(BIT_RATE/1e6)}MHz_{timestamp}.png', bbox_inches='tight')
plt.show()
# very important to close AWG, otherwise the buffers will not get properly flushed
print("closing AWG and DAQ")
awg_in.stop()
daq.stop()

if DEBUG_BER_PLOT:
    print("plotting BER results")
    tvec_daq = np.linspace(0,(DAQ_BUFFER_SIZE*TEST_CYCLES-1)/DAQ_FSAMP/t_units,DAQ_BUFFER_SIZE*TEST_CYCLES)
    tvec_awg = np.linspace(0,(AWG_BUFFER_SIZE*TEST_CYCLES-1)/AWG_FSAMP/t_units,AWG_BUFFER_SIZE*TEST_CYCLES)
    colors = list(mcolors.TABLEAU_COLORS.keys())
    fig, axs = plt.subplots(3,1)
    axs[0].plot(tvec_daq, daq_data[2], label = 'daq data', color=colors[0])
    axs[0].plot(tvec_awg, np.tile(input_signal,TEST_CYCLES), label = 'awg output', color=colors[1])
    axs[0].legend()
    axs[0].set_xlabel(f"t [{siprefix[t_units]}s]")
    axs[0].set_ylabel("V [V]")
    #daq_data_delayed = np.zeros(DAQ_BUFFER_SIZE*TEST_CYCLES)
    #daq_data_delayed[:len(daq_data[2])-lags[2]] = daq_data[2,lags[2]:]
    axs[1].plot(tvec_daq, daq_data_delayed[2], label = 'daq data', color=colors[0])
    axs[1].plot(tvec_awg[::2], 0.5+np.tile(input_signal[::2],TEST_CYCLES), label = 'awg output', color=colors[1])
    axs[1].plot(np.linspace(0,(SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES-1)/BIT_RATE/t_units,SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES), 0.6*np.tile(bits,TEST_CYCLES), '.', label = 'bits', color=colors[2])
    axs[1].legend()
    axs[1].set_xlabel(f"t [{siprefix[t_units]}s]")
    axs[1].set_ylabel("V [V]")
    axs[2].step(np.linspace(0,(SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES-1),SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES), np.tile(bits, TEST_CYCLES), label = 'bits', color=colors[0])
    #axs[2].step(np.linspace(0,(SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES-1),SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES), 5*energy, label = 'energy', color=colors[1])
    axs[2].step(np.linspace(0,(SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES-1),SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES), 0.5*daq_symbols, label = 'syms', color=colors[1])
    axs[2].legend()
    plt.show()

if PLOT_RESULTS:
    print("plotting results")
    plt.figure()
    tvec = np.linspace(0,(DAQ_BUFFER_SIZE*TEST_CYCLES-1)/DAQ_FSAMP/t_units,DAQ_BUFFER_SIZE*TEST_CYCLES)
    colors = list(mcolors.TABLEAU_COLORS.keys())
    for n,channel in enumerate(DAQ_CHANNELS):
        plt.plot(tvec, n + daq_data[n], label = f'ch{channel}', color=colors[n])
    plt.legend()
    plt.xlabel(f"t [{siprefix[t_units]}s]")
    plt.ylabel("V [V]")
    plt.show()
