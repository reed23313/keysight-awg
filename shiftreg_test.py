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

###########################################
# AWG/DAQ constants (do not change)
###########################################
AWG_FSAMP = int(1e9) # 1GS/s
DAQ_FSAMP = int(500e6) # 500MS/s
# oversampling rate of AWG vs DAQ
# symbol period in AWG samples must be a multiple of AWG_OSR
# so that the number of DAQ samples in a symbol period is an integer
AWG_OSR = AWG_FSAMP//DAQ_FSAMP # 2
# set the delay until capturing samples from when the trigger condition is met
# 184 samples seems to be an intrinsic delay between when the DAQ and AWG start given the PXI bus trigger from the DAQ
DAQ_TRIG_DELAY = 184

###########################################
# Experiment parameters
###########################################

# debug options
DEBUG_OPTIMAL_AMPLITUDES = False
SAVE_TRANSIENT = False
DEBUG_DELAY_CALC = False
DEBUG_BER_PLOT = False
DEBUG_BER = 0.01
RAND_SEED = 1000*((ord('S')*256 + ord('P'))*256 + ord('G')) + 717 # if None, use current system time as random seed
random.seed(RAND_SEED)

# Vinput, Vclkin, Vclksh, Vclkro
# "current source" impedances:
# 1.96k, 1.96k, 1.96k, 11k
R_SHIFT = 1.96e3
R_READOUT = 11e3
MAGNET_CURRENT = 0.0 # current through magnet in A (field scales as 0.01284 T/A)
# determine these by setting DEBUG_OPTIMAL_AMPLITUDES = True
AWG_OPTIMAL_AMPLITUDES = [0.20, 0.18, 0.25, 0.90] # optimal for short pulses at 10MHz
#AWG_OPTIMAL_AMPLITUDES = [0.25, 0.20, 0.25, 0.60] # optimal for short pulses at 25MHz
#AWG_OPTIMAL_AMPLITUDES = [0.20, 0.19, 0.25, 0.60] # optimal for short pulses at 50MHz
#AWG_OPTIMAL_AMPLITUDES = [0.16, 0.12, 0.23, 0.40] # optimal for short pulses at 75MHz
#AWG_OPTIMAL_AMPLITUDES = [0.20, 0.15, 0.20, 0.60] # optimal for short pulses at 10MHz (+100mA magnet current)
#AWG_OPTIMAL_AMPLITUDES = [0.22, 0.19, 0.22, 0.60] # optimal for short pulses at 10MHz (-100mA magnet current)
#AWG_OPTIMAL_AMPLITUDES = [0.12, 0.05, 0.15, 0.50] # optimal for short pulses at 10MHz (+500mA magnet current)
#AWG_OPTIMAL_AMPLITUDES = [0.18, 0.12, 0.21, 0.60] # optimal for short pulses at 10MHz (-500mA magnet current)
# clock amplitude sweeps
#N_VCLKSH = 16
#N_VCLKRO = 16
#N_VCLKSH = 26
#N_VCLKRO = 26
N_VCLKSH = 51
N_VCLKRO = 51
#V_CLKSH_SWEEP = np.linspace(AWG_OPTIMAL_AMPLITUDES[2], AWG_OPTIMAL_AMPLITUDES[2], N_VCLKSH)
#V_CLKRO_SWEEP = np.linspace(AWG_OPTIMAL_AMPLITUDES[3], AWG_OPTIMAL_AMPLITUDES[3], N_VCLKRO)
V_CLKSH_SWEEP = np.linspace(0.1, 0.4, N_VCLKSH) # sweeps for 10MHz
V_CLKRO_SWEEP = np.linspace(0.2, 1.5, N_VCLKRO)
#V_CLKSH_SWEEP = np.linspace(0.1, 0.4, N_VCLKSH) # sweeps for 25MHz
#V_CLKRO_SWEEP = np.linspace(0.2, 1.5, N_VCLKRO)
#V_CLKSH_SWEEP = np.linspace(0.1, 0.4, N_VCLKSH) # sweeps for 50MHz
#V_CLKRO_SWEEP = np.linspace(0.2, 0.9, N_VCLKRO)
#V_CLKSH_SWEEP = np.linspace(0.1, 0.4, N_VCLKSH) # sweeps for 75MHz
#V_CLKRO_SWEEP = np.linspace(0.2, 0.5, N_VCLKRO)
#V_CLKSH_SWEEP = np.linspace(0.0, 0.4, N_VCLKSH) # sweeps for 10MHz (+100mA magnet current)
#V_CLKRO_SWEEP = np.linspace(0.2, 1.5, N_VCLKRO)
#V_CLKSH_SWEEP = np.linspace(0.0, 0.4, N_VCLKSH) # sweeps for 10MHz (-100mA magnet current)
#V_CLKRO_SWEEP = np.linspace(0.2, 1.5, N_VCLKRO)
#V_CLKSH_SWEEP = np.linspace(0.0, 0.3, N_VCLKSH) # sweeps for 10MHz (+500mA magnet current)
#V_CLKRO_SWEEP = np.linspace(0.2, 1.5, N_VCLKRO)
#V_CLKSH_SWEEP = np.linspace(0.0, 0.3, N_VCLKSH) # sweeps for 10MHz (-500mA magnet current)
#V_CLKRO_SWEEP = np.linspace(0.2, 1.5, N_VCLKRO)

# generate configs
N_CONFIGS = N_VCLKSH*N_VCLKRO
TEST_CONFIGURATIONS = np.zeros((N_CONFIGS,4))
i = 0
for v_clk_sh in V_CLKSH_SWEEP:
    for v_clk_ro in V_CLKRO_SWEEP:
        TEST_CONFIGURATIONS[i,:] = np.array([AWG_OPTIMAL_AMPLITUDES[0], AWG_OPTIMAL_AMPLITUDES[1], v_clk_sh, v_clk_ro])
        i = i + 1

BIT_RATE = int(10e6) # sym/s (maybe up to 150MHz; 100MHz looks good however)
# test vector parameters
# BER calculation takes roughly 120ms per bias point per channel per 1Msamp
if DEBUG_OPTIMAL_AMPLITUDES:
    TEST_CYCLES = 1
    NUM_WORDS = 10
else:
    # NUM_WORDS should be suitably large enough so that we get good word entropy
    # TEST_CYCLES = 1, NUM_WORDS = 1250, SYMS_PER_WORD = 8 is 10ksamp
    TEST_CYCLES = 1
    NUM_WORDS = 1250
SYMS_PER_WORD = 8 # word size in symbols
TOTAL_SYMBOLS = SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES
DAQ_SYMBOL_SIZE = DAQ_FSAMP//BIT_RATE # symbol size in samples
AWG_SYMBOL_SIZE = AWG_OSR*DAQ_SYMBOL_SIZE # symbol size in samples
AWG_WORD_SIZE = SYMS_PER_WORD*AWG_SYMBOL_SIZE # word size in samples
DAQ_WORD_SIZE = SYMS_PER_WORD*DAQ_SYMBOL_SIZE # word size in samples
AWG_BUFFER_SIZE = NUM_WORDS*AWG_WORD_SIZE # buffer size in samples
DAQ_BUFFER_SIZE = NUM_WORDS*DAQ_WORD_SIZE
TRUE_BIT_RATE = DAQ_FSAMP/DAQ_SYMBOL_SIZE # true bit rate based on quantization of symbol size in DAQ samples
PULSE_SHAPE = AWG_SHORT_PULSES[0]
#PULSE_SHAPE = AWG_LONG_PULSE

# AWG channel setup
AWG_CHANNELS = [1, 2, 3, 4]
AWG_DELAYS = [0, 0, int(round(1e9/(2*TRUE_BIT_RATE))), int(round(1e9/TRUE_BIT_RATE))] # delay in ns

# DAQ setup
#DAQ_CHANNELS = [3]
DAQ_CHANNELS = [1, 2, 3]
SHIFTREG_OUTPUT_CHANNEL = 3
DAQ_FULL_SCALE = [1, 1, 1] # full scale in V
DAQ_IMPEDANCE = [imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50]
DAQ_COUPLING = [cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC]

def make_word(value, bits, pulse, symbol_size):
    if len(pulse) > symbol_size:
        raise ValueError(f"symbol size ({symbol_size} samples) is too small for the selected pulse length of {len(pulse)}")
    word = []
    for i in range(bits):
        if value & 1:
            # binary 1
            word += pulse + [0]*(symbol_size - len(pulse))
        else:
            # binary 0
            word += [0]*symbol_size
        value >>= 1
    return word

# filter parameters for bit error rate calculations
FILT_B, FILT_A = sigproc.butter(1, DAQ_FSAMP/5, fs=DAQ_FSAMP, btype='low', analog=False)

input_signal = np.zeros(AWG_BUFFER_SIZE)
clock_signal = np.zeros(AWG_BUFFER_SIZE)
bits = np.zeros(NUM_WORDS*SYMS_PER_WORD)
clock_word = make_word(2**SYMS_PER_WORD-1, SYMS_PER_WORD, PULSE_SHAPE, AWG_SYMBOL_SIZE)
# generate test vectors
for i in range(NUM_WORDS):
    val = random.randint(0, 2**(SYMS_PER_WORD-1))
    if i == NUM_WORDS - 1:
        # send 0 for the last word so next cycle doesn't have an accidental bit flip at the beginning due to a held over circulating current
        # i.e., make sure that there are an equal number of pulses for each channel
        val = 0
    input_signal[i*AWG_WORD_SIZE:(i+1)*AWG_WORD_SIZE] = make_word(val, SYMS_PER_WORD, PULSE_SHAPE, AWG_SYMBOL_SIZE)
    for b in range(SYMS_PER_WORD):
        bits[i*SYMS_PER_WORD + b] = val & 1
        val >>= 1
    # clocks just get 1111....
    clock_signal[i*AWG_WORD_SIZE:(i+1)*AWG_WORD_SIZE] = clock_word

# add AWG and DAQ channels
awg = pxi_modules.AWG("M3202A", 1, 5, AWG_BUFFER_SIZE)
daq = pxi_modules.DAQ("M3102A", 1, 7, DAQ_BUFFER_SIZE, TEST_CYCLES)
awg.add_channels(AWG_CHANNELS)
daq.add_channels(DAQ_CHANNELS, DAQ_FULL_SCALE, DAQ_IMPEDANCE, DAQ_COUPLING)

# set up triggering modes
# by default, use EXTTRIG (or EXTTRIG_CYCLE in the case of the AWG) to allow
# for synchronization of the AWG and DAQ so no samples are lost/missed
awg.set_trigger_mode(trg.EXTTRIG_CYCLE)
daq.set_trigger_mode(trg.EXTTRIG, trigger_delay = DAQ_TRIG_DELAY)

noise_stddev = np.zeros(len(DAQ_CHANNELS))
noise_maxabs = np.zeros(len(DAQ_CHANNELS))
energy_stddev = np.zeros(len(DAQ_CHANNELS))
lags = np.zeros(len(DAQ_CHANNELS), dtype=np.int32)

# get noise level of DAQ inputs and measure cable delay
try:
    daq.set_capture_cycles(1)
    # first send a bunch of zeros to check noise level on DAQ inputs
    for n,c in enumerate(AWG_CHANNELS):
        # set capture sizes and cycle counts so we don't gather too much data
        # weird bugs happen when changing the DAQ capture points,
        # so just decrease the number of capture cycles
        # set delays and amplitudes to zero; we don't want to send any data
        awg.set_channel_delay(c, AWG_DELAYS[n])
        awg.set_channel_amplitude(c, 0)
        awg.set_buffer_contents(c, np.zeros(AWG_BUFFER_SIZE))
    # actually capture noise data
    awg.launch_channels(sum(2**(c-1) for c in AWG_CHANNELS), 1)
    daq_data = daq.capture(sum(2**(c-1) for c in DAQ_CHANNELS))
    # unit conversion to V and get noise properties
    for n in range(len(DAQ_CHANNELS)):
        daq_data[n] = (daq_data[n]/2**15)*DAQ_FULL_SCALE[n]
        daq_lpf = sigproc.lfilter(FILT_B, FILT_A, daq_data[n])
        energy = np.trapz(np.square(np.reshape(daq_lpf, (NUM_WORDS*SYMS_PER_WORD,DAQ_SYMBOL_SIZE))), axis=1)
        noise_stddev[n] = np.std(daq_data[n])
        noise_maxabs[n] = np.max(np.abs(daq_data[n]))
        energy_stddev[n] = np.std(energy)
    # then send a frame of random samples to determine the lag
    # for each channel (due to lightspeed delay of coax)
    for n,c in enumerate(AWG_CHANNELS):
        awg.set_channel_delay(c, AWG_DELAYS[n])
        awg.set_channel_amplitude(c, AWG_OPTIMAL_AMPLITUDES[n])
        if c == 1:
            awg.set_buffer_contents(c, input_signal)
        else:
            awg.set_buffer_contents(c, clock_signal)
    # capture lag data
    awg.launch_channels(sum(2**(c-1) for c in AWG_CHANNELS), 1)
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
    awg.stop()
    daq.stop()
    exit()
except LagException as e:
    print("Encountered fatal exception when calculating lag")
    print(e)
    awg.stop()
    daq.stop()
    exit()

if DEBUG_OPTIMAL_AMPLITUDES and SAVE_TRANSIENT:
    with open(f'csvs/shiftreg/shiftreg_transient_{int(TRUE_BIT_RATE/1e6)}MHz_{timestamp}.npy', 'wb') as f:
        np.savez(
            f,
            timestamp=timestamp,
            bitrate=BIT_RATE,
            true_bitrate=TRUE_BIT_RATE,
            awg_fsamp=AWG_FSAMP,
            daq_fsamp=DAQ_FSAMP,
            v_input=AWG_OPTIMAL_AMPLITUDES[0],
            v_clkin=AWG_OPTIMAL_AMPLITUDES[1],
            v_clksh=AWG_OPTIMAL_AMPLITUDES[2],
            v_clkro=AWG_OPTIMAL_AMPLITUDES[3],
            awg_waveforms=awg.waveforms,
            daq_data=daq_data,
            magnet_current=MAGNET_CURRENT
            )

# plot AWG and DAQ channel data and cross correlation between DAQ channels and input_signal
if DEBUG_OPTIMAL_AMPLITUDES or DEBUG_DELAY_CALC:
    tvec_daq = np.linspace(0,(DAQ_BUFFER_SIZE-1)/DAQ_FSAMP,DAQ_BUFFER_SIZE)
    tvec_awg = np.linspace(0,(AWG_BUFFER_SIZE-1)/AWG_FSAMP,AWG_BUFFER_SIZE)
    tvec_energy = np.linspace(0,(SYMS_PER_WORD*NUM_WORDS-1)/BIT_RATE,NUM_WORDS*SYMS_PER_WORD)
    corr_fig, corr_ax = plt.subplots()
    time_fig, time_ax = plt.subplots(2,1,sharex=True)
    for n in range(len(AWG_CHANNELS)):
        time_ax[0].plot(tvec_awg/t_units, (len(AWG_CHANNELS) - n - 1) + awg.waveforms[n]*0.9, label = f"AWG {AWG_CHANNELS[n]}")
    daq_offset = 0
    for n in range(len(DAQ_CHANNELS)):
        time_ax[1].plot(tvec_daq/t_units, daq_offset - max(daq_data[n]) + daq_data[n], label = f"DAQ {DAQ_CHANNELS[n]}")
        daq_offset -= 1.2*(max(daq_data[n]) - min(daq_data[n]))
        corr_ax.plot(corr_lags, n + corr / np.max(corr), label = f"corr(DAQ {n+1}, AWG {DAQ_CHANNELS[n]})")
    corr_ax.legend()
    time_ax[0].legend()
    time_ax[1].legend()
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
    awg.stop()
    daq.stop()
    exit()
# reset capture cycles to original value
daq.set_capture_cycles(TEST_CYCLES)

one_to_zero = np.zeros((N_CONFIGS,len(DAQ_CHANNELS)))
zero_to_one = np.zeros((N_CONFIGS,len(DAQ_CHANNELS)))
ber = np.zeros((N_CONFIGS,len(DAQ_CHANNELS)))
# set up channel delays and waveforms
for n,c in enumerate(AWG_CHANNELS):
    awg.set_channel_delay(c, AWG_DELAYS[n])
    if c == 1:
        awg.set_buffer_contents(c, input_signal)
    else:
        awg.set_buffer_contents(c, clock_signal)
try:
    t0 = time.time()
    for cfg in range(N_CONFIGS):
        amplitudes = TEST_CONFIGURATIONS[cfg]
        # set up amplitudes
        for n,c in enumerate(AWG_CHANNELS):
            awg.set_channel_amplitude(c, amplitudes[n])
        # load waveforms and start AWGs
        if cfg == 0:
            # only load waveforms into AWG at the beginning of the sweep so it runs faster
            awg.launch_channels(sum(2**(c-1) for c in AWG_CHANNELS), TEST_CYCLES)
        else:
            awg.launch_channels_start_only(sum(2**(c-1) for c in AWG_CHANNELS))
        # capture data
        daq_data = daq.capture(sum(2**(c-1) for c in DAQ_CHANNELS))
        daq_data_delayed = np.zeros((len(DAQ_CHANNELS),DAQ_BUFFER_SIZE*TEST_CYCLES))
        # unit conversion to V
        for n in range(len(DAQ_CHANNELS)):
            daq_data[n] = (daq_data[n]/2**15)*DAQ_FULL_SCALE[n]
            daq_data_delayed[n,:DAQ_BUFFER_SIZE*TEST_CYCLES-lags[n]] = daq_data[n,lags[n]:]
        # calculate bit error rate
        for n,channel in enumerate(DAQ_CHANNELS):
            if channel != SHIFTREG_OUTPUT_CHANNEL:
                # only check output channel for errors
                continue
            daq_symbols = np.zeros(TOTAL_SYMBOLS)
            # if threshold is too high, then there will be incorrectly many 1->0 errors
            # if threshold is too low, then there will be incorrectly many 0->1 errors
            threshold = 0.6*np.max(daq_data[n])
            # get peaks and then bin each peak into a symbol location/time
            # allow 3/4 of a clock/symbol period of separation between peaks
            peaks, _ = sigproc.find_peaks(daq_data[n], height=threshold, distance=(7*DAQ_FSAMP)//(8*BIT_RATE))
            daq_symbols[(peaks - lags[n]) // DAQ_SYMBOL_SIZE] = 1
            long_bits = np.tile(bits, TEST_CYCLES)
            zero_to_one[cfg,n] = np.sum(long_bits < daq_symbols)
            one_to_zero[cfg,n] = np.sum(long_bits > daq_symbols)
            ber[cfg,n] = (one_to_zero[cfg,n] + zero_to_one[cfg,n])/(TOTAL_SYMBOLS)
        if (cfg % max(1,round(N_CONFIGS/20))) == 0:
            print(round(cfg*100/N_CONFIGS), "% done")
    t1 = time.time()
    ber_calc_time = t1 - t0
    print("it took ", ber_calc_time, "s to calculate the bit error rate")

    print("num (1->0 errors) = ", one_to_zero)
    print("num (0->1 errors) = ", zero_to_one)
    print("BER = ", ber)
    with open(f'csvs/shiftreg/shiftreg_ber_{int(TRUE_BIT_RATE/1e6)}MHz_{timestamp}.npy', 'wb') as f:
        np.savez(
            f,
            timestamp=timestamp,
            bitrate=BIT_RATE,
            true_bitrate=TRUE_BIT_RATE,
            awg_fsamp=AWG_FSAMP,
            daq_fsamp=DAQ_FSAMP,
            v_input=AWG_OPTIMAL_AMPLITUDES[0],
            v_clkin=AWG_OPTIMAL_AMPLITUDES[1],
            v_clksh=V_CLKSH_SWEEP,
            v_clkro=V_CLKRO_SWEEP,
            test_bits=TOTAL_SYMBOLS,
            ber=ber,
            one_to_zero=one_to_zero,
            zero_to_one=zero_to_one,
            magnet_current=MAGNET_CURRENT
            )
    
    ber_fig, ber_ax = plt.subplots()
    z2o_o2z_fig, z2o_o2z_ax = plt.subplots(2,1)
    axs = [ber_ax, *z2o_o2z_ax]
    # reshape ber into a 2D plot based on Vclksh and Vclkro
    ber_shro = np.zeros((3, N_VCLKRO, N_VCLKSH))
    for cfg in range(N_CONFIGS):
        ampl = TEST_CONFIGURATIONS[cfg]
        ber_shro[0, np.where(V_CLKRO_SWEEP == ampl[3])[0][0], np.where(V_CLKSH_SWEEP == ampl[2])[0][0]] = ber[cfg,DAQ_CHANNELS.index(3)]
        ber_shro[1, np.where(V_CLKRO_SWEEP == ampl[3])[0][0], np.where(V_CLKSH_SWEEP == ampl[2])[0][0]] = one_to_zero[cfg,DAQ_CHANNELS.index(3)]/TOTAL_SYMBOLS
        ber_shro[2, np.where(V_CLKRO_SWEEP == ampl[3])[0][0], np.where(V_CLKSH_SWEEP == ampl[2])[0][0]] = zero_to_one[cfg,DAQ_CHANNELS.index(3)]/TOTAL_SYMBOLS
    ber_min = 1/TOTAL_SYMBOLS
    for i in range(3):
        im = axs[i].imshow(np.clip(ber_shro[i,:,:],ber_min,1), origin='lower', norm=mcolors.LogNorm(vmin=ber_min, vmax=1))
        t = np.logspace(np.ceil(np.log10(ber_min)), 0, 5)
        if i == 0:
            cb = ber_fig.colorbar(im, ax=axs[i], pad=0.15)
        else:
            cb = z2o_o2z_fig.colorbar(im, ax=axs[i], pad=0.15)
        cb.set_label('BER')
        # only use 5 ticks
        axs[i].set_xticks((N_VCLKSH//5)*np.arange(len(V_CLKSH_SWEEP[::(N_VCLKSH//5)])), np.round(V_CLKSH_SWEEP[::(N_VCLKSH//5)]*1e3).astype(np.int32))
        axs[i].set_yticks((N_VCLKRO//5)*np.arange(len(V_CLKRO_SWEEP[::(N_VCLKRO//5)])), np.round(V_CLKRO_SWEEP[::(N_VCLKRO//5)]*1e3).astype(np.int32))
        x2 = axs[i].secondary_xaxis('top', functions=(lambda v: v, lambda i: i))
        y2 = axs[i].secondary_yaxis('right', functions=(lambda v: v, lambda i: i))
        x2.set_xticks((N_VCLKSH//5)*np.arange(len(V_CLKSH_SWEEP[::(N_VCLKSH//5)])), np.round(V_CLKSH_SWEEP[::(N_VCLKSH//5)]*1e6/R_SHIFT).astype(np.int32))
        y2.set_yticks((N_VCLKRO//5)*np.arange(len(V_CLKRO_SWEEP[::(N_VCLKRO//5)])), np.round(V_CLKRO_SWEEP[::(N_VCLKRO//5)]*1e6/R_READOUT).astype(np.int32))
        x2.set_xlabel("I_clksh [uA]")
        y2.set_ylabel("I_clkro [uA]")
        axs[i].set_xlabel("V_clksh [mV]")
        axs[i].set_ylabel("V_clkro [mV]")
    total_sym_string = f'{round(TOTAL_SYMBOLS/1e6,3)}MS' if TOTAL_SYMBOLS >= 1e6 else f'{round(TOTAL_SYMBOLS/1e3,3)}kS'
    bit_rate_string = f'{round(TRUE_BIT_RATE/1e6,3)}MS/s' if TRUE_BIT_RATE>= 1e6 else f'{round(TRUE_BIT_RATE/1e3,3)}kS/s'
    axs[0].set_title(f'Bit error rate (total) {total_sym_string} at {bit_rate_string} ({timestamp})')
    axs[1].set_title(f'Bit error rate (1 -> 0) {total_sym_string} at {bit_rate_string} ({timestamp})')
    axs[2].set_title(f'Bit error rate (0 -> 1) {total_sym_string} at {bit_rate_string} ({timestamp})')
    ber_fig.tight_layout()
    z2o_o2z_fig.tight_layout()
    ber_fig.savefig(f'csvs/shiftreg/shiftreg_ber_total_{int(TRUE_BIT_RATE/1e6)}MHz_{timestamp}.png', bbox_inches='tight')
    z2o_o2z_fig.savefig(f'csvs/shiftreg/shiftreg_ber_o2z_z2o_{int(TRUE_BIT_RATE/1e6)}MHz_{timestamp}.png', bbox_inches='tight')
    plt.show()
    # very important to close AWG, otherwise the buffers will not get properly flushed
    print("closing AWG and DAQ")
    awg.stop()
    daq.stop()
except pxi_modules.KeysightException as e:
    print("Encountered fatal exception when commanding Keysight equipment, exiting now")
    print("please reboot the PC and re-open the Keysight Connection Expert after rebooting to re-enable the PXI trigger forwarding between busses 1 and 2")
    print("Exception: ", e)
    print("Stacktrace:")
    traceback.print_exc()
    awg.stop()
    daq.stop()
    exit()
except (Exception, KeyboardInterrupt) as e:
    print("Caught generic exception, closing AWG/DAQ now")
    print("Exception: ", e)
    print("Stacktrace:")
    traceback.print_exc()
    awg.stop()
    daq.stop()
    exit()

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
