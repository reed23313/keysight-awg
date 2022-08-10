# This program has been adapted from Keysight's SD1 docs
# It allows for synchronous capture of a device's response given an arbitrary stimulus
# It is designed to test a superconducting shift register made with nanocryotrons

# test pattern: 500 random bytes
# send pattern 250 times to 1Mbit of data (BER down to 1e-6) ---> this should be reasonable
# send pattern 250000 times to get 1Gbit of data (BER down to 1e-9) ---> only do this if lower BER is really needed

# sweep Vinput, Vclk_in, Vclk_shift, Vclk_readout

import sys
import time
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
DEBUG_BER_PLOT = False
DEBUG_BER = 0.01
RAND_SEED = None # if None, use current system time as random seed
random.seed(RAND_SEED)

# will need an amplifier going into cryostat since AWG can only do +/-1.5V
# likely will want voltages in the range of 0.5-5V
# 25dB amplifier and 10dB attenuator should be good (voltage gain of 5.6)
# 20dB and 3 dB would also be appropriate (voltage gain of 7.1)
N_VINPUT = 2
N_VCLKIN = 1
N_VCLKSH = 1
N_VCLKRO = 1
V_INPUT_SWEEP = np.linspace(0.7, 1.5, N_VINPUT)
V_CLKIN_SWEEP = np.linspace(0.7, 1.5, N_VCLKIN)
V_CLKSH_SWEEP = np.linspace(0.7, 1.5, N_VCLKSH)
V_CLKRO_SWEEP = np.linspace(0.7, 1.5, N_VCLKRO)

N_CONFIGS = N_VINPUT*N_VCLKIN*N_VCLKSH*N_VCLKRO
TEST_CONFIGURATIONS = np.zeros((N_CONFIGS,4))
i = 0
for v_input in V_INPUT_SWEEP:
    for v_clk_in in V_CLKIN_SWEEP:
        for v_clk_sh in V_CLKSH_SWEEP:
            for v_clk_ro in V_CLKRO_SWEEP:
                TEST_CONFIGURATIONS[i,:] = np.array([v_input, v_clk_in, v_clk_sh, v_clk_ro])
                i = i + 1

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
DAQ_NOISE_CAL_BUFFER_LENGTH = 0#1000
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

# test parameters
TEST_CYCLES = 100
SYMS_PER_WORD = 8
NUM_WORDS = 1250 # should be suitably large enough so that we get good word entropy
BIT_RATE = int(100e6) # sym/s (maybe up to 150MHz; 100MHz looks good however)
AWG_WORD_SIZE = SYMS_PER_WORD*(AWG_FSAMP//BIT_RATE)
DAQ_WORD_SIZE = SYMS_PER_WORD*(DAQ_FSAMP//BIT_RATE)
AWG_BUFFER_SIZE = NUM_WORDS*AWG_WORD_SIZE
DAQ_BUFFER_SIZE = NUM_WORDS*DAQ_WORD_SIZE + DAQ_NOISE_CAL_BUFFER_LENGTH
input_signal = np.zeros(AWG_BUFFER_SIZE)
error_signal = np.zeros(AWG_BUFFER_SIZE) # generate copy of input signal with some errors based on DEBUG_BER
clock_signal = np.zeros(AWG_BUFFER_SIZE)
bits = np.zeros(NUM_WORDS*SYMS_PER_WORD)
clock_word = make_word(2**SYMS_PER_WORD-1, SYMS_PER_WORD, AWG_SHORT_PULSES[0], BIT_RATE, AWG_FSAMP)
# generate test vectors
for i in range(NUM_WORDS):
    val = random.randint(0, NUM_WORDS-1)
    input_signal[i*AWG_WORD_SIZE:(i+1)*AWG_WORD_SIZE] = make_word(val, SYMS_PER_WORD, AWG_SHORT_PULSES[0], BIT_RATE, AWG_FSAMP)
    new_val = []
    for b in range(SYMS_PER_WORD):
        if random.random() < DEBUG_BER:
            # flip each bit w.p. DEBUG_BER
            new_val.append(1 - (val & 1))
        else:
            new_val.append(val & 1)
        bits[i*SYMS_PER_WORD + b] = val & 1
        val >>= 1
    error_signal[i*AWG_WORD_SIZE:(i+1)*AWG_WORD_SIZE] = make_word(sum(i*2**n for n,i in enumerate(new_val)), SYMS_PER_WORD, AWG_SHORT_PULSES[0], BIT_RATE, AWG_FSAMP)
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

# first send a bunch of zeros to check noise level on DAQ inputs
daq.set_capture_cycles(1)
for c in AWG_CHANNELS:
    # set capture sizes and cycle counts so we don't gather too much data
    # weird bugs happen when changing the DAQ capture points,
    # so just decrease the number of capture cycles
    # set delays and amplitudes to zero; we don't want to send any data
    awg.set_channel_delay(c, 0)
    awg.set_channel_amplitude(c, 0)
    awg.set_buffer_contents(c, np.zeros(AWG_BUFFER_SIZE))
# actually capture noise data
try:
    awg.launch_channels(sum(2**(c-1) for c in AWG_CHANNELS), 1)
    daq_data = daq.capture(sum(2**(c-1) for c in DAQ_CHANNELS))
    # unit conversion to V and get noise properties
    for n in range(len(DAQ_CHANNELS)):
        daq_data[n] = (daq_data[n]/2**15)*DAQ_FULL_SCALE[n]
        noise_stddev[n] = np.std(daq_data[n])
        noise_maxabs[n] = np.max(np.abs(daq_data[n]))
except pxi_modules.KeysightException as e:
    print("Encountered fatal exception when commanding Keysight equipment, exiting now")
    print(e)
    awg.stop()
    daq.stop()
    exit()
print("channel ", c, " noise levels (stddev) = ", noise_stddev)
print("channel ", c, " noise levels (maxabs) = ", noise_maxabs)
# reset capture sizes to original params
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
        if DEBUG_BER_CHECKER and c == 2:
            awg.set_buffer_contents(c, error_signal)
        else:
            awg.set_buffer_contents(c, clock_signal)
for cfg in range(N_CONFIGS):
    amplitudes = TEST_CONFIGURATIONS[cfg]
    # set up amplitudes
    for n,c in enumerate(AWG_CHANNELS):
        awg.set_channel_amplitude(c, amplitudes[n])
    # capture data
    try:
        awg.launch_channels(sum(2**(c-1) for c in AWG_CHANNELS), TEST_CYCLES)
        daq_data = daq.capture(sum(2**(c-1) for c in DAQ_CHANNELS))
        # unit conversion to V
        for n in range(len(DAQ_CHANNELS)):
            daq_data[n] = (daq_data[n]/2**15)*DAQ_FULL_SCALE[n]
    except pxi_modules.KeysightException as e:
        print("Encountered fatal exception when commanding Keysight equipment, exiting now")
        print(e)
        awg.stop()
        daq.stop()
        exit()
    except Exception as e:
        print("Caught generic exception, closing AWG/DAQ now")
        print(e)
        awg.stop()
        daq.stop()
        exit()

    # calculate bit error rate
    t0 = time.process_time()
    for n,channel in enumerate(DAQ_CHANNELS):
        if DEBUG_BER_CHECKER and channel != 2:
            continue
        # zscore reflects the SNR of the readout circuit --- 200-500 is suitable for loopback,
        # but most likely a lower threshold (e.g. 20sigma) is needed for a real test
        #zscore = 50
        #threshold = zscore*np.std(daq_data[n,:1000])
        #threshold = max(daq_data[n])*0.4
        # if threshold is too high, then there will be incorrectly many 1->0 errors
        # if threshold is too low, then there will be incorrectly many 0->1 errors
        threshold = 200*noise_stddev[n]
        # get cross correlation between daq data and 2x decimated awg signal to help with alignment for BER calc
        corr = sigproc.correlate(daq_data[n], np.tile(input_signal[::2],TEST_CYCLES))
        # corr /= np.max(corr)
        lags = sigproc.correlation_lags(len(daq_data[n]), len(input_signal[::2])*TEST_CYCLES)
        best_lag = lags[np.argmax(corr)]
        daq_data_delayed = np.zeros(DAQ_BUFFER_SIZE*TEST_CYCLES)
        if best_lag > 0:
            daq_data_delayed[:len(daq_data[n])-best_lag] = daq_data[n,best_lag:]
        else:
            print("WARNING: got unphysical optimal lag, BER rate estimate has been set to 1 (indicating > 0.5 BER)")
            one_to_zero[cfg,n] = -1
            zero_to_one[cfg,n] = -1
            ber[cfg,n] = 1
            continue
        for bit in range(SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES):
            # each symbol is DAQ_FSAMP//BIT_RATE samples
            search_min = max(0, bit*(DAQ_FSAMP//BIT_RATE) - DAQ_FSAMP//(2*BIT_RATE))
            search_max = min(DAQ_BUFFER_SIZE*TEST_CYCLES, bit*(DAQ_FSAMP//BIT_RATE) + DAQ_FSAMP//(2*BIT_RATE))
            if bits[bit % (SYMS_PER_WORD*NUM_WORDS)]:
                # make sure that there wasn't a 1 -> 0 error
                # check that there's a peak in the vicinity
                if max(daq_data_delayed[search_min:search_max]) < threshold:
                    one_to_zero[cfg,n] += 1
            else:
                # make sure there wasn't a 0 -> 1 error
                # check that there aren't any peaks in the vicinity
                if max(daq_data_delayed[search_min:search_max]) > threshold:
                    zero_to_one[cfg,n] += 1
        ber[cfg,n] = (one_to_zero[cfg,n] + zero_to_one[cfg,n])/(SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES)
    t1 = time.process_time()
    ber_calc_time = t1 - t0
    print("it took ", ber_calc_time, " to calculate the bit error rate")

print("num (1->0 errors) = ", one_to_zero)
print("num (0->1 errors) = ", zero_to_one)
print("BER = ", ber)

# very important to close AWG, otherwise the buffers will not get properly flushed
print("closing AWG and DAQ")
awg.stop()
daq.stop()

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
    axs[0].plot(tvec_awg, np.tile(input_signal,TEST_CYCLES), label = 'awg output', color=colors[1])
    axs[0].legend()
    axs[0].set_xlabel(f"t [{siprefix[t_units]}s]")
    axs[0].set_ylabel("V [V]")
    axs[1].plot(lags/DAQ_FSAMP/t_units, corr, label = "corr(daq_data, awg_data)", color=colors[0])
    axs[1].legend()
    axs[1].set_xlabel(f"t [{siprefix[t_units]}s]")
    axs[2].plot(tvec_daq, daq_data_delayed, label = 'daq data', color=colors[0])
    axs[2].plot(tvec_awg[::2], np.tile(input_signal[::2],TEST_CYCLES), label = 'awg output', color=colors[1])
    axs[2].plot(np.linspace(0,(SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES-1)/BIT_RATE/t_units,SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES), 0.6*np.tile(bits,TEST_CYCLES), '.', label = 'bits', color=colors[2])
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
