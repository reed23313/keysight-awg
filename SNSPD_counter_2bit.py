# This program has been adapted from Keysight's SD1 docs
# It allows for synchronous capture of a device's response given an arbitrary stimulus
# channels 1 and 2 emulate ntron/SNSPD pulses
# pulse height of SNSPD waveforms is inconsistent
# this is mostly due to the low sample rate of the 500MS/s digitizer
# ----------

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

SAVE = 1
PLOT = 0

CHASSIS = 1

AWG_TSAMP = 1e-9 # 1GS/s
AWG_CHANNELS = [1, 2, 3, 4]
#AWG_DELAYS = [0, 2, 7, 12] # delay in ns
AWG_DELAYS = [0, 0, 0, 0] # delay in ns
# AWG_BUFFER_LENGTH = 1000 # number of samples to hold in buffer# AWG constants

# DAQ constants
DAQ_TSAMP = 2e-9 # 500MS/s
DAQ_CHANNELS = [1, 2, 3, 4]
DAQ_CYCLES = 1 # number of acquisition cycles/frames
DAQ_TRIG_DELAY = 180 # set the delay until capturing samples from when the trigger condition is met
DAQ_FULL_SCALE = [1, 1, 1, 1] # full scale in V
DAQ_IMPEDANCE = [imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50]
DAQ_COUPLING = [cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC]


TEST_CYCLES = 1 # 1Msamp
NUM_WORDS = 50 # should be suitably large enough so that we get good word entropy
SYMS_PER_WORD = 32 # word size in symbols
AWG_FSAMP = int(1e9) # 1GS/s
DAQ_FSAMP = int(500e6) # 500MS/s
BIT_RATE = int(1e6) # sym/s (maybe up to 150MHz; 100MHz looks good however)
AWG_WORD_SIZE = SYMS_PER_WORD*(AWG_FSAMP//BIT_RATE) # word size in samples
DAQ_WORD_SIZE = SYMS_PER_WORD*(DAQ_FSAMP//BIT_RATE) # word size in samples
AWG_BUFFER_SIZE = NUM_WORDS*AWG_WORD_SIZE # buffer size in samples
DAQ_BUFFER_SIZE = NUM_WORDS*DAQ_WORD_SIZE
PULSE_SHAPE = AWG_SHORT_PULSES[0]

BASE = 2


N_VIN = 1
N_VBIAS1 = 10
N_VBIAS2 = 10
N_VBIAS3 = 1
V_VIN_MIN = 0
V_VIN_MAX = 0
V_BIAS1_MIN = 0.2
V_BIAS1_MAX = 0.5
V_BIAS2_MIN = 0.2
V_BIAS2_MAX = 0.5
V_BIAS3_MIN = 0.28
V_BIAS3_MAX = 0.28
V_VIN_SWEEP = np.linspace(V_VIN_MIN,V_VIN_MAX,N_VIN)
V_BIAS1_SWEEP = np.linspace(V_BIAS1_MIN,V_BIAS1_MAX,N_VBIAS1)
V_BIAS2_SWEEP = np.linspace(V_BIAS2_MIN,V_BIAS2_MAX,N_VBIAS2)
V_BIAS3_SWEEP = np.linspace(V_BIAS3_MIN,V_BIAS3_MAX,N_VBIAS3)

tvec_daq = np.linspace(0,(DAQ_BUFFER_SIZE*TEST_CYCLES-1)/DAQ_FSAMP,DAQ_BUFFER_SIZE*TEST_CYCLES)
tvec_awg = np.linspace(0,(AWG_BUFFER_SIZE*TEST_CYCLES-1)/AWG_FSAMP,AWG_BUFFER_SIZE*TEST_CYCLES)

i1 = 0
i2 = 0
i3 = 0
i4 = 0

# generate configs
N_CONFIGS = N_VIN*N_VBIAS1*N_VBIAS2*N_VBIAS3
peaks1 = np.zeros((N_CONFIGS, DAQ_BUFFER_SIZE*TEST_CYCLES//4),dtype=np.int32)
peaks2 = np.zeros((N_CONFIGS, DAQ_BUFFER_SIZE*TEST_CYCLES//4),dtype=np.int32)
peaks3 = np.zeros((N_CONFIGS, DAQ_BUFFER_SIZE*TEST_CYCLES//4),dtype=np.int32)
peaks4 = np.zeros((N_CONFIGS, DAQ_BUFFER_SIZE*TEST_CYCLES//4),dtype=np.int32)
TEST_CONFIGURATIONS = np.zeros((N_CONFIGS,4))
i = 0
for v_in in V_VIN_SWEEP:
    for v_bias1 in V_BIAS1_SWEEP:
        for v_bias2 in V_BIAS2_SWEEP:
            for v_bias3 in V_BIAS3_SWEEP:
                TEST_CONFIGURATIONS[i,:] = np.array([v_in, v_bias1, v_bias2, v_bias3])
                i = i + 1

# add AWG and DAQ channels, set AWG buffer contents
awg = pxi_modules.AWG("M3202A", 1, 5, AWG_BUFFER_SIZE)
awg.add_channels(AWG_CHANNELS)
    # awg.set_channel_amplitude(c, AWG_AMPLITUDE[n])
    # awg.set_channel_offset(c, 0)
    # awg.set_buffer_contents(c, awg_data[n])

daq = pxi_modules.DAQ("M3102A", 1, 7, DAQ_BUFFER_SIZE, DAQ_CYCLES)
daq.add_channels(DAQ_CHANNELS, DAQ_FULL_SCALE, DAQ_IMPEDANCE, DAQ_COUPLING)

# set up triggering modes
# by default, use EXTTRIG (or EXTTRIG_CYCLE in the case of the AWG) to allow
# for synchronization of the AWG and DAQ so no samples are lost/missed
awg.set_trigger_mode(trg.EXTTRIG_CYCLE)
daq.set_trigger_mode(trg.EXTTRIG, trigger_delay = DAQ_TRIG_DELAY)

lags = np.zeros(len(DAQ_CHANNELS), dtype=np.int32)


    # clocks just get 1111....
    # clock_signal[i*AWG_WORD_SIZE:(i+1)*AWG_WORD_SIZE] = clock_word

bias_zerobits = np.zeros(AWG_BUFFER_SIZE)
for n,c in enumerate(AWG_CHANNELS):
    if c != 1:
        awg.set_buffer_contents(c, bias_zerobits)
    # set delay
    awg.set_channel_delay(c, AWG_DELAYS[n])


try:
    t0 = time.time()


    for cfg in range(N_CONFIGS):

        bits = np.zeros(NUM_WORDS*SYMS_PER_WORD)
        input_signal = np.zeros(AWG_BUFFER_SIZE)
        # for i in range(NUM_WORDS):
        #     val = random.randint(0, 2**(SYMS_PER_WORD//BASE**2-1))
        #     if i == NUM_WORDS - 1:
        #         # send 0 for the last word so next cycle doesn't have an accidental bit flip at the beginning due to a held over circulating current
        #         # i.e., make sure that there are an equal number of pulses for each channel
        #         val = 0
        #     seq = make_word(val, SYMS_PER_WORD//BASE**2, PULSE_SHAPE, BIT_RATE, AWG_FSAMP)
        #     input_signal[i*AWG_WORD_SIZE:(i+1)*AWG_WORD_SIZE] = seq*BASE**2
        #     for b in range(SYMS_PER_WORD):
        #         bits[i*SYMS_PER_WORD + b] = val & 1
        #         val >>= 1
        awg.set_buffer_contents(1, input_signal)

        amplitudes = TEST_CONFIGURATIONS[cfg]
        # set up amplitudes
        for n,c in enumerate(AWG_CHANNELS):
            if c == 1:
                awg.set_channel_amplitude(c, amplitudes[n])
                awg.set_channel_offset(c, 0)
            else:
                awg.set_channel_offset(c, amplitudes[n])
        # capture data
        #awg.launch_channels(sum(2**(c-1) for c in AWG_CHANNELS), TEST_CYCLES)
        # if cfg == 0:
        if cfg == 0:
            awg.launch_channels(sum(2**(c-1) for c in AWG_CHANNELS), TEST_CYCLES)
        else:
            awg.launch_channels_start_only(sum(2**(c-1) for c in AWG_CHANNELS))
        # else:
        # awg.launch_channels_start_only(sum(2**(c-1) for c in AWG_CHANNELS))

        daq_data = daq.capture(sum(2**(c-1) for c in DAQ_CHANNELS))


        ##### change it 
        # if abs(max(input_signal)) > 0.1:  ###############chnage this value, it should be lower because lower voltages 
        #     min_peak_height_1 = 0.5*max(input_signal)
        # else:
        #     min_peak_height_1 = 1000*max(input_signal)
       
        if abs(max(daq_data[0])) > 1500:
            min_peak_height_1 = 0.7*max(daq_data[0])
        else:
            min_peak_height_1 = 1000*max(daq_data[0])
        if abs(max(daq_data[1])) > 800:
            min_peak_height_2 = 0.7*max(daq_data[1])
        else:
            min_peak_height_2 = 1000*max(daq_data[1])
        if abs(max(daq_data[2])) > 800:
            min_peak_height_3 = 0.7*max(daq_data[2])
        else:
            min_peak_height_3 = 1000*max(daq_data[2])
        if abs(max(daq_data[3])) > 800:  ###############chnage this value, it should be lower because lower voltages 
            min_peak_height_4 = 0.5*max(daq_data[3])
        else:
            min_peak_height_4 = 1000*max(daq_data[3])
        

        c1 = daq_data[0]
        c2 = daq_data[1]
        c3 = daq_data[2]
        c4 = daq_data[3]
        
        # distance=(7*DAQ_FSAMP)//(8*BIT_RATE)

        p1, _ = sigproc.find_peaks(c1, prominence=min_peak_height_1, distance=4)
        p2, _ = sigproc.find_peaks(c2, prominence=min_peak_height_2, distance=4)
        p3, _ = sigproc.find_peaks(c3, prominence=min_peak_height_3, distance=4)
        p4, _ = sigproc.find_peaks(c4, prominence=min_peak_height_4, distance=4)

        if PLOT == 1:
            print("plotting results")
            fig, axs = plt.subplots(4,1)
            axs[0].plot(c1)
            axs[0].plot(p1,c1[p1])
            axs[1].plot(c2)
            axs[1].plot(p2,c2[p2])
            axs[2].plot(c3)
            axs[2].plot(p3,c3[p3])
            axs[3].plot(c4)
            axs[3].plot(p4,c4[p4])
            fig.tight_layout()
            plt.show()
        
        if SAVE == 1:
            peaks1[i1,0:len(p1)] = p1 
            i1 += 1
            peaks2[i2,0:len(p2)] = p2
            i2 += 1
            peaks3[i3,0:len(p3)] = p3 
            i3 += 1
            peaks4[i4,0:len(p4)] = p4 
            i4 += 1
            # peaks2 = np.append(peaks2,p2)
            # peaks3 = np.append(peaks3,p3)
            # peaks4 = np.append(peaks4,p4)
            # peaks1 = np.append(peaks1,[-2])
            # peaks2 = np.append(peaks2,[-2])
            # peaks3 = np.append(peaks3,[-2])
            # peaks4 = np.append(peaks4,[-2])

        # t_peak_1 = tvec[peaks1]
        # t_peak_2 = tvec[peaks2]
        # t_peak_3 = tvec[peaks3]
        # t_peak_4 = tvec[peaks4]


        # TIME_TOL = 1e-8
        
        # ### detection of first nTron
        # seq = np.zeros(1, len(t_peak_1))
        
        # for i in range(len(t_peak_1)):
        #     same_t = (abs(t_peak_2 - t_peak_1[i]) < TIME_TOL)
        #     seq[i] = any(same_t);
            


    
    # if not(any(seq)) || numel(t_peak_1) == 0
    #     right = 0; %finireeeee
    #     no_peaks = 1;
    # else
    #     right = mean(seq);
    # end

    # #### counting
    # seq = np.zeros(1, len(t_peak_1));
    # seq2 = np.zeros(1, len(t_peak_1));

    # for i in range(len(t_peak_1)):
    #     same_t = (abs(t_peak_2 - t_peak_1[i]) < TIME_TOL)
    #     seq[i] = any(same_t)
        
    #     same_t_2 = (abs(t_peak_3 - t_peak_1[i]) < TIME_TOL); 
    #     seq2[i] = any(same_t_2);

    

#     err_one2 = 0;
#     N = numel(t_peak_1);
#     right_vect = zeros(1, N);
#     right_vect2 = zeros(1, N);
#     base2 = base^2;
#     %|| sum(seq) >= 1.1*numel(seq)/base || sum(seq) <= 0.9*numel(seq)/base
#     if not(any(seq)) && not(any(seq2))
#         right = 0; %finireeeee
#     else
#         right_seq = zeros(1, N);
#         right_seq(base : base : numel(right_seq)) = 1;
#         right_seq2 = zeros(1, N);
#         right_seq2(base2 : base2 : numel(right_seq)) = 1;
        
#         num_one = 0;
#         for i = 1 : numel(right_seq)
            
#             if num_one == base
#                 num_one = 0;
#             end
#             if seq(i) == 1
#                 num_one = num_one + 1;
#             end
            
#             right_vect(i) = not(xor(right_seq(i), seq(i)));
            
#             if right_vect(i) == 0 && seq(i) == 1
#                 right_seq(i+1:end) = zeros(1, N-i);
#                 right_seq(i+base : base : numel(right_seq)) = 1;
                
#                 if err_one2 == 0
#                     right_seq2(i:end) = zeros(1, N-i+1);
#                     right_seq2(i+base*(base-num_one) : base2 : numel(right_seq2)) = 1;  
#                 end
#             end
            
#             if right_vect(i) == 0 && seq(i) == 0
#                 right_seq(i+1:end) = zeros(1, N-i);
#                 right_seq(i+1 : base : numel(right_seq)) = 1;
                
#                 if err_one2 == 0
#                     right_seq2(i:end) = zeros(1, N-i+1);
#                     right_seq2(i+1+base*(base-num_one-1) : base2 : numel(right_seq)) = 1;
#                 end
                
#             end
            
#             right_vect2(i) = not(xor(right_seq2(i), seq2(i)));
            
#             err_one2 = 0;
#             if right_vect2(i) == 0 && seq2(i) == 1
#                 right_seq2(i+1:end) = zeros(1, N-i);
#                 right_seq2(i+base2 : base2 : numel(right_seq2)) = 1;
#             end
            
#             if right_vect2(i) == 0 && seq2(i) == 0
#                 right_seq2(i+1:end) = zeros(1, N-i);
#                 right_seq2(i+1 : base2 : numel(right_seq2)) = 1;
#                 err_one2 = 1;
#             end
                    
#         end
#     end
#     if not(any(seq))
#         right_vect = zeros(1, N);
#     end
#     if not(any(seq2)) 
#         right_vect2 = zeros(1, N);
#     end
#     right = mean([right_vect right_vect2]);
# end
        # daq_data_delayed = np.zeros((len(DAQ_CHANNELS),DAQ_BUFFER_SIZE*TEST_CYCLES))
        # unit conversion to V
        # for n in range(len(DAQ_CHANNELS)):
        #     daq_data[n] = (daq_data[n]/2**15)*DAQ_FULL_SCALE[n]
            # daq_data_delayed[n,:DAQ_BUFFER_SIZE*TEST_CYCLES-lags[n]] = daq_data[n,lags[n]:]
        # calculate bit error rate
        #t0 = time.process_time()
        # for n,channel in enumerate(DAQ_CHANNELS):
        #     if channel != 3:
        #         continue
        #     # zscore reflects the SNR of the readout circuit --- 200-500 is suitable for loopback,
        #     # but most likely a lower threshold (e.g. 20 sigma) will be needed for a real test
        #     zscore = 20
        #     # if threshold is too high, then there will be incorrectly many 1->0 errors
        #     # if threshold is too low, then there will be incorrectly many 0->1 errors
        #     daq_symbols = np.zeros(SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES)
        #     #threshold = zscore*noise_stddev[n]
        #     threshold = 0.3*np.max(daq_data[n])
        #     #threshold = 20*energy_stddev[n]
        #     #daq_lpf = sigproc.lfilter(FILT_B, FILT_A, daq_data_delayed[n])
        #     #energy = np.trapz(np.square(np.reshape(daq_lpf, (NUM_WORDS*SYMS_PER_WORD*TEST_CYCLES,DAQ_FSAMP//BIT_RATE))), axis=1)
        #     #daq_symbols = energy > threshold
        #     # get peaks and then bin each peak into a symbol location/time
        #     # allow 3/4 of a clock/symbol period of separation between peaks
        #     peaks, _ = sigproc.find_peaks(daq_data[n], height=threshold, distance=(7*DAQ_FSAMP)//(8*BIT_RATE))
        #     daq_symbols[(peaks - lags[n]) // (DAQ_FSAMP//BIT_RATE)] = 1
        #     long_bits = np.tile(bits, TEST_CYCLES)
        #     zero_to_one[cfg,n] = np.sum(long_bits < daq_symbols)
        #     one_to_zero[cfg,n] = np.sum(long_bits > daq_symbols)
        #     ber[cfg,n] = (one_to_zero[cfg,n] + zero_to_one[cfg,n])/(SYMS_PER_WORD*NUM_WORDS*TEST_CYCLES)
        # #t1 = time.process_time()
        # #ber_calc_time = t1 - t0
        # #print("it took ", ber_calc_time, " to calculate the bit error rate")
        if (cfg % max(1,(N_CONFIGS//100))) == 0:
            print(cfg*100/N_CONFIGS, "% done")
    t1 = time.time()
    ber_calc_time = t1 - t0
    print("it took ", ber_calc_time, "s to calculate the bit error rate")
except pxi_modules.KeysightException as e:
    print("Encountered fatal exception when commanding Keysight equipment, exiting now")
    print(e)
    awg.stop()
    daq.stop()
    exit()
except (Exception, KeyboardInterrupt) as e:
    print("Caught generic exception, closing AWG/DAQ now")
    print(e)
    traceback.print_exc()
    awg.stop()
    daq.stop()
    exit()


# very important to close AWG, otherwise the buffers will not get properly flushed
print("closing AWG and DAQ")
awg.stop()
daq.stop()

# peaks1 = tvec[peaks1]
# peaks2 = tvec[peaks2]
# peaks3 = tvec[peaks3]
# peaks4 = tvec[peaks4]

    ##save in mat
if SAVE == 1:
    print("saving results")
    # save_data = {"ch1": daq_data[0], "ch2": daq_data[1], "ch3": daq_data[2], "ch4": daq_data[3], "tvec": tvec, "peaks1": peaks1, "peaks2": peaks2, "peaks3": peaks3, "peaks4": peaks4}
    save_data = {"tvec_awg": tvec_awg, "tvec_daq": tvec_daq, "peaks1": peaks1, "peaks2": peaks2, "peaks3": peaks3, "peaks4": peaks4, "vin_bias123": TEST_CONFIGURATIONS}
    mat = sio.savemat("SNSPDbias210_120_counter2bit_sweep.mat", save_data)




