# This sample program shows how Python uses the Keysight SD1 Programming Libraries.
# ----------
# Import required system components
import sys
import time
import matplotlib.pyplot as plt
import numpy as np

sys.path.append('C:\Program Files (x86)\Keysight\SD1\Libraries\Python')
import keysightSD1
from keysightSD1 import AIN_Impedance as imp
from keysightSD1 import AIN_Coupling as cpl
from keysightSD1 import SD_TriggerModes as trg

class AWG:

    def __init__(self, model_name, chassis, slot, buffer_length):
        self.model_name = model_name
        self.chassis = chassis
        self.slot = slot
        self.awg = keysightSD1.SD_AOU()
        resp = self.awg.openWithSlot(self.model_name, self.chassis, self.slot)
        if resp < 0:
            raise Exception(f"Couldn't connect to AWG {self.model_name} in slot {self.slot}")
        self.channels = {}
        self.channel_delays = {}
        self.waveforms = np.array([])
        self.buffer_length = buffer_length

    def add_channels(self, channels):
        for c in channels:
            if len(self.channels) == 0:
                self.waveforms = np.array([[0]*self.buffer_length])
            else:
                self.waveforms = np.append(self.waveforms, [np.zeros(self.buffer_length)], axis=0)
            self.channels[c] = len(self.channels)
            self.channel_delays[c] = 0

    def set_channel_delay(self, channel, delay):
        self.channel_delays[channel] = delay

    def set_buffer_length(self, buffer_length):
        # zero pads buffer
        if buffer_length < self.buffer_length:
            self.waveforms = self.waveforms[:,0:buffer_length-1]
        else:
            self.waveforms = np.append(self.waveforms, np.zeros(shape=(self.waveforms.shape[0],buffer_length-self.buffer_length)), axis=1)
        self.buffer_length = buffer_length
    
    def set_buffer_contents(self, channel, data, amplitude):
        # do fine-grained (i.e. Tsamp resolution) delay adjustment, since delay parameter to AWGqueue is in multiples of 5 samples
        data_delayed = np.zeros(self.buffer_length)
        start_idx = self.channel_delays[channel] % 5
        data_delayed[start_idx:] = data[:self.buffer_length-start_idx]
        self.waveforms[self.channels[channel],:] = np.array(data_delayed)
        # set up channel
        self.awg.channelWaveShape(channel, keysightSD1.SD_Waveshapes.AOU_AWG)
        self.awg.channelAmplitude(channel, amplitude)
        # very important, behavior seems wonky if the AWG queue fills up
        self.awg.AWGflush(channel)

    def launch_channels(self, awg_mask):
        for channel in self.channels.keys():
            wave = keysightSD1.SD_Wave()
            waveformID = channel
            wave.newFromArrayDouble(0, self.waveforms[self.channels[channel]])
            error = self.awg.waveformLoad(wave, waveformID, 0)
            if error < 0:
                raise Exception(f"AWG load waveform error (channel {channel}):", error)
            error = self.awg.AWGqueueWaveform(channel, waveformID, 0, self.channel_delays[channel]//5, 0, 0)
            if error < 0:
                raise Exception(f"AWG queue waveform error (channel {channel}):", error)
        # start all AWGs simultaneously
        error = self.awg.AWGstartMultiple(awg_mask)
        if error < 0:
            raise Exception(f"AWG start error:", error)

    def stop(self):
        for channel in self.channels.keys():
            self.awg.AWGstop(channel)
        self.awg.close()

class DAQ:
    
    def __init__(self, model_name, chassis, slot, points_per_cycle, cycles):
        self.model_name = model_name
        self.chassis = chassis
        self.slot = slot
        self.daq = keysightSD1.SD_AIN()
        resp = self.daq.openWithSlot(self.model_name, self.chassis, self.slot)
        if resp < 0:
            raise Exception(f"Couldn't connect to DAQ {self.model_name} in slot {self.slot}")
        self.channels = {}
        self.full_scale = {}
        self.points_per_cycle = points_per_cycle
        self.cycles = cycles

    def wait_until_points_read(self, channel, points, timeout):
        t0 = time.time()
        time_elapsed = 0
        points_read = 0
        while (points_read < points) and (time_elapsed < timeout):
            points_read = self.daq.DAQcounterRead(channel)
            if points_read < points:
                time.sleep(.01) #input argument of time.sleep is in seconds
            time_elapsed = time.time() - t0
    
    def add_channels(self, channels, full_scale, impedance, coupling):
        for n,c in enumerate(channels):
            self.channels[c] = len(self.channels)
            self.full_scale[c] = full_scale[n]
            self.daq.channelInputConfig(c, self.full_scale[c], impedance[n], coupling[n])
    
    def swcapture(self, daq_mask):
        for c in self.channels.keys():
            self.daq.DAQconfig(c, self.points_per_cycle, self.cycles, 0, trg.SWHVITRIG)
        self.daq.DAQflushMultiple(daq_mask)
        self.daq.DAQstartMultiple(daq_mask)
        data = np.zeros((len(self.channels),self.cycles,self.points_per_cycle))
        for cycle in range(self.cycles):
            self.daq.DAQtriggerMultiple(daq_mask)
            for channel in self.channels.keys():
                timeout = 0.1 # 100ms
                self.wait_until_points_read(channel, self.points_per_cycle, timeout)
                read_timeout = 100 # 0 is inifinite
                data[self.channels[channel], cycle] = self.daq.DAQread(channel, self.points_per_cycle, read_timeout)
        self.daq.DAQstopMultiple(daq_mask)
        return data
    
    def set_trigger_mode(self, channel, trigger_mode):
        self.daq.DAQconfig(channel, self.points_per_cycle, self.cycles, 0, trigger_mode[self.channels[channel]])

    def stop(self):
        self.awg.close()

CHASSIS = 1

# AWG constants
AWG_CHANNELS = [1,2]
AWG_DELAYS = [0,0] # delay in ns
AWG_AMPLITUDE = [1.0, 1.0] # full scale in V
AWG_BUFFER_LENGTH = 1000

AWG_PULSES = [
    [0.5, 1, 0.7, 0.2] # best (low ripple pulse) we have so far: 2.3ns FWHM
    #[0.75, 1, 0.3, 0.1], # this is pretty good, moderate ripple, 1.9ns FWHM
    #[0.6, 1, 0.5, 0.1], # this is pretty good, low ripple, 2.1ns FWHM
    #[0.7, 1, 0.4, 0.0], # also pretty good, moderate ripple
]

# DAQ constants
TSAMP = 2e-9
DAQ_CHANNELS = [1, 2]
DAQ_POINTS_PER_CYCLE = 500
DAQ_CYCLES = 1
DAQ_TRIG_DELAY = 0
DAQ_FULL_SCALE = [1, 1]
DAQ_IMPEDANCE = [imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50]
DAQ_COUPLING = [cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC]
#DAQ_TRIG = [trg.SWHVITRIG, trg.SWHVITRIG]
data = [
    ((AWG_PULSES[0] + [0]*16)*3 + [0]*40)*10,
    ([1]*48 + [0.5, 0.1, -0.1, -0.5] + [-1]*48)*10
]
awg = AWG("M3202A", 1, 5, AWG_BUFFER_LENGTH)
awg.add_channels(AWG_CHANNELS)
for n,c in enumerate(AWG_CHANNELS):
    awg.set_channel_delay(c, AWG_DELAYS[n])
    awg.set_buffer_contents(c, data[n], AWG_AMPLITUDE[n])

daq = DAQ("M3102A", 1, 7, DAQ_POINTS_PER_CYCLE, DAQ_CYCLES)
daq.add_channels(DAQ_CHANNELS, DAQ_FULL_SCALE, DAQ_IMPEDANCE, DAQ_COUPLING)

awg.launch_channels(3)
data = daq.swcapture(3)

print("plotting results")
plt.figure()
for n,channel in enumerate(DAQ_CHANNELS):
    for cycle in range(DAQ_CYCLES):
        plt.plot(np.linspace(0,(DAQ_POINTS_PER_CYCLE-1)*TSAMP*1e6,DAQ_POINTS_PER_CYCLE), (data[n, cycle]/2**15)*DAQ_FULL_SCALE[n], label = f'ch{channel} - {cycle}')
plt.legend()
plt.xlabel("t [us]")
plt.ylabel("V [V]")
plt.show()
# ----------
# © Keysight Technologies, 2020
# All rights reserved.
# You have a royalty-free right to use, modify, reproduce and distribute
# this Sample Application (and/or any modified # version) in any way you find useful,
# provided that you agree that Keysight Technologies has no warranty, obligations or liability 
# for any Sample Application Files.
#
# Keysight Technologies provides programming examples for illustration only.
# This sample program assumes that you are familiar with the programming language 
# being demonstrated and the tools used to create and debug procedures.
# Keysight Technologies support engineers can help explain the functionality 
# of Keysight Technologies software components and associated commands,
# but they will not modify these samples to provide added functionality or 
# construct procedures to meet your specific needs.
# ----------
