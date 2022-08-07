# This program has been adapted from Keysight's SD1 docs
# It allows for synchronous capture of a device's response given an arbitrary stimulus
# ----------

import sys
import time
import matplotlib.pyplot as plt
import numpy as np

sys.path.append('C:\Program Files (x86)\Keysight\SD1\Libraries\Python')
import keysightSD1
from keysightSD1 import AIN_Impedance as imp
from keysightSD1 import AIN_Coupling as cpl
from keysightSD1 import SD_TriggerModes as trg
from keysightSD1 import SD_AIN_TriggerMode as atrg
from keysightSD1 import SD_TriggerBehaviors as trgb
from keysightSD1 import SD_TriggerExternalSources as trgext

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
        self.trigger_mode = trg.SWHVITRIG

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
        start_idx = self.channel_delays[channel] % 10
        data_delayed[start_idx:] = data[:self.buffer_length-start_idx]
        self.waveforms[self.channels[channel],:] = np.array(data_delayed)
        # set up channel
        self.awg.channelWaveShape(channel, keysightSD1.SD_Waveshapes.AOU_AWG)
        self.awg.channelAmplitude(channel, amplitude)
        # very important, behavior seems wonky if the AWG queue fills up
        self.awg.AWGflush(channel)

    def set_trigger_mode(self, trigger_mode):
        self.trigger_mode = trigger_mode
        if trigger_mode == trg.EXTTRIG or trigger_mode == trg.EXTTRIG_CYCLE:
            for channel in self.channels.keys():
                sync = 0
                error = self.awg.AWGtriggerExternalConfig(channel, trgext.TRIGGER_PXI, trgb.TRIGGER_RISE, sync)
                if error < 0:
                    raise Exception(f"AWG external trigger configuration failed (channel {channel}): {error}")

    def launch_channels(self, awg_mask, num_cycles):
        for channel in self.channels.keys():
            self.awg.AWGflush(channel)
            wave = keysightSD1.SD_Wave()
            waveformID = channel
            wave.newFromArrayDouble(0, self.waveforms[self.channels[channel]])
            error = self.awg.waveformLoad(wave, waveformID, 0)
            if error < 0:
                raise Exception(f"AWG load waveform error (channel {channel}): {error}")
            error = self.awg.AWGqueueWaveform(channel, waveformID, self.trigger_mode, self.channel_delays[channel]//10, num_cycles, 0)
            if error < 0:
                raise Exception(f"AWG queue waveform error (channel {channel}): {error}")
        # start all AWGs simultaneously
        error = self.awg.AWGstartMultiple(awg_mask)
        if error < 0:
            raise Exception(f"AWG start error: {error}")

    def stop(self):
        self.awg.AWGstopMultiple(sum(2**(c-1) for c in self.channels.keys()))
        time.sleep(0.1)
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
        self.trigger_mode = trg.SWHVITRIG

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
    
    def set_trigger_mode(self, trigger_mode, trigger_delay=0, analog_trig_chan=1, analog_trig_type=atrg.RISING_EDGE, threshold=0):
        self.trigger_mode = trigger_mode
        for c in self.channels.keys():
            self.daq.DAQconfig(c, self.points_per_cycle, self.cycles, trigger_delay, trigger_mode)
            # analog trigger requires setting the proper mask for each channel and
            # configuring the analog triggering condition for the analog trigger channel.
            # a non-one-hot mask (e.g. 0011 or 0101) will trigger when any of the enabled
            # channels in the mask meet the trigger condition as configured by channelTriggerConfig()
            if trigger_mode == trg.ANALOGTRIG:
                self.daq.DAQanalogTriggerConfig(c, 2**(analog_trig_chan-1))
                if c == analog_trig_chan:
                    self.daq.channelTriggerConfig(c, analog_trig_type, threshold)
            if trigger_mode == trg.EXTTRIG:
                self.daq.DAQdigitalTriggerConfig(c, trgext.TRIGGER_PXI, trgb.TRIGGER_RISE)
    
    def capture(self, daq_mask):
        data = np.zeros((len(self.channels),self.cycles*self.points_per_cycle))
        self.daq.DAQflushMultiple(daq_mask)
        self.daq.DAQstartMultiple(daq_mask)
        capture_timeout = 0.1 # 100ms
        # if trigger mode is SW/HVI or PXI/EXT, then manually send a trigger for each capture cycle
        if self.trigger_mode == trg.SWHVITRIG:
            for cycle in range(self.cycles):
                self.daq.DAQtriggerMultiple(daq_mask)
                for channel in self.channels.keys():
                    self.wait_until_points_read(channel, self.points_per_cycle*(cycle+1), capture_timeout)
        if self.trigger_mode == trg.EXTTRIG:
            # by default use PXI0 for everything
            PXI = 0
            for cycle in range(self.cycles):
                # PXI triggers are active low (i.e. 0 represents a triggered state)
                self.daq.PXItriggerWrite(PXI, 1)
                self.daq.PXItriggerWrite(PXI, 0)
                self.daq.PXItriggerWrite(PXI, 1)
                for channel in self.channels.keys():
                    self.wait_until_points_read(channel, self.points_per_cycle*(cycle+1), capture_timeout)
        # capture data
        for channel in self.channels.keys():
            if self.trigger_mode != trg.SWHVITRIG:
                self.wait_until_points_read(channel, self.points_per_cycle*self.cycles, capture_timeout)
        self.daq.DAQstopMultiple(daq_mask)
        read_timeout = 100 # timeout in ms; 0 is infinite
        for channel in self.channels.keys():
            data[self.channels[channel]] = self.daq.DAQread(channel, self.cycles*self.points_per_cycle, read_timeout)
        return data

    def stop(self):
        # just in case it's not stopped, stop the DAQ
        self.daq.DAQstopMultiple(sum(2**(c-1) for c in self.channels.keys()))
        time.sleep(0.1)
        self.daq.close()

CHASSIS = 1

# AWG constants
AWG_CHANNELS = [1, 2, 3, 4]
AWG_DELAYS = [0, 5, 10, 15] # delay in ns
AWG_AMPLITUDE = [0.7, 0.7, 0.7, 0.7] # full scale in V
AWG_BUFFER_LENGTH = 1000

AWG_PULSES = [
    [0.5, 1, 0.7, 0.2], # best (low ripple pulse) we have so far: 2.3ns FWHM
    [0.75, 1, 0.3, 0.1], # this is pretty good, moderate ripple, 1.9ns FWHM
    [0.6, 1, 0.5, 0.1], # this is pretty good, low ripple, 2.1ns FWHM
    [0.7, 1, 0.4, 0.0] # also pretty good, moderate ripple
]

# DAQ constants
TSAMP = 2e-9
DAQ_CHANNELS = [1, 2, 3, 4]
DAQ_POINTS_PER_CYCLE = 1000
DAQ_CYCLES = 3
DAQ_TRIG_DELAY = 0
DAQ_FULL_SCALE = [1, 1, 1, 1] # full scale in V
DAQ_IMPEDANCE = [imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50, imp.AIN_IMPEDANCE_50]
DAQ_COUPLING = [cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC, cpl.AIN_COUPLING_DC]

awg_data = [
    ((AWG_PULSES[0] + [0]*16)*4 + [0]*20)*10,
    ((AWG_PULSES[0] + [0]*16)*4 + [0]*20)*10,
    ((AWG_PULSES[0] + [0]*16)*4 + [0]*20)*10,
    #((AWG_PULSES[0] + [0]*16)*4 + [0]*20)*10
    ([1]*98 + [0.5, 0.1, -0.1, -0.5] + [-1]*98)*5
]

# add AWG and DAQ channels, set AWG buffer contents
awg = AWG("M3202A", 1, 5, AWG_BUFFER_LENGTH)
awg.add_channels(AWG_CHANNELS)
for n,c in enumerate(AWG_CHANNELS):
    awg.set_channel_delay(c, AWG_DELAYS[n])
    awg.set_buffer_contents(c, awg_data[n], AWG_AMPLITUDE[n])

daq = DAQ("M3102A", 1, 7, DAQ_POINTS_PER_CYCLE, DAQ_CYCLES)
daq.add_channels(DAQ_CHANNELS, DAQ_FULL_SCALE, DAQ_IMPEDANCE, DAQ_COUPLING)

# set up triggering modes
# by default, use EXTTRIG (or EXTTRIG_CYCLE in the case of the AWG) to allow
# for synchronization of the AWG and DAQ so no samples are lost/missed
awg.set_trigger_mode(trg.EXTTRIG_CYCLE)
daq.set_trigger_mode(trg.EXTTRIG)
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

awg.launch_channels(sum(2**(c-1) for c in AWG_CHANNELS), DAQ_CYCLES)
daq_data = daq.capture(sum(2**(c-1) for c in DAQ_CHANNELS))

# very important to close AWG, otherwise the buffers will not get properly flushed
print("closing AWG and DAQ")
awg.stop()
daq.stop()

print("plotting results")
plt.figure()
for n,channel in enumerate(DAQ_CHANNELS):
    plt.plot(np.linspace(0,(DAQ_POINTS_PER_CYCLE*DAQ_CYCLES-1)*TSAMP*1e9,DAQ_POINTS_PER_CYCLE*DAQ_CYCLES), (daq_data[n]/2**15)*DAQ_FULL_SCALE[n], label = f'ch{channel}')
plt.legend()
plt.xlabel("t [ns]")
plt.ylabel("V [V]")
plt.show()