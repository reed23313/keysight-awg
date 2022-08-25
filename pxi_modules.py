import sys
import time
import numpy as np

sys.path.append('C:\Program Files (x86)\Keysight\SD1\Libraries\Python')
import keysightSD1
from keysightSD1 import SD_TriggerModes as trg
from keysightSD1 import SD_AIN_TriggerMode as atrg
from keysightSD1 import SD_TriggerBehaviors as trgb
from keysightSD1 import SD_TriggerExternalSources as trgext

class KeysightException(Exception):
    pass

class AWG:

    def __init__(self, model_name, chassis, slot, buffer_length):
        self.model_name = model_name
        self.chassis = chassis
        self.slot = slot
        self.awg = keysightSD1.SD_AOU()
        resp = self.awg.openWithSlot(self.model_name, self.chassis, self.slot)
        if resp < 0:
            raise KeysightException(f"Couldn't connect to AWG {self.model_name} in slot {self.slot}")
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
            self.waveforms = self.waveforms[:,:buffer_length]
        else:
            self.waveforms = np.append(self.waveforms, np.zeros(shape=(self.waveforms.shape[0],buffer_length-self.buffer_length)), axis=1)
        self.buffer_length = buffer_length
    
    def set_buffer_contents(self, channel, data):
        # do fine-grained (i.e. Tsamp resolution) delay adjustment, since delay parameter to AWGqueue is in multiples of 10 samples
        data_delayed = np.zeros(self.buffer_length)
        #start_idx = self.channel_delays[channel] % 10
        # setting the delay in the AWGqueueWaveform method actually doesn't work properly when queuing multiple cycles
        # so instead, we'll just set the delay in samples here
        start_idx = self.channel_delays[channel]
        data_delayed[start_idx:] = data[:self.buffer_length-start_idx]
        if np.max(np.abs(data_delayed)) > 1:
            raise KeysightException(f"AWG buffer contents in channel ({channel}) cannot exceed |1| full scale")
        self.waveforms[self.channels[channel],:] = np.array(data_delayed)
        # set up channel
        self.awg.channelWaveShape(channel, keysightSD1.SD_Waveshapes.AOU_AWG)
        # very important, behavior seems wonky if the AWG queue fills up
        self.awg.AWGflush(channel)
    
    def set_channel_amplitude(self, channel, amplitude):
        self.awg.channelAmplitude(channel, amplitude)
    
    def set_channel_offset(self, channel, offset):
        self.awg.channelOffset(channel, offset)

    def set_trigger_mode(self, trigger_mode):
        self.trigger_mode = trigger_mode
        if trigger_mode == trg.EXTTRIG or trigger_mode == trg.EXTTRIG_CYCLE:
            for channel in self.channels.keys():
                sync = 0
                error = self.awg.AWGtriggerExternalConfig(channel, trgext.TRIGGER_PXI, trgb.TRIGGER_RISE, sync)
                if error < 0:
                    raise KeysightException(f"AWG external trigger configuration failed (channel {channel}): {error}")

    def launch_channels(self, awg_mask, num_cycles):
        for channel in self.channels.keys():
            self.awg.AWGflush(channel)
            wave = keysightSD1.SD_Wave()
            waveformID = channel
            error = wave.newFromArrayDouble(0, self.waveforms[self.channels[channel]])
            if error < 0:
                raise KeysightException(f"AWG create waveform error (channel {channel}): {error}")
            error = self.awg.waveformLoad(wave, waveformID, 0)
            if error < 0:
                raise KeysightException(f"AWG load waveform error (channel {channel}): {error}")
            # setting the delay in the AWGqueueWaveform method actually doesn't work properly when queuing multiple cycles
            #error = self.awg.AWGqueueWaveform(channel, waveformID, self.trigger_mode, self.channel_delays[channel]//10, num_cycles, 0)
            error = self.awg.AWGqueueWaveform(channel, waveformID, self.trigger_mode, 0, num_cycles, 0)
            if error < 0:
                raise KeysightException(f"AWG queue waveform error (channel {channel}): {error}")
        # start all AWGs simultaneously
        error = self.awg.AWGstartMultiple(awg_mask)
        if error < 0:
            raise KeysightException(f"AWG start error: {error}")

    def setup_single_channel(self, channel, num_cycles):
        self.awg.AWGflush(channel)
        wave = keysightSD1.SD_Wave()
        waveformID = channel
        error = wave.newFromArrayDouble(0, self.waveforms[self.channels[channel]])
        if error < 0:
            raise KeysightException(f"AWG create waveform error (channel {channel}): {error}")
        error = self.awg.waveformLoad(wave, waveformID, 0)
        if error < 0:
            raise KeysightException(f"AWG load waveform error (channel {channel}): {error}")
        # setting the delay in the AWGqueueWaveform method actually doesn't work properly when queuing multiple cycles
        #error = self.awg.AWGqueueWaveform(channel, waveformID, self.trigger_mode, self.channel_delays[channel]//10, num_cycles, 0)
        error = self.awg.AWGqueueWaveform(channel, waveformID, self.trigger_mode, 0, num_cycles, 0)
        if error < 0:
            raise KeysightException(f"AWG queue waveform error (channel {channel}): {error}")
        
    def launch_channels_start_only(self, awg_mask):
        # start all AWGs simultaneously
        error = self.awg.AWGstartMultiple(awg_mask)
        if error < 0:
            raise KeysightException(f"AWG start error: {error}")

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
            raise KeysightException(f"Couldn't connect to DAQ {self.model_name} in slot {self.slot}")
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
    
    def set_buffer_length(self, buffer_size):
        self.points_per_cycle = buffer_size
    
    def set_capture_cycles(self, capture_cycles):
        self.cycles = capture_cycles
    
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
        read_timeout = 100 # timeout in ms; 0 is infinite (DO NOT SET TO ZERO OR YOUR PROGRAM WILL HANG)
        for channel in self.channels.keys():
            data[self.channels[channel]] = self.daq.DAQread(channel, self.cycles*self.points_per_cycle, read_timeout)
        return data

    def stop(self):
        # just in case it's not stopped, stop the DAQ
        self.daq.DAQstopMultiple(sum(2**(c-1) for c in self.channels.keys()))
        time.sleep(0.1)
        self.daq.close()