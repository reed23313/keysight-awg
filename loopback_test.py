# This sample program shows how Python uses the Keysight SD1 Programming Libraries.
# ----------
# Import required system components
import sys
import time
import matplotlib.pyplot as plt
import numpy as np

sys.path.append('C:\Program Files (x86)\Keysight\SD1\Libraries\Python')
import keysightSD1

def waitUntilPointsRead(module, DAQchannel, totalPoints, timeOut):
    t0 = time.time()
    timeElapsed = 0
    totalPointsRead = 0
    while (totalPointsRead < totalPoints) and (timeElapsed < timeOut):
        totalPointsRead = module.DAQcounterRead(DAQchannel)
        if(totalPointsRead < totalPoints):
            time.sleep(.01) #input argument of time.sleep is in seconds
        timeElapsed = time.time() - t0

CHASSIS = 1
# awg CONSTANTS
AWG = ["M3202A",5]
AWG_CHANNELS = [2]
AWG_DELAYS = [0] # delay in ns
AWG_AMPLITUDE = 1.0

AWG_WAVEFORM = [
    [0.5, 1, 0.7, 0.2] # best (low ripple pulse) we have so far: 2.3ns FWHM
    #[0.75, 1, 0.3, 0.1], # this is pretty good, moderate ripple, 1.9ns FWHM
    #[0.6, 1, 0.5, 0.1], # this is pretty good, low ripple, 2.1ns FWHM
    #[0.7, 1, 0.4, 0.0]
]

# DAQ constants
DAQ = ["M3102A",7]
DAQ_CHANNELS = [2]
DAQ_POINTS_PER_CYCLE = 5000
DAQ_CYCLES = 2
DAQ_TOTAL_POINTS = DAQ_POINTS_PER_CYCLE*DAQ_CYCLES
DAQ_TRIG_DELAY = 0
DAQ_FULL_SCALE = [1]

# set up AWG buffers
AWG_BUFFER_LENGTH = 1000
waveform_data_list = []
for channel in AWG_CHANNELS:
    # pad beginning of pulses
    # do fine delay adjustments
    # modulo 5 since delay argument to AWGfromArray is in multiples of 5 samples
    waveform_data_list.append([])
    for d in range(AWG_DELAYS[AWG_CHANNELS.index(channel)] % 5):
        waveform_data_list[AWG_CHANNELS.index(channel)].append(0)
    # optimal (?) fast pulse waveform with minimal ripple
    # (this seems to act like some sort of FIR filter to cancel the super ripply impulse response of analog front end in AWG)
    waveform_data_list[AWG_CHANNELS.index(channel)] += AWG_WAVEFORM[AWG_CHANNELS.index(channel)]
    # pad end of pulse with zeros
    waveform_data_list[AWG_CHANNELS.index(channel)] += [0]*(AWG_BUFFER_LENGTH - len(waveform_data_list[AWG_CHANNELS.index(channel)]))
    # store multiple pulses in buffer
    # better to read through a long buffer than have to recycle through it at a high frequency; the buffering in the AWG can't keep up
    if AWG_BUFFER_LENGTH < 50:
        waveform_data_list[AWG_CHANNELS.index(channel)] = waveform_data_list[AWG_CHANNELS.index(channel)] * 5

# CREATE AND OPEN awg
awg = keysightSD1.SD_AOU()
# open AWG
awgID = awg.openWithSlot(AWG[0], CHASSIS, AWG[1])
if awgID < 0:
    print("awg open error:", awgID)
else:
    print("awg opened:", awgID)
    print("awg name:", awg.getProductName())
    print("slot:", awg.getSlot())
    print("Chassis:", awg.getChassis())
    print()
    awg_mask = 0
    for channel in AWG_CHANNELS:
        awg_mask += 2**(channel-1)
        awg.channelWaveShape(channel, keysightSD1.SD_Waveshapes.AOU_AWG)
        awg.channelAmplitude(channel, AWG_AMPLITUDE)
        # very important, behavior seems wonky if the AWG queue fills up
        awg.AWGflush(channel)

        # instead of using AWGfromArray (which will start the channels at different times),
        # separately load the data into each channel's queue and then start all channels simultaneously)
        waveformID = channel
        wave = keysightSD1.SD_Wave()
        res = wave.newFromArrayDouble(0, waveform_data_list[AWG_CHANNELS.index(channel)])
        error = awg.waveformLoad(wave, waveformID, 0)
        if error < 0:
            print(f"AWG load waveform error (channel {channel}):", error)
        else:
            print(f"AWG waveform loaded in channel {channel} successfully")
        error = awg.AWGqueueWaveform(channel, waveformID, 0, AWG_DELAYS[AWG_CHANNELS.index(channel)]//5, 0, 0)
        if error < 0:
            print(f"AWG queue waveform error (channel {channel}):", error)
        else:
            print(f"AWG waveform queued in channel {channel} successfully")
    # start all AWGs simultaneously
    error = awg.AWGstartMultiple(awg_mask)
    if error < 0:
        print(f"AWG start error:", error)
    else:
        print("AWG started successfully")

# set up DAQ/Digitizer
digitizer = keysightSD1.SD_AIN()
digitizerID = digitizer.openWithSlot(DAQ[0], CHASSIS, DAQ[1])
data = np.zeros(shape=(len(DAQ_CHANNELS),DAQ_CYCLES,DAQ_POINTS_PER_CYCLE))
if digitizerID < 0:
    print("digitizer open error:", digitizerID)
else:
    print("digitizer opened:", digitizerID)
    print("digitizer name:", digitizer.getProductName())
    print("slot:", digitizer.getSlot())
    print("Chassis:", digitizer.getChassis())
    print()
    daq_mask = 0
    for channel in DAQ_CHANNELS:
        daq_mask += 2**(channel-1)
        digitizer.channelInputConfig(channel, DAQ_FULL_SCALE[DAQ_CHANNELS.index(channel)], keysightSD1.AIN_Impedance.AIN_IMPEDANCE_50, keysightSD1.AIN_Coupling.AIN_COUPLING_DC)
        digitizer.DAQconfig(channel, DAQ_POINTS_PER_CYCLE, DAQ_CYCLES, 0, keysightSD1.SD_TriggerModes.SWHVITRIG)
        #digitizer.DAQconfig(channel, DAQ_POINTS_PER_CYCLE, DAQ_CYCLES, 0, keysightSD1.SD_TriggerModes.ANALOGTRIG)
        #digitizer.channelTriggerConfig(channel, keysightSD1.SD_AIN_TriggerMode.RISING_EDGE, 0.1)
    # measure all channels at once  
    digitizer.DAQflushMultiple(daq_mask)
    digitizer.DAQstartMultiple(daq_mask)
    for cycle in range(DAQ_CYCLES):
        digitizer.DAQtriggerMultiple(daq_mask)
        for channel in DAQ_CHANNELS:
            timeout = 0.1 # 100ms
            waitUntilPointsRead(digitizer, channel, DAQ_TOTAL_POINTS, timeout)
            print(f"read {digitizer.DAQcounterRead(channel)} point(s) from channel {channel}")
            read_timeout = 100 # 0 is inifinite
            data[DAQ_CHANNELS.index(channel), cycle] = digitizer.DAQread(channel, DAQ_TOTAL_POINTS, read_timeout)
    digitizer.DAQstopMultiple(daq_mask)
    print("plotting results")
    plt.figure()
    for channel in DAQ_CHANNELS:
        for cycle in range(DAQ_CYCLES):
            plt.plot(data[DAQ_CHANNELS.index(channel), cycle], label = f'ch{channel} - {cycle}')
    plt.legend()
    plt.show()

if awgID > 0:
    # exiting...
    input("Press any key to stop AWG")
    for channel in AWG_CHANNELS:
        awg.AWGstop(channel)
    awg.close()
    print()
    print("AOU closed")
if digitizerID > 0:
    input("Press any key to stop DAQ")
    digitizer.close()
    print()
    print("AIN closed")
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