# This sample program shows how Python uses the Keysight SD1 Programming Libraries.
# ----------
# Import required system components
import sys

sys.path.append('C:\Program Files (x86)\Keysight\SD1\Libraries\Python')
import keysightSD1

# MODULE CONSTANTS
PRODUCT = "M3202A"
CHASSIS = 1

# change slot number to your value
SLOT = 5
CHANNELS = [1, 2, 3, 4]
#CHANNELS = [1, 2]
DELAYS = [0, 0, 0, 0] # delay in ns
#DELAYS = [0, 0] # delay in ns
AMPLITUDE = [0.1, -2.5, 0.5, 0.5]

WAVEFORM = [
    [0.5, 1, 0.7, 0.2], # best (low ripple pulse) we have so far: 2.3ns FWHM
    [0.75, 1, 0.3, 0.1], # this is pretty good, moderate ripple, 1.9ns FWHM
    [0.6, 1, 0.5, 0.1], # this is pretty good, low ripple, 2.1ns FWHM
    [0.7, 1, 0.4, 0.0]
]

#waveform_data_list = [i/100 for i in range(100)]
BUFFER_LENGTH = 1000
waveform_data_list = []
for channel in CHANNELS:
    # pad beginning of pulses
    # do fine delay adjustments
    # modulo 5 since delay argument to AWGfromArray is in multiples of 5 samples
    waveform_data_list.append([])
    for d in range(DELAYS[channel-1] % 5):
        waveform_data_list[channel-1].append(0)
    # optimal (?) fast pulse waveform with minimal ripple
    # (this seems to act like some sort of FIR filter to cancel the super ripply impulse response of analog front end in AWG)
    waveform_data_list[channel-1] += WAVEFORM[channel-1]
    # pad end of pulse with zeros
    waveform_data_list[channel-1] += [0]*(BUFFER_LENGTH - len(waveform_data_list[channel-1]))
    # store multiple pulses in buffer
    # better to read through a long buffer than have to recycle through it at a high frequency; the buffering in the AWG can't keep up
    if BUFFER_LENGTH < 50:
        waveform_data_list[channel-1] = waveform_data_list[channel-1] * 5

# CREATE AND OPEN MODULE
module = keysightSD1.SD_AOU()
moduleID = module.openWithSlot(PRODUCT, CHASSIS, SLOT)
if moduleID < 0:
    print("Module open error:", moduleID)
else:
    error = module.waveformFlush()
    if error < 0:
        print("failed to flush waveforms")
    print("Module opened:", moduleID)
    print("Module name:", module.getProductName())
    print("slot:", module.getSlot())
    print("Chassis:", module.getChassis())
    print()
    awg_mask = 0
    for n,channel in enumerate(CHANNELS):
        awg_mask += 2**(channel-1)
        module.channelWaveShape(channel, keysightSD1.SD_Waveshapes.AOU_AWG)
        module.channelAmplitude(channel, -AMPLITUDE[n])
        if channel == 2:
            module.channelOffset(channel, -1.5)
        # very important, behavior seems wonky if the AWG queue fills up
        module.AWGflush(channel)

        # WAVEFORM FROM ARRAY/LIST
        # This function is equivalent to create a waveform with new,
        # and then to call waveformLoad, AWGqueueWaveform and AWGstart
        # AWGfromArray(channel, triggerMode, startDelay, cycles, prescaler, waveformType, waveformDataA)
        # channel - AWG channel number (e.g. 1-4)
        # triggerMode - 0: AUTOTRIG (immediate, upon AWGStart or when previous waveform in the queue finishes)
        #               1: SWHVITRIG (software trigger, upon AWGtrigger) 
        #               5: SWHVITRIG_CYCLE (identical to SWHVITRIG, except a trigger is required to each cycle of the waveform)
        #               2: EXTTRIG (hardware trigger, set by externalSource in AWGtriggerExternalConfig)
        #               6: ExTTRIG_CYCLE
        # startDelay - delay in multiples of 5ns for M3202A, multiples of 10ns for M3201A, max values 2^16 - 1
        #               to get finer resolution, change the delay manually by changing the data stored in the buffer
        #module.AWGfromArray(channel, 0, DELAYS[channel-1]//5,0,0,0,waveform_data_list[channel-1])

        # instead of using AWGfromArray (which will start the channels at different times),
        # separately load the data into each channel's queue and then start all channels simultaneously)
        waveformID = channel
        wave = keysightSD1.SD_Wave()
        res = wave.newFromArrayDouble(0, waveform_data_list[channel-1])
        error = module.waveformLoad(wave, waveformID, 0)
        if error < 0:
            print(f"AWG load waveform error (channel {channel}):", error)
        else:
            print(f"AWG waveform loaded in channel {channel} successfully")
        error = module.AWGqueueWaveform(channel, waveformID, keysightSD1.SD_TriggerModes.AUTOTRIG, DELAYS[channel-1]//5, 0, 0)
        if error < 0:
            print(f"AWG queue waveform error (channel {channel}):", error)
        else:
            print(f"AWG waveform queued in channel {channel} successfully")
    # start all AWGs simultaneously
    error = module.AWGstartMultiple(awg_mask)
    if error < 0:
        print(f"AWG start error:", error)
    else:
        print("AWG started successfully")
    # exiting...
    input("Press any key to stop AWG")
    module.AWGstopMultiple(awg_mask)
    module.close()
    print()
    print("AOU closed")
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

