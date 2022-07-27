# This sample program shows how Python uses the Keysight SD1 Programming Libraries.
# ----------
# Import required system components
import sys
import keysightSD1

# MODULE CONSTANTS
PRODUCT = "M3202A"
CHASSIS = 1

# change slot number to your value
SLOT = 5
CHANNELS = [1, 2, 3]
DELAYS = [0, 0, 0]
AMPLITUDE = 1.0

#waveform_data_list = [i/100 for i in range(100)]
waveform_data_list = [[0]*100 for channel in CHANNELS]
for channel in CHANNELS:
    if channel == 1 or channel == 2:
        waveform_data_list[channel-1] = [0, 0, 0.5, 1, 0.7, 0.2] + [0]*4
    else:
        waveform_data_list[channel-1] = [0.5, 1, 0.7, 0.2, 0] + [0]*5
    waveform_data_list[channel-1] = waveform_data_list[channel-1] * 5

# CREATE AND OPEN MODULE
module = keysightSD1.SD_AOU()
moduleID = module.openWithSlot(PRODUCT, CHASSIS, SLOT)
if moduleID < 0:
    print("Module open error:", moduleID)
else:
    print("Module opened:", moduleID)
    print("Module name:", module.getProductName())
    print("slot:", module.getSlot())
    print("Chassis:", module.getChassis())
    print()
    wave = keysightSD1.SD_Wave()
    for channel in CHANNELS:
        module.channelWaveShape(channel, keysightSD1.SD_Waveshapes.AOU_AWG)
        module.channelAmplitude(channel, AMPLITUDE)

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
        error = module.AWGfromArray(channel, 0, DELAYS[channel-1], 0, 0, 0, waveform_data_list[channel-1])
        if error < 0:
            print(f"AWG from array error (channel {channel}):", error)
        else:
            print("AWG from array started successfully")
    # exiting...
    input("Press any key to stop AWG")
    for channel in CHANNELS:
        module.AWGstop(channel)
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

