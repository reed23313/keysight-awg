# This sample program shows how Python uses the Keysight SD1 Programming Libraries.
# ----------
# Import required system components
import sys
import keysightSD1

waveform_data_list = [i/100 for i in range(100)]

# MODULE CONSTANTS
PRODUCT = "M3202A"
CHASSIS = 1

# change slot number to your value
SLOT = 5
CHANNEL = 1
AMPLITUDE = 1.0

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
    module.channelWaveShape(CHANNEL, keysightSD1.SD_Waveshapes.AOU_AWG)
    module.channelAmplitude(CHANNEL, AMPLITUDE)

    # WAVEFORM FROM FILE
    WAVE_NUMBER = 0

    # create, open from file, load to module RAM and queue for execution 
    wave = keysightSD1.SD_Wave()

    # set path to your file here
    wave.newFromFile("C:/Users/Public/Documents/keysightSD1/Examples/Waveforms/Gaussian.csv")
    module.waveformLoad(wave, WAVE_NUMBER)
    module.AWGqueueWaveform(CHANNEL, WAVE_NUMBER, 0, 0, 0, 0)
    error = module.AWGstart(CHANNEL)
    if error < 0:
        print("AWG from file error:", error)
    else:
        print("AWG from file started successfully")
    input("Press any key to stop AWG")
    module.AWGflush(CHANNEL)
    module.AWGstop(CHANNEL)
    input("Press any key to start AWG from array")

    # WAVEFORM FROM ARRAY/LIST
    # This function is equivalent to create a waveform with new,
    # and then to call waveformLoad, AWGqueueWaveform and AWGstart
    error = module.AWGfromArray(CHANNEL, 0, 0, 0, 0, 0, waveform_data_list)
    if error < 0:
        print("AWG from array error:", error)
    else:
        print("AWG from array started successfully")
    # exiting...
    input("Press any key to stop AWG")
    module.AWGstop(CHANNEL)
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

