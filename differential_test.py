# ----------
# Python - Sample Application to set up the AWG to output an array that was created with numpy.
# ----------
# Import required system components
import sys
# ----------
# Append the system path to include the location of Keysight SD1 Programming Libraries
# then import the library
sys.path.append('C:\Program Files (x86)\Keysight\SD1\Libraries\Python')
import keysightSD1 as sdlib # Import Python SD1 library and AWG/Digitizer commands.
# ----------
# Specify values for variables
product = 'M3202A' # Product's model number.
chassis = 1
slot = 5
# Opening the M3202A module
awg = sdlib.SD_AOU()
awgId = awg.openWithSlotCompatibility(product, chassis, slot, 1)
print(awgId)
# Perform Waveform flush
error = awg.waveformFlush()
print(error)
# Setup channels for differential partner mode
error=awg.channelWaveShape(1, sdlib.SD_Waveshapes.AOU_AWG)
print(error)
error=awg.channelWaveShape(2, sdlib.SD_Waveshapes.AOU_PARTNER)
print(error)
error=awg.channelWaveShape(3, sdlib.SD_Waveshapes.AOU_AWG)
print(error)
error=awg.channelWaveShape(4, sdlib.SD_Waveshapes.AOU_PARTNER)
print(error)
# Set amplitude on the Differential Channel pairs
error=awg.channelAmplitude(1, 1)
print(error)
error=awg.channelAmplitude(2, -1)
print(error)
error=awg.channelAmplitude(3, 1)
print(error)
error=awg.channelAmplitude(4, -1)
print(error)
# Configure Waveform and Trigger Characteristics to play the waveform
zeroDelay = 0
infiniteIterations = 0;
preScaler = 0
wave = sdlib.SD_Wave()
# Loading the AWG waveform type
wfm0 = 0
gaussianFilepath = r"C:\Users\Public\Documents\Keysight\SD1\Examples\Waveforms\Gaussian.csv"
# "\Waveforms\Gaussian.csv"
error=wave.newFromFile(gaussianFilepath)
print(error)
error=awg.waveformLoad(wave, wfm0)
print(error)
# Setting Trigger modes for the Queued AWG Waveforms on each Channel
triggerMode = sdlib.SD_TriggerModes.SWHVITRIG
print("TriggerMode set to SWHVITRIG. ")
# Queuing AWG waveforms on each AWG Channel
error=awg.AWGqueueWaveform(1, wfm0, triggerMode, zeroDelay, infiniteIterations, preScaler)
print(error)
error=awg.AWGqueueWaveform(3, wfm0, triggerMode, zeroDelay, infiniteIterations, preScaler)
print(error)
# Setting mask and perfoming Start operation on multiple AWGs
AWGmask = 5 # Start and trigger channels 1 & 3. Channels 2 & 4 will follow suit.
error=awg.AWGstartMultiple(AWGmask)
print(error)
# Perfoming Trigger operation on multiple AWGs for playing Queued Waveforms
error = awg.AWGtriggerMultiple(AWGmask)
print(error)
# Perfoming Stop operation on multiple AWGs
input("Press any key to stop : ")
error=awg.AWGstopMultiple(AWGmask)
print(error)
# Closing module
error = awg.close()
print(error)