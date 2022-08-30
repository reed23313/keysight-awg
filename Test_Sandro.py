from gettext import find
import sys
import traceback
import time
import datetime
import random
import matplotlib.pyplot as pltslack
import matplotlib.colors as mcolors
import numpy as np
import pxi_modules
import scipy.io as sio
import scipy.signal as sigproc
from pulse_lib import *

vread = np.linspace(0, 100, 11)
vin = vread[:-1]


period = 10

print(vread)
print(vin)