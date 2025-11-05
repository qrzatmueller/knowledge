'''
With this script I will attempt to produce data that can be nicely plotted in tex
The objective is to do simple simulaitons with defined eta values and to identify those 
eta values

Simulation used: 
the yacs auto bw2 simulation, it accepts E, eta, and thickness. 
To just do the identification no script was needed other than the BW_aster_yacs.py, 
maybe I copy it here. 

what is here bellow is just for visualization


'''

import os
import sys
import pandas as pd
import numpy as np
import olf
import BWoptimization as BW
import SolutionSurface as SS
import BeamWaves as beamW
import olf
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import EinfachBelaege as cremer

current_file_path = os.path.abspath(__file__)
mt_eq_path = os.path.abspath(os.path.join(current_file_path, '..', '..', '..', '..', '..'))
'''
case: 1mm steel with 1.5mm of viscous material 
'''
file_path = os.path.join(mt_eq_path, 'CodeAsterModels/YACSauto/YACS_BW/autoRes/res3.txt') 
file_path2 = os.path.join(mt_eq_path, 'CodeAsterModels/YACSauto/YACS_BW/autoRes/res0.txt') 
label1 = "EB identified"
'''
Identification with Euler Bernoulli
'''
case = "BW"
plot = "a"
remark = "test"
frequencies, responses = SS.openFEMresults(file_path)
frequencies2, responses2 = SS.openFEMresults(file_path2)

if (case == "BW"):
        #Xpos for usual BW, I have not retested this
        start = 0.070 #it used to be 0.17, but that is completely wrong, although no real effect.  
        step = 0.020
        x_size = 21
        stop = start + ((x_size-1) * step)
        Xpos = np.linspace(start, stop, x_size)
        #beam = BW.EBBeam(b=0.05) # for initial guess use steel beam

#plotting frequency
n = frequencies.size-1
frequency = frequencies[n]
#plotting response (real imaginary, magnitude, phase)
fig = plt.figure(figsize=(8,6))
ax0 = fig.add_subplot(221)
ax0.plot(Xpos, responses[n,:].real, label = "d = 10mm, h=1.5mm")
ax0.plot(Xpos, responses2[n,:].real, label = "d = 0, h = 1mm")
# plt.axvline(0.27, color="black", linestyle="--", linewidth=0.8)
# plt.axvline(0.28, color="black", linestyle="--", linewidth=0.8)
#ax.plot(Xpos, responsesP3[n,:].real, label = "thck1")
# Create reference plane at x = 0.27

ax0.set_xlabel("Position (x)")
ax0.set_ylabel("Real Part")
plt.title(str(frequency)+ "Hz, eta=0.1")
plt.legend()

ax1 = fig.add_subplot(222)
ax1.plot(Xpos, responses[n,:].imag, label = "d = 10mm, h=1.5mm")
ax1.plot(Xpos, responses2[n,:].imag, label = "d = 0, h = 1mm")
#ax.plot(Xpos, responsesP3[n,:].imag, label = "thck1")
# plt.axvline(0.27, color="black", linestyle="--", linewidth=0.8)
# plt.axvline(0.28, color="black", linestyle="--", linewidth=0.8)
# Create reference plane at x = 0.27
ax1.set_xlabel("Position (x)")
ax1.set_ylabel("Imag Part")
plt.legend()


ax2 = fig.add_subplot(223)
ax2.plot(Xpos, np.abs(responses[n,:]), label = "d = 10mm, h=1.5mm")
ax2.plot(Xpos, np.abs(responses2[n,:]), label = "d = 0, h = 1mm")
#ax2.plot(Xpos, np.abs(responsesP3[n,:]), label = "thck1")
# plt.axvline(0.27, color="black", linestyle="--", linewidth=0.8)
# plt.axvline(0.28, color="black", linestyle="--", linewidth=0.8)
ax2.set_xlabel("Position (x)")
ax2.set_ylabel("Magnitude")
plt.legend()

ax3 = fig.add_subplot(224)
ax3.plot(Xpos, np.angle(responses[n,:]), label = "d = 10mm, h=1.5mm")
ax3.plot(Xpos, np.angle(responses2[n,:]), label = "d = 0, h = 1mm")
#ax3.plot(Xpos, np.angle(responsesP3[n,:]), label = "thck1")
# plt.axvline(0.27, color="black", linestyle="--", linewidth=0.8)
# plt.axvline(0.28, color="black", linestyle="--", linewidth=0.8)
ax3.set_xlabel("Position (x)")
ax3.set_ylabel("Phase")
plt.legend()
plt.show()
