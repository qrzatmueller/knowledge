'''
In this file I will attempt to investigate the effects of having a thickness of one half of the value
for the second half of the BW beam. What happens to the wavefield, what happens to the identification.
What happens if I only use the first half of measurement points, what about the second half, all?
'''

import os
import sys
import pandas as pd
import numpy as np
import olf
import BWoptimization as BW
import SolutionSurface as SS
import BeamWaves as beamW
import matplotlib.pyplot as plt

current_file_path = os.path.abspath(__file__)
mt_eq_path = os.path.abspath(os.path.join(current_file_path, '..', '..', '..', '..', '..'))
file_path = os.path.join(mt_eq_path, 'CodeAsterModels/YACSauto/YACS_BW_TH/autoRes/InHvsEta/resInHsmall_F0.txt')
file_path2 = os.path.join(mt_eq_path, 'CodeAsterModels/YACSauto/YACS_BW_TH/autoRes/InHvsEta/resInHsmall_F1.txt')
file_path3 = os.path.join(mt_eq_path, 'CodeAsterModels/YACSauto/YACS_BW_TH/autoRes/InHvsEta/resInHsmall_F2.txt')
#etas,cbs,Es = BW.IdentifyEta_YACS(file_path, "BW", "Plot", "HalfHalf")
case = "fullBW"
plot = "a"
remark = "test"
frequencies, responses = SS.openFEMresults(file_path)
frequencies2, responses2 = SS.openFEMresults(file_path2)
frequencies3, responses3 = SS.openFEMresults(file_path3)

#make Xpos for the first half, and the second half separate.
#is the Xpos for BW completely wrong??
if (case == "firstHalf"):
    start = 0.07
    step = 0.02
    x_size = 8
    stop = start + ((x_size-1) * step)
    Xpos = np.linspace(start, stop, x_size)
    beam = beamW.EBBeam(b=0.05) # for initial guess use steel beam
    #cut data to only the first 8 values
    responsesP = responses[:,:8]
if (case == "secondHalf"):
    start = 0.33
    step = 0.02
    x_size = 8
    stop = start + ((x_size-1) * step)
    Xpos = np.linspace(start, stop, x_size)
    beam = beamW.EBBeam(b=0.05) # for initial guess use steel beam

    responsesP = responses[:,-8:]
if (case == "BW"):
    #Xpos for usual BW, I have not retested this
    start = 0.070
    step = 0.020
    x_size = 21
    stop = start + ((x_size-1) * step)
    Xpos = np.linspace(start, stop, x_size)
    responsesP = responses
    
    beam = beamW.EBBeam(b=0.05) # for initial guess use steel beam
if (case == "OLF"):
    #Xpos for OLF, and guessing beam:
    start = 0.170
    step = 0.01 #it cannot work like this because the steps are not regelmäßig
    x_size = 12 #meaning one step but 2 datapoints
    stop = start + ((x_size-1) * step)
    Xpos = np.linspace(start, stop, x_size)
    Xpos2 = np.linspace(start+0.2, stop+0.2, x_size)
    Xpos = np.concatenate((Xpos, Xpos2), axis = 0)
    beam = beamW.EBBeam(b = 0.03) # for initial guess use steel beam
if (case == "fullBW"):
    start = 0.0
    step = 0.005
    stop = 0.48
    x_size = (stop-start)/step
    Xpos = np.linspace(start, stop, int(x_size)+1)
    responsesP = responses
    responsesP2 = responses2
    responsesP3 = responses3
    beam = beamW.EBBeam(b=0.05) # for initial guess use steel beam
    
n = 400 #which frequency to plot
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(221)

ax.plot(Xpos, responsesP[n,:].real, label = "eta01")
ax.plot(Xpos, responsesP2[n,:].real, label = "eta03")
ax.plot(Xpos, responsesP3[n,:].real, label = "eta05")
# Create reference plane at x = 0.27

ax.set_xlabel("Position (x)")
ax.set_ylabel("Real Part")
plt.legend()

ax = fig.add_subplot(222)

ax.plot(Xpos, responsesP[n,:].imag, label = "eta01")
ax.plot(Xpos, responsesP2[n,:].imag, label = "eta03")
ax.plot(Xpos, responsesP3[n,:].imag, label = "eta05")
# Create reference plane at x = 0.27

ax.set_xlabel("Position (x)")
ax.set_ylabel("Imag Part")
plt.legend()

ax2 = fig.add_subplot(223)
ax2.plot(Xpos, np.abs(responsesP[n,:]), label = "eta01")
ax2.plot(Xpos, np.abs(responsesP2[n,:]), label = "eta03")
ax2.plot(Xpos, np.abs(responsesP3[n,:]), label = "eta05")
ax.set_xlabel("Position (x)")
ax.set_ylabel("Magnitude")
plt.legend()
ax3 = fig.add_subplot(224)
ax3.plot(Xpos, np.angle(responsesP[n,:]), label = "eta01")
ax3.plot(Xpos, np.angle(responsesP2[n,:]), label = "eta03")
ax3.plot(Xpos, np.angle(responsesP3[n,:]), label = "eta05")
ax.set_xlabel("Position (x)")
ax.set_ylabel("Phase")
plt.legend()
plt.show()

fig2 = plt.figure(figsize=(8,6))
ax = fig2.add_subplot(111,projection = '3d')

ax.plot(Xpos, responsesP[n,:].real, responsesP[n,:].imag, label = "eta01")
ax.plot(Xpos, responsesP2[n,:].real, responsesP2[n,:].imag, label = "eta03")
ax.plot(Xpos, responsesP3[n,:].real, responsesP3[n,:].imag, label = "eta05")
# Create reference plane at x = 0.27
x_plane = np.full((10, 10), 0.27)  # A constant x-plane
y_plane = np.linspace(-0.01, 0.01, 10)  # Y-axis (real part range)
z_plane = np.linspace(-0.01, 0.01, 10)  # Z-axis (imaginary part range)
Y, Z = np.meshgrid(y_plane, z_plane)  # Create a grid
ax.plot_surface(x_plane, Y, Z, color='gray', alpha=0.3)  # Add semi-transparent plane
ax.set_xlabel("Position (x)")
ax.set_ylabel("Real Part")
ax.set_zlabel("Imag Part")
plt.legend()
plt.legend()
plt.show()
# etas = np.zeros(frequencies.shape)
# cbs = np.zeros(frequencies.shape)
# for i, freq in enumerate(frequencies):
#     meas = responsesP[i,:]
#     result = BW.Optimize(freq, Xpos, meas, BW.ObjectiveNormal)
#     k_ = result.x[0]+1j*result.x[1]
#     if (result.x[1]>0):
#         print("imaginary part is positive at frequency "+str(freq)+ ", check formulation")
#     E_ = beam.mu * freq*freq / beam.I /(pow(k_,4))
#     etas[i] = E_.imag / E_.real
#     omega = 2*np.pi*freq
#     cbs[i] = pow(E_.real*beam.I*omega*omega/beam.mu,0.25)
# if (plot != ""):
#     olf.save_scatter_plot(file_path,"BWID"+remark, frequencies, etas, plot)
