'''
The aim of this file is to take a simulation of a two layered beam, 
identify its loss factor and stiffness with the tools from BiegeWelle
verfahren and after that implement the formulas from Cremer to identify
eta2 and G2. This value can then be compared with the given values for 
the simulation.
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
file_path = os.path.join(mt_eq_path, 'CodeAsterModels/BeamVs3D/twoLayers_Files/resFiles/desp.txt') 
label1 = "EB identified"
'''
Identification with Euler Bernoulli
'''
case = "2Layers"
plot = "a"
remark = "test"
frequencies, responses = SS.openFEMresults(file_path)
firstColumn = 25 *4 #this eight is because I refined the mesh 400 elements now
lastColumn =75 * 4
cStep = 1 * 4
responsesR = responses[:,firstColumn:lastColumn:cStep]
# responsesR = responsesR[:,0:15]
elementSize = 0.0048 / 4
start =firstColumn * elementSize 
step = elementSize * cStep
x_size = int((lastColumn-firstColumn)/cStep)
stop = start + ((x_size-1) * step)
Xpos = np.linspace(start, stop, x_size)

# #plotting frequency
# n = 1
# frequency = frequencies[n]
# #plotting response (real imaginary, magnitude, phase)
# fig = plt.figure(figsize=(8,6))
# ax0 = fig.add_subplot(221)
# ax0.plot(Xpos, responsesR[n,:].real, label = "d = 10mm, h=1.5mm")
# # ax0.plot(Xpos, responsesP2[n,:].real, label = "d = 0, h = 1mm")
# # plt.axvline(0.27, color="black", linestyle="--", linewidth=0.8)
# # plt.axvline(0.28, color="black", linestyle="--", linewidth=0.8)
# #ax.plot(Xpos, responsesP3[n,:].real, label = "thck1")
# # Create reference plane at x = 0.27

# ax0.set_xlabel("Position (x)")
# ax0.set_ylabel("Real Part")
# plt.title(str(frequency)+ "Hz, eta=0.1")
# plt.legend()

# ax1 = fig.add_subplot(222)
# ax1.plot(Xpos, responsesR[n,:].imag, label = "d = 10mm, h=1.5mm")
# # ax1.plot(Xpos, responsesP2[n,:].imag, label = "d = 0, h = 1mm")
# #ax.plot(Xpos, responsesP3[n,:].imag, label = "thck1")
# # plt.axvline(0.27, color="black", linestyle="--", linewidth=0.8)
# # plt.axvline(0.28, color="black", linestyle="--", linewidth=0.8)
# # Create reference plane at x = 0.27
# ax1.set_xlabel("Position (x)")
# ax1.set_ylabel("Imag Part")
# plt.legend()


# ax2 = fig.add_subplot(223)
# ax2.plot(Xpos, np.abs(responsesR[n,:]), label = "d = 10mm, h=1.5mm")
# # ax2.plot(Xpos, np.abs(responsesP2[n,:]), label = "d = 0, h = 1mm")
# #ax2.plot(Xpos, np.abs(responsesP3[n,:]), label = "thck1")
# # plt.axvline(0.27, color="black", linestyle="--", linewidth=0.8)
# # plt.axvline(0.28, color="black", linestyle="--", linewidth=0.8)
# ax2.set_xlabel("Position (x)")
# ax2.set_ylabel("Magnitude")
# plt.legend()

# ax3 = fig.add_subplot(224)
# ax3.plot(Xpos, np.angle(responsesR[n,:]), label = "d = 10mm, h=1.5mm")
# # ax3.plot(Xpos, np.angle(responsesP2[n,:]), label = "d = 0, h = 1mm")
# #ax3.plot(Xpos, np.angle(responsesP3[n,:]), label = "thck1")
# # plt.axvline(0.27, color="black", linestyle="--", linewidth=0.8)
# # plt.axvline(0.28, color="black", linestyle="--", linewidth=0.8)
# ax3.set_xlabel("Position (x)")
# ax3.set_ylabel("Phase")
# plt.legend()
# plt.show()

#identify eta
etas = np.zeros(frequencies.shape)
cbs = np.zeros(frequencies.shape)
k_s = np.zeros(frequencies.shape, dtype = "complex")
for i, freq in enumerate(frequencies):
    meas = responsesR[i,:]
    result = BW.Optimize(freq, Xpos, meas, BW.ObjectiveNormal)
    k_ = result.x[0]+1j*result.x[1]
    k_s[i] = k_
    if (result.x[1]>0):
        print("imaginary part is positive at frequency "+str(freq)+ ", check formulation")
    # E_ = beam.mu * freq*freq / beam.I /(pow(k_,4))
    # etas[i] = E_.imag / E_.real
    k4 = k_**4
    etas[i] = -k4.imag / k4.real #todo I havent used this yet so I need to double check it
    omega = 2*np.pi*freq
    # cbs[i] = pow(E_.real*beam.I*omega*omega/beam.mu,0.25)
freqKs = np.zeros([frequencies.size,2],dtype = "complex")
freqKs[:,0] = frequencies
freqKs[:,1] = k_s
print("EB ident etas: ", etas)
# #plot identified eta
# fig4 = plt.figure(figsize=(8,6))
# ax4 = fig4.add_subplot(111)
# ax4.plot(frequencies,etas, label = label1)
# ax4.set_xlabel("Frequency [Hz]")
# ax4.set_ylabel("Loss Factor eta")
# plt.legend()
# plt.show()

'''
First, forward calculation then identification with cremer
'''
E1 = 206e9
d1 = 0.001
rho1 = 7850
E2 = 1e9
d2 = 0.0015
rho2 = 2600
eta2 = 1
system = cremer.EinfachBeläge(d1,E1,rho1,d2,E2,rho2,eta2)
etaTotsimp, etaTotex, etaTotchat = system.CalcEtaTot()
Bforward = system.CalcBTot()
print("forward calc etas: ",[etaTotsimp, etaTotex, etaTotchat])

EIeff, y = system.ChatStiffness()
BTotCalcd = system.CalcBTotEx()
BTot_Calcd = system.CalcBTot_Ex()
print("Bsimp, StiffExact, StiffComplx, StiffChat: ", [Bforward, BTotCalcd, BTot_Calcd, EIeff])
# print("d2/d1:", d2/d1)
# print("E2/E1", E2/E1)

# for i, freq in enumerate(frequencies):
#     i = 1

#identification:
system2 = cremer.EinfachBeläge(d1, E1, rho1, d2, 0, rho2, 0)
i = 0
etaMeas = etas[i]
print("etaMeas:", etaMeas)
freq = frequencies[i]
print("frequency", freq)
k_ = k_s[i] 
k = k_.real
om = freq * 2 * np.pi
m_ = system2.d1*system2.rho1 + system.d2*system2.rho2
Bmeas_ = om**2 * m_ / k_**4
print("Bmeas: ", Bmeas_.real)
E2id, eta2id, E2guess, eta2guess = system2.IdentifyParameters(Bmeas_.real,etaMeas)
print("Eta2 id: ", eta2id)
print("Eta2 guess: ", eta2guess)
print("E2 id: ", E2id)
print("E2 guess: ", E2guess)

system3 = cremer.EinfachBeläge(d1,E1,rho1,d2,E2id,rho2,eta2id)
EtaTot = system3.CalcEtaTot()
print("EtaTot", EtaTot)