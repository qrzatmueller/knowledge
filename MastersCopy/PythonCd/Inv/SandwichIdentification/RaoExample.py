# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 16:09:07 2025

@author: prqrz

In this document I will try to identify sandwiches first using the EB assumption
used at MBBM and later with the 6th order developed identification.
you left at trying to make it work lol ideas:
try only with numpy stuff, once you are already reducing the matrix. 
if its true that all roots can be defined in terms of the other 2, can you change
the problem? probably not if for that description you already need the other parameters 
to be identified right?
If this in general does not work, would it work if from the euler bernoulli you get the
general damping and then you optimize for the quantities G and eta2, so only guessing 3 parameters???

question 2 would it work if I had more waves within one length? why are they so few?
"""
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
import RAOExacttt as rao

current_file_path = os.path.abspath(__file__)
mt_eq_path = os.path.abspath(os.path.join(current_file_path, '..', '..', '..', '..', '..'))
'''
case: Rao example 7.1. 
'''
# file_path = os.path.join(mt_eq_path, 'CodeAsterModels/BeamVs3D/RAOexampleLong_Files/resFiles/RAOlong_vites.txt') #raoLongCase
file_path = os.path.join(mt_eq_path, 'CodeAsterModels/BeamVs3D/Silka625_BW_Files/resFiles/desp.txt') #Sika635 case
label1 = "EB identified"

'''
Identification with Euler Bernoulli
'''
case = "shortBar"
plot = "a"
remark = "test"
refine = 4
frequencies, responses = SS.openFEMresults(file_path)
firstColumn = 25 *refine
lastColumn =75 *refine
cStep = 1 *refine
responsesR = responses[:,firstColumn:lastColumn:cStep]
# responsesR = responsesR[:,0:15]
elementSize = 0.0048 / refine
start =firstColumn * elementSize 
step = elementSize * cStep
x_size = int((lastColumn-firstColumn)/cStep)
stop = start + ((x_size-1) * step)
Xpos = np.linspace(start, stop, x_size)

# #calculation frequency
# n = 19 #n=8 for 100hz, n=159 for 1620
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

# #identify eta
# etas = np.zeros(frequencies.shape)
# cbs = np.zeros(frequencies.shape)
# k_s = np.zeros(frequencies.shape, dtype = "complex")
# for i, freq in enumerate(frequencies):
#     meas = responsesR[i,:]
#     result = BW.Optimize(freq, Xpos, meas, BW.ObjectiveNormal)
#     k_ = result.x[0]+1j*result.x[1]
#     k_s[i] = k_
#     if (result.x[1]>0):
#         print("imaginary part is positive at frequency "+str(freq)+ ", check formulation")
#     # E_ = beam.mu * freq*freq / beam.I /(pow(k_,4))
#     # etas[i] = E_.imag / E_.real
#     k4 = k_**4
#     etas[i] = -k4.imag / k4.real #todo I havent used this yet so I need to double check it
#     omega = 2*np.pi*freq
#     # cbs[i] = pow(E_.real*beam.I*omega*omega/beam.mu,0.25)
# freqKs = np.zeros([frequencies.size,2],dtype = "complex")
# freqKs[:,0] = frequencies
# freqKs[:,1] = k_s
# #plot identified eta
# fig4 = plt.figure(figsize=(8,6))
# ax4 = fig4.add_subplot(111)
# ax4.plot(frequencies,etas, label = label1)
# ax4.set_xlabel("Frequency [Hz]")
# ax4.set_ylabel("Loss Factor eta")
# plt.legend()
# plt.show()

'''
Identification mit Taranto
'''

'''
TrantoParameters
'''
# L = 0.48 #to match the length of a bw beam
# rho1 = 7850
# rho3 = rho1
# E1 = 206e9
# E3 = E1

# rho2 = 2600
# eta2 = 1 #originally 0.1
# G2 = 0.98e10

# b = 0.05 #this is not given in the paper. check if it affects anything. 
# h1 = 0.001 / 2 #todo changed the heights to check something
# h2 = 0.005 /2 #originally 0.005
# h3 = 0.001 / 2
'''
Sika Parameters
'''
L = 0.48 #to match the length of a bw beam
rho1 = 2700
E1 = 69e9

rho3 = 7850
E3 = 206e9

rho2 = 2600
eta2 = 1.5 #should 1.5 #originally 0.1 identified 0.248212 total 0.51395
G2 = 0.98e10 #0.98e10 #2000000000 #should 0.98e10 #identified 1.0853e8

b = 1 #this is not given in the paper. check if it affects anything. 
h1 = 0.00015/2 #originally 0.00015 / 2 #todo changed the heights to check something
h2 = 0.0015 /2 #originally 0.005
h3 = 0.001 / 2

#Create Sandwich object
example = rao.SandwichBeam(h1,E1, rho1, h3, E3, rho3, h2, G2, eta2, rho2, b, L)
XposN = Xpos/example.L #normalization of Xpos
example.ForwardCalc(file_path, "Custom", Xposi = XposN, resp = [firstColumn, lastColumn, cStep])
example.IdentifyG_(file_path, "Custom", Xposi = XposN, resp = [firstColumn, lastColumn, cStep])
# G2s = np.logspace(7.0,12.0,5)
# for i, Gs in enumerate(G2s): 
#     example = rao.SandwichBeam(h1,E1, rho1, h3, E3, rho3, h2, Gs, eta2, rho2, b, L)
#     XposN = Xpos/example.L #normalization of Xpos
#     example.ForwardCalc(file_path, "Custom", Xposi = XposN, resp = [firstColumn, lastColumn, cStep])

# ETAS = np.linspace(0.5,2,3)
# for i, etas in enumerate(ETAS): 
#     example = rao.SandwichBeam(h1,E1, rho1, h3, E3, rho3, h2, G2, etas, rho2, b, L)
#     XposN = Xpos/example.L #normalization of Xpos
#     example.ForwardWave(file_path, "Custom", Xposi = XposN, resp = [firstColumn, lastColumn, cStep])