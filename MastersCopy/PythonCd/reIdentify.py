# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 16:31:14 2024

@author: prqrz
""" 
import numpy as np
import matplotlib.pyplot as plt
import BeamWaves as fW #I changed this from the olf fabricateWaves.py hopefully it works 16.01.25
import SolutionSurface as SS
import os

current_file_path = os.path.abspath(__file__)
mt_eq_path = os.path.abspath(os.path.join(current_file_path, '..', '..', '..'))
#file_path = os.path.join(mt_eq_path, 'CodeAsterModels/OLF/EBlong_Files/resF/middlepoints.txt')
file_path = os.path.join(mt_eq_path, 'Python/PythonCd/input/res0p800_1mm_nu0.txt')

#frequencies, responses = SS.openFEMresults("res0p800_1mm_nu0.txt")
frequencies, responses = SS.openFEMresults(file_path)
num= 350
freq = frequencies[num]
meas = responses[num,:]

#Define X position vector
start = 0.070
step = 0.020
x_size = 21
stop = start + (x_size * step)
Xpos = np.arange(start,stop,step)

plt.figure()
plt.plot(Xpos, np.real(meas), marker='o')  # Add marker for clarity if desired
plt.xlabel('X-axis Label')
plt.ylabel('Y-axis Label')
plt.title("real FEM")
plt.show()

plt.figure()
plt.plot(Xpos, np.imag(meas), marker='o')  # Add marker for clarity if desired
plt.xlabel('X-axis Label')
plt.ylabel('Y-axis Label')
plt.title("imag FEM")
plt.show()

min_cb, min_eta,coefs = SS.exploreSpace(meas, freq, 0.001,0.05,7850, (50,90), (0.6,0.81), 500,500, "0p800")

# freq = 50
VMs = fW.fabricateVibration(eta=0.68122, cb=69.23847, freq=freq, Vp=coefs[0], Vpj=coefs[1], Vm=coefs[2], Vmj=coefs[3], plot="Identified from Meas")
# def fabricateVibration(eta=0.001, cb=65, freq=500, Vp=1, Vm=0, Vpj=0, Vmj=0, plot=""):
VMs = VMs.flatten()
min_cb, min_eta, coefs2 = SS.exploreSpace(VMs, freq, 0.001, (50,90), (0.6,0.81), 500,500, "Re-Identified Error Space")