# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 16:31:14 2024
This is the only bending wave method python script that works at the moment
17.01.2025 and it requires some manual adjustments to the data it works with
meaning number of data points to use. 
@author: prqrz
""" 
import numpy as np
import matplotlib.pyplot as plt
import BeamWaves as fW
import SolutionSurface as SS
import os


# freq = 50
#VMs = fW.fabricateVibration(eta=0.68122, cb=69.23847, freq=freq, Vp=coefs[0], Vpj=coefs[1], Vm=coefs[2], Vmj=coefs[3], plot="Identified from Meas")
# def fabricateVibration(eta=0.001, cb=65, freq=500, Vp=1, Vm=0, Vpj=0, Vmj=0, plot=""):
#VMs = VMs.flatten()
#min_cb, min_eta, coefs2 = SS.exploreSpace(VMs, freq, 0.001, (50,90), (0.6,0.81), 500,500, "Re-Identified Error Space")

current_file_path = os.path.abspath(__file__)
mt_eq_path = os.path.abspath(os.path.join(current_file_path, '..', '..', '..'))


file_path = os.path.join(mt_eq_path, 'CodeAsterModels/OLF/EBlong_Files/resF/res0p2.txt')

case = "OlfMeasurement" #to see if the BW procedure can identify from OLF data correctly. 

if (case == "BWsimulation"):
    frequencies, responses = SS.openFEMresults("input/EBres0p300.txt")
    num= 350
    freq = frequencies[num]
    meas = responses[num,:]
    print("freq = ", freq)

    #Define X position vector, this should coincide with the amount of measurements in the .txt
    start = 0.170
    step = 0.020
    x_size = 21
    stop = start + (x_size * step)
    Xpos = np.arange(start,stop,step)

    # plt.figure()
    # plt.plot(Xpos, np.real(meas), marker='o')  # Add marker for clarity if desired
    # plt.xlabel('X-axis Label')
    # plt.ylabel('Y-axis Label')
    # plt.title("real FEM")
    # plt.show()

    # plt.figure()
    # plt.plot(Xpos, np.imag(meas), marker='o')  # Add marker for clarity if desired
    # plt.xlabel('X-axis Label')
    # plt.ylabel('Y-axis Label')
    # plt.title("imag FEM")
    # plt.show()

    #min_cb, min_eta, coefs = SS.exploreSpace(meas, freq, 0.001,0.05, 7850, (5,150), (0.0,1.5), 500,500, Xpos, "plot")

    freqlist = np.arange(300,320)
    results = np.zeros((freqlist.size, 3))

    for i, freqindex in enumerate(freqlist):
        meas = responses[freqindex,:] #the two dots could be changed to match the desired data.
        freq = frequencies[freqindex]
        min_cb, min_eta, coefs = SS.exploreSpace(meas, freq, 0.001,0.05, 7850, (40,110), (0.2,0.4), 500,300, Xpos, "")
        results[i,0] = freq
        results[i,1] = min_cb
        results[i,2] = min_eta

    print(results)
    

if (case == "OlfMeasurement"):
    frequencies, responses = SS.openFEMresults(file_path)
    num= 330 
    freq = frequencies[num]
    meas = responses[num,:]
    print("freq = ", freq)

    
    #Define X position vector, this should coincide with the amount of measurements in the .txt
    start = 0.170
    step = 0.01 #it cannot work like this because the steps are not regelmäßig
    x_size = 12 #meaning one step but 2 datapoints
    stop = start + ((x_size-1) * step)
    Xpos = np.linspace(start, stop, x_size)
    Xpos2 = np.linspace(start+0.2, stop+0.2, x_size)
    Xpos = np.concatenate((Xpos, Xpos2), axis = 0)
    min_cb, min_eta, coefs, residuals = SS.exploreSpace(meas, freq, 0.001, 0.03, 7850, (50,200), (0.05,1.5), 200,200, Xpos, "plot freq 930Hz(100)")
    
    # freqlist = np.arange(50,80)
    # results = np.zeros((freqlist.size, 3))

    # for i, freqindex in enumerate(freqlist):
    #     meas = responses[freqindex,:] #the two dots could be changed to match the desired data.
    #     freq = frequencies[freqindex]
    #     min_cb, min_eta, coefs = SS.exploreSpace(meas, freq, 0.001,0.03, 7850, (50,150), (0.1,0.4), 500,300, Xpos, "")
    #     results[i,0] = freq
    #     results[i,1] = min_cb
    #     results[i,2] = min_eta

    # print(results)
    