# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 14:31:56 2025

@author: prqrz
In this file, A results file of the OLF simulation is taken, noise is added
to the results. and the identification is done with both OLF and BW procedures
Additionally, the BW simulation for a material with the same properties is 
carried out and noise is also added to this and then the loss factor identified. 

The idea is to compare the robustness to noisy data of both procedures on the same
data set, and the BW with its appropriate simulation. 3 Data set types in total

Also the difference of phase noise to noise on the real and imaginary part is 
checked. 

All the resulting pngs are saved in the folder of the first "file_path", I cut 
the results from the file_path2 and moved them to the file_path path.

the long bar with BW identified has res0p2 while the short bar BW identified has res0p200
"""
import os
import sys
import pandas as pd
import numpy as np
import olf
import BWoptimization as BW
import matplotlib.pyplot as plt
import SolutionSurface as SS
from scipy.ndimage import maximum_filter1d, gaussian_filter1d


current_file_path = os.path.abspath(__file__)
mt_eq_path = os.path.abspath(os.path.join(current_file_path, '..', '..', '..', '..', '..'))
#file_path = os.path.join(mt_eq_path, 'CodeAsterModels/OLF/EBlong_Files/resF/res0p2.txt')
#



# h = 0.001
# b = 0.03
# rho = 7850
# pointDistance = 0.2
# eta3, cb3, E3, frequencies = olf.olf_improved(b, h, rho, file_path, pointDistance, windowSize=1, noise_coef= 0.01, noise_type="phase")
# olf.save_scatter_plot(file_path, x_vector= frequencies,y_vector= eta3, title=  "Loss Factor vs Freq (OLF) 1% Ph-Noise", remark = "_OLF_Noise01pct_Ph")
# eta3, cb3, E3, frequencies = olf.olf_improved(b, h, rho, file_path, pointDistance, windowSize=1, noise_coef = 0.01, noise_type="notphase")
# olf.save_scatter_plot(file_path, x_vector= frequencies,y_vector= eta3, title=  "Loss Factor vs Freq (OLF) 1% Noise", remark = "_OLF_Noise01pct_")

# eta3, cb3, E3, frequencies = olf.olf_improved(b, h, rho, file_path, pointDistance, windowSize=1, noise_coef= 0.05, noise_type="phase")
# olf.save_scatter_plot(file_path, x_vector= frequencies,y_vector= eta3, title=  "Loss Factor vs Freq (OLF) 5% Ph-Noise", remark = "_OLF_Noise05pct_Ph")
# eta3, cb3, E3, frequencies = olf.olf_improved(b, h, rho, file_path, pointDistance, windowSize=1, noise_coef = 0.05, noise_type="notphase")
# olf.save_scatter_plot(file_path, x_vector= frequencies,y_vector= eta3, title=  "Loss Factor vs Freq (OLF) 5%Noise", remark = "_OLF_Noise05pct_")

# etas,cbs,_ = BW.IdentifyEta_YACS_Noise(file_path,0.01, case = "OLF", plot = "identified Eta from 1% Noisy Simulation", remark= "_Noise01pct", noise_type="")
# etas,cbs,_ = BW.IdentifyEta_YACS_Noise(file_path,0.01, case = "OLF", plot = "identified Eta from 1% PH-Noisy Simulation", remark= "_Noise01pct_Ph")
# etas,cbs,_ = BW.IdentifyEta_YACS_Noise(file_path,0.05, case = "OLF", plot = "identified Eta from 5% Noisy Simulation", remark= "_Noise05pct", noise_type="")
# etas,cbs,_ = BW.IdentifyEta_YACS_Noise(file_path,0.05, case = "OLF", plot = "identified Eta from 5% PH-Noisy Simulation", remark= "_Noise05pct_Ph")

# file_path2 = os.path.join(mt_eq_path, 'CodeAsterModels/AnregungUnt/PoissonRatio_Files/resF/EBres0p200.txt')
# etas,cbs,_ = BW.IdentifyEta_YACS_Noise(file_path2,0.01, case = "BW", plot = "identified Eta from 1% Noisy Simulation", remark= "_Noise01pct", noise_type="")
# etas,cbs,_ = BW.IdentifyEta_YACS_Noise(file_path2,0.01, case = "BW", plot = "identified Eta from 1% PH-Noisy Simulation", remark= "_Noise01pct_Ph")
# etas,cbs,_ = BW.IdentifyEta_YACS_Noise(file_path2,0.05, case = "BW", plot = "identified Eta from 5% Noisy Simulation", remark= "_Noise05pct", noise_type="")
# etas,cbs,_ = BW.IdentifyEta_YACS_Noise(file_path2,0.05, case = "BW", plot = "identified Eta from 5% PH-Noisy Simulation", remark= "_Noise05pct_Ph")

'''
Section to investigate the relative error curve w.r.t increasing added error 
what are the frequencies used in this?
The maximum error added seems to be the same for high frequencies, but 10 times as high for low frequencies. 
For example, a 0.005 added error is also 0.005 loss factor error at high frequencies but 0.05 error for low ones!
 the randomness of the points is also observed, maybe trying with non random error could also give an 
insight

tried with noise type=max_noise to see if it would give a nice engulfing graph
so using this, the same pattern of greater problems at low frequencies appear, but
there is no super nice contour. the error is mostly positive at high frequencies, 
and negative at low frequencies. so maybe when there are many wavelengths in the measurement aperture
the error is positive, but having only one waelength that is most likely in the positive side, gives an artificially large
amplitude maybe? this insigh is not that usefull, I will leae it there

I lied. 
There are some frequencies at which the error is maximum, there are 54 hz, 70Hz, 136 hz, 249 hz, do they correspond
to eigenfrequencies? No they do not, the eigenfrequencies were 23, 63,124,132,206,265,307
using the sign formula showd that the frequencies are because of this, so lets let it gooo
'''
# file_path2 = os.path.join(mt_eq_path, 'CodeAsterModels/AnregungUnt/PoissonRatio_Files/resF/EBres0p200.txt')
# eta_ideal = 0.2
# etas,cbs,_,frequencies = BW.IdentifyEta_YACS_Noise(file_path2,0.02, case = "BW", plot = "identified Eta from 0.5% Noisy Simulation", remark= "_Noise0p5pct", noise_type="max_noise")

# loss_factor_error0 = np.abs((etas - eta_ideal) / eta_ideal)
# loss_factor_error1 = ((etas - eta_ideal) / eta_ideal)

# etas,cbs,_,_ = BW.IdentifyEta_YACS_Noise(file_path2,0.01, case = "BW", plot = "identified Eta from 1% Noisy Simulation", remark= "_Noise01pct", noise_type="")
# loss_factor_error1 = np.abs((etas - eta_ideal) / eta_ideal)
# etas,cbs,_,_ = BW.IdentifyEta_YACS_Noise(file_path2,0.02, case = "BW", plot = "identified Eta from 2% Noisy Simulation", remark= "_Noise02pct", noise_type="")
# loss_factor_error2 = np.abs((etas - eta_ideal) / eta_ideal)
# etas,cbs,_,_ = BW.IdentifyEta_YACS_Noise(file_path2,0.04, case = "BW", plot = "identified Eta from 4% Noisy Simulation", remark= "_Noise04pct", noise_type="")
# loss_factor_error3 = np.abs((etas - eta_ideal) / eta_ideal)

# # ----- Save output -----
# output = np.column_stack((frequencies, loss_factor_error0, loss_factor_error1, loss_factor_error2, loss_factor_error3))
# np.savetxt('Inv/NoiseInv/NoiseData.dat', output, header='NonDimLength  RelLossFactorError', fmt='%.6e')

# # ------ Retrieve previously calculated data data ------
# data = np.loadtxt('Inv/NoiseInv/NoiseData.dat')
# frequencies = data[:, 0]       # Hz
# loss_factor_error0 = data[:, 1]      # dimensionless
# loss_factor_error1 = data[:, 2] 
# loss_factor_error2 = data[:, 3] 
# loss_factor_error3 = data[:, 4] 

# #plot
# plt.figure(figsize = (8,5))
# plt.plot(frequencies, loss_factor_error0, label='Relative Error')
# plt.plot(frequencies, loss_factor_error1, label='Relative Error1')
# #plt.plot(frequencies, loss_factor_error2, label='Relative Error')
# #plt.plot(frequencies, loss_factor_error3, label='Relative Error')
# plt.axhline(0, color='gray', linestyle='--', linewidth=1)
# plt.xlabel('freq')
# plt.ylabel('Relative Loss Factor Error')
# plt.title('Loss Factor Error vs. Non-dimensional Length')
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.show()

'''
In this section I will explore to see if there is a pattern between the loss factor 
of a simulation and its sensitivity to error using the files with 0.1, 0.2 and 0.3 eta  and 1pct of error
it seems the lower the damping, the more sensitive to added error. I cannot think of a nice graph to put 
this all in a visible manner. 
'''
window_size = 50  #for moving maximum, then its smoothed with gaussian filter with sigma = 15 (10 to 30 smooth a lot 5 keeps true)

file_path2 = os.path.join(mt_eq_path, 'CodeAsterModels/AnregungUnt/PoissonRatio_Files/resF/EBres0p100.txt')
eta_ideal = 0.1
etas,cbs,_,frequencies = BW.IdentifyEta_YACS_Noise(file_path2,0.01, case = "BW", plot = "identified Eta from 01% Noisy Simulation", remark= "_Noise1pctEtap1", noise_type="")
loss_factor_error0 = np.abs((etas - eta_ideal) / eta_ideal)
loss_factor_error0 = maximum_filter1d(loss_factor_error0, size=window_size)
loss_factor_error0 = gaussian_filter1d(loss_factor_error0, sigma = 15)

file_path2 = os.path.join(mt_eq_path, 'CodeAsterModels/AnregungUnt/PoissonRatio_Files/resF/EBres0p200.txt')
eta_ideal = 0.2
etas,cbs,_,frequencies = BW.IdentifyEta_YACS_Noise(file_path2,0.01, case = "BW", plot = "identified Eta from 01% Noisy Simulation", remark= "_Noise1pctEtap2", noise_type="")
loss_factor_error1 = np.abs((etas - eta_ideal) / eta_ideal)
loss_factor_error1 = maximum_filter1d(loss_factor_error1, size=window_size)
loss_factor_error1 = gaussian_filter1d(loss_factor_error1, sigma = 15)

file_path2 = os.path.join(mt_eq_path, 'CodeAsterModels/AnregungUnt/PoissonRatio_Files/resF/EBres0p300.txt')
eta_ideal = 0.3
etas,cbs,_,frequencies = BW.IdentifyEta_YACS_Noise(file_path2,0.01, case = "BW", plot = "identified Eta from 01% Noisy Simulation", remark= "_Noise1pctEtap3", noise_type="")
loss_factor_error2 = np.abs((etas - eta_ideal) / eta_ideal)
loss_factor_error2 = maximum_filter1d(loss_factor_error2, size=window_size)
loss_factor_error2 = gaussian_filter1d(loss_factor_error2, sigma = 15)

file_path2 = os.path.join(mt_eq_path, 'CodeAsterModels/AnregungUnt/PoissonRatio_Files/resF/EBres0p800.txt')
eta_ideal = 0.8
etas,cbs,_,frequencies = BW.IdentifyEta_YACS_Noise(file_path2,0.01, case = "BW", plot = "identified Eta from 01% Noisy Simulation", remark= "_Noise1pctEtap5", noise_type="")
loss_factor_error3 = np.abs((etas - eta_ideal) / eta_ideal)
loss_factor_error3 = maximum_filter1d(loss_factor_error3, size=window_size)
loss_factor_error3 = gaussian_filter1d(loss_factor_error3, sigma = 15)
#plot
plt.figure(figsize = (8,5))
plt.plot(frequencies, loss_factor_error0, label='eta = 0.1')
plt.plot(frequencies, loss_factor_error1, label='eta = 0.2')
plt.plot(frequencies, loss_factor_error2, label='eta = 0.3')
plt.plot(frequencies, loss_factor_error3, label='eta = 0.8')
plt.axhline(0, color='gray', linestyle='--', linewidth=1)
plt.xlabel('freq')
plt.ylabel('Relative Loss Factor Error')
plt.title('Loss Factor Error vs. Non-dimensional Length')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# # ----- Save output -----
# output = np.column_stack((frequencies, loss_factor_error0, loss_factor_error1, loss_factor_error2, loss_factor_error3))
# np.savetxt('Inv/NoiseInv/NoiseDataEtaSmooth.dat', output, header='freq  0 1 2 3', fmt='%.6e')

'''
For plotting the real and iamginary parts
I cannot see something super clear here either with respect to the high and lows  :(
maybe they already had a smaller amplitude, maaybe
'''

# case = "BW"
# plot = "a"
# remark = "test"
# frequencies, responses = SS.openFEMresults(file_path2)
# responses2 = responses
# # frequencies2, responses2 = SS.openFEMresults(file_path2)

# if (case == "BW"):
        # #Xpos for usual BW, I have not retested this
        # start = 0.070 #it used to be 0.17, but that is completely wrong, although no real effect.  
        # step = 0.020
        # x_size = 21
        # stop = start + ((x_size-1) * step)
        # Xpos = np.linspace(start, stop, x_size)
        # #beam = BW.EBBeam(b=0.05) # for initial guess use steel beam

# #plotting frequency
# n = 97 #70 hz peak
# n2 = 60 #54Hz peak
# #n2 = 187 #136hz peak
# n3 = 150 #106hz node
# frequency = frequencies[n]
# print(frequency)
# print(frequencies[n2])
# print(frequencies[n3])
# #plotting response (real imaginary, magnitude, phase)
# fig = plt.figure(figsize=(8,6))
# ax0 = fig.add_subplot(221)
# ax0.plot(Xpos, responses[n,:].real, label = "n1")
# ax0.plot(Xpos, responses2[n2,:].real, label = "n2")
# ax0.plot(Xpos, responses2[n3,:].real, label = "n3")
# # plt.axvline(0.27, color="black", linestyle="--", linewidth=0.8)
# # plt.axvline(0.28, color="black", linestyle="--", linewidth=0.8)
# #ax.plot(Xpos, responsesP3[n,:].real, label = "thck1")
# # Create reference plane at x = 0.27

# ax0.set_xlabel("Position (x)")
# ax0.set_ylabel("Real Part")
# plt.title(str(frequency)+ "Hz, eta=0.1")
# plt.legend()

# ax1 = fig.add_subplot(222)
# ax1.plot(Xpos, responses[n,:].imag, label = "d = 10mm, h=1.5mm")
# ax1.plot(Xpos, responses2[n2,:].imag, label = "d = 0, h = 1mm")
# ax1.plot(Xpos, responses2[n3,:].imag, label = "d = 0, h = 1mm")
# #ax.plot(Xpos, responsesP3[n,:].imag, label = "thck1")
# # plt.axvline(0.27, color="black", linestyle="--", linewidth=0.8)
# # plt.axvline(0.28, color="black", linestyle="--", linewidth=0.8)
# # Create reference plane at x = 0.27
# ax1.set_xlabel("Position (x)")
# ax1.set_ylabel("Imag Part")
# plt.legend()


# ax2 = fig.add_subplot(223)
# ax2.plot(Xpos, np.abs(responses[n,:]), label = "d = 10mm, h=1.5mm")
# ax2.plot(Xpos, np.abs(responses2[n2,:]), label = "d = 0, h = 1mm")
# ax2.plot(Xpos, np.abs(responses2[n3,:]), label = "d = 0, h = 1mm")
# #ax2.plot(Xpos, np.abs(responsesP3[n,:]), label = "thck1")
# # plt.axvline(0.27, color="black", linestyle="--", linewidth=0.8)
# # plt.axvline(0.28, color="black", linestyle="--", linewidth=0.8)
# ax2.set_xlabel("Position (x)")
# ax2.set_ylabel("Magnitude")
# plt.legend()

# ax3 = fig.add_subplot(224)
# ax3.plot(Xpos, np.angle(responses[n,:]), label = "d = 10mm, h=1.5mm")
# ax3.plot(Xpos, np.angle(responses2[n2,:]), label = "d = 0, h = 1mm")
# ax3.plot(Xpos, np.angle(responses2[n3,:]), label = "d = 0, h = 1mm")
# #ax3.plot(Xpos, np.angle(responsesP3[n,:]), label = "thck1")
# # plt.axvline(0.27, color="black", linestyle="--", linewidth=0.8)
# # plt.axvline(0.28, color="black", linestyle="--", linewidth=0.8)
# ax3.set_xlabel("Position (x)")
# ax3.set_ylabel("Phase")
# plt.legend()
# plt.show()
