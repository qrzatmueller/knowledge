# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 08:34:01 2024
this file uses the functions defined in olf.py to identify a loss factor given 
fem results from code aster
@author: prqrz
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import olf 
import SolutionSurface as SS

current_file_path = os.path.abspath(__file__)
mt_eq_path = os.path.abspath(os.path.join(current_file_path, '..', '..', '..'))
#file_path = os.path.join(mt_eq_path, 'CodeAsterModels/OLF/EBlong_Files/resF/res0p2.txt')
file_path = os.path.join(mt_eq_path, 'CodeAsterModels/OLF/autoRes/res0.txt')

frequencies, responses = SS.openFEMresults(file_path)
# n = responses.shape[1]/2
# n=int(n)
# first12 = responses[:,:n]
# lastqw = responses[:,-1*n:]
# complex_matrix = responses[:,-1*n:]/responses[:,:n]
# test = complex_matrix[:,0]-responses[:,n]/responses[:,0]

# window_size =1
# smoothed_matrix = olf.moving_average_complex_matrix(complex_matrix, window_size)

# average_vector, abs_average_vector, level_vector = olf.process_smoothed_data(smoothed_matrix)


# phase_lag = olf.compute_phase_lag(average_vector)
# nopi = phase_lag/np.pi

# eta = olf.compute_eta(level_vector, nopi)

# coefficients = olf.fit_polynomial(frequencies, eta, 3, freq_range=(0,6000))
# evaluated_values = np.polyval(coefficients, frequencies)


# #file_path2 = r"/home/pqrz/Documents/MasterThesis/mt_eq/CodeAsterModels/OLF/EBlong_Files/resF/res.txt"
# # file_path2 = os.path.join(mt_eq_path, 'CodeAsterModels/OLF/EBlong_Files/resF/res1p0.txt')
# # frequencies2, responses2 = SS.openFEMresults(file_path2)
# # n = responses2.shape[1]/2
# # n=int(n)
# # first12 = responses2[:,:n]
# # lastqw = responses2[:,-1*n:]
# # complex_matrix2 = responses2[:,-1*n:]/responses2[:,:n]


# # window_size2 = 1  # Adjust the window size as needed (e.g., 17th to 37th frequency)
# # smoothed_matrix2 = olf.moving_average_complex_matrix(complex_matrix2, window_size2)
# # # Example usage
# # average_vector2, abs_average_vector2, level_vector2 = olf.process_smoothed_data(smoothed_matrix2)

# # # Example usage
# # phase_lag2 = olf.compute_phase_lag(average_vector2)
# # nopi2 = phase_lag2/np.pi

# # eta2 = olf.compute_eta(level_vector2, nopi2)

h = 0.001
b = 0.03
rho = 7850
pointDistance = 0.2
eta3, cb3, E3, freqs = olf.olf_improved(b, h, rho, file_path, pointDistance, windowSize=5, noise_coef= 0, noise_type="")
olf.save_scatter_plot(file_path, x_vector= frequencies,y_vector= eta3, title=  "Loss Factor vs Freq (OLF)", remark = "_OLF_PhNoise01pct_")
# eta3, cb3, E3, freqs = olf.olf_improved(b, h, rho, file_path, pointDistance, windowSize=1, noise_coef= 0.01, noise_type="")
# olf.save_scatter_plot(file_path, x_vector= frequencies,y_vector= eta3, title=  "Loss Factor vs Freq (OLF)", remark = "_OLF_PhNoise01pct_")
# eta3, cb3, E3, freqs = olf.olf_improved(b, h, rho, file_path, pointDistance, windowSize=1, noise_coef = 0.01, noise_type="phase")
# olf.save_scatter_plot(file_path, x_vector= frequencies,y_vector= eta3, title=  "Loss Factor vs Freq (OLF)", remark = "_OLF_Noise01pct_")




aaa= plt.figure(figsize=(12, 6))
#plt.subplot(1, 2, 1)
#plt.scatter(frequencies, eta,s = 1, c="blue", label="FEM_OLF_loss Factor ")
#plt.scatter(frequencies2, eta2,s = 1, c="red", label="FEM_OLFlong_loss Factor ")
plt.scatter(frequencies, eta3,s = 1, c="black", label="improved ")
plt.axhline(0, color="black", linestyle="--", linewidth=0.8)
plt.title("Loss Factor vs Frequency")
plt.xlabel("Frequency")
plt.ylabel("Loss Factor")
plt.legend()
plt.show()



# aaa2 = plt.figure(figsize = (12,6))
# xPos = np.linspace(170,280,12)
# xPos = np.append(xPos,np.linspace(370,480,12))
# plt.subplot(1, 2, 1)
# plt.plot(xPos, np.abs(responses[0,:]), c="blue", label="response1 real ")
# plt.plot(xPos, np.abs(responses[1,:]), c="red", label="response2 real ")
# plt.axhline(0, color="black", linestyle="--", linewidth=0.8)
# plt.title("Loss Factor vs Frequency")
# plt.xlabel("Frequency")
# plt.ylabel("Loss Factor")
# plt.legend()
# plt.show()


# plt.subplot(1, 2, 2)
# plt.plot(xPos, np.angle(responses[0,:]), c="blue", label="response1 imag ")
# plt.plot(xPos, np.angle(responses[1,:]), c="red", label="response2 imag ")
# plt.axhline(0, color="black", linestyle="--", linewidth=0.8)
# plt.title("Loss Factor vs Frequency")
# plt.xlabel("Frequency")
# plt.ylabel("Loss Factor")
# plt.legend()
# plt.show()

# aaa3 = plt.figure(figsize = (12,6))

# level_matrix = 20 * np.log10(np.abs(complex_matrix))
# level_matrix2 = 20 * np.log10(np.abs(complex_matrix2))
# # Scatter plot for the real part
# plt.subplot(1, 2, 1)
# plt.plot(xPos[:12], level_matrix[0,:], c="blue", label="response1 real ")
# plt.plot(xPos[:12], level_matrix[100,:], c="red", label="response2 real ")
# plt.axhline(0, color="black", linestyle="--", linewidth=0.8)
# plt.title("Level vs Position")
# plt.xlabel("Position")
# plt.ylabel("Velocity Level dB")
# plt.legend()
# plt.show()


# plt.subplot(1, 2, 2)
# plt.plot(frequencies, level_vector, c="blue", label="level_vector")
# plt.plot(xPos[:12], np.angle(complex_matrix[100,:]), c="red", label="response2 imag ")
# plt.axhline(0, color="black", linestyle="--", linewidth=0.8)
# plt.title("Loss Factor vs Frequency")
# plt.xlabel("Frequency")
# plt.ylabel("Loss Factor")
# plt.legend()
# plt.show()

# aaa4 = plt.figure(figsize = (12,6))
# plt.subplot(1, 2, 1)
# plt.imshow(level_matrix, cmap='viridis', interpolation='nearest', aspect = 'auto')

# # Add a color bar to show the scale
# plt.colorbar(label="Values")

# # Add labels and title
# plt.xlabel("Column Index")
# plt.ylabel("Row Index")
# plt.title("Matrix as Color Plot")

# # Show the plot
# plt.show()
# plt.subplot(1, 2, 2)
# plt.imshow(level_matrix2, cmap='viridis', interpolation='nearest', aspect = 'auto')

# # Add a color bar to show the scale
# plt.colorbar(label="Values")

# # Add labels and title
# plt.xlabel("Column Index")
# plt.ylabel("Row Index")
# plt.title("Matrix as Color Plot")

# # Show the plot
# plt.show()