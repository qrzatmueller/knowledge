# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 13:11:28 2025
from an exel file from Sika reproduce their calculations of the loss factor. 
It evaluates all files in a folder.
@author: prqrz
"""
import os
import olf
import numpy as np
import matplotlib.pyplot as plt

def evaluateOlfMeasurement(file_path):
    if os.path.exists(file_path):
        print("File found! "+file_path)
    else:
        print("File not found. Please check the file path.")
        

    frequencies, complex_matrix = olf.read_measurement_data(file_path)
    frequencies = frequencies[-1591:]

    # Display the results
    print("Frequencies:", frequencies)
    print("Complex Data Matrix Shape:", complex_matrix.shape)

    # Example usage
    window_size = 21  # Adjust the window size as needed (e.g., 17th to 37th frequency)
    smoothed_matrix = olf.moving_average_complex_matrix(complex_matrix, window_size)
    smoothed_matrix = smoothed_matrix[-1591:,:]
    # Example usage
    average_vector, abs_average_vector, level_vector = olf.process_smoothed_data(smoothed_matrix)

    # Example usage
    phase_lag = olf.compute_phase_lag(average_vector)
    nopi = phase_lag/np.pi

    eta = olf.compute_eta(level_vector, nopi)

    coefficients = olf.fit_polynomial(frequencies, eta, 2, freq_range=(215,6000))


    evaluated_values = np.polyval(coefficients, frequencies)

    plt.figure(figsize=(12, 6))
    plt.figure
    
    # Scatter plot for the real part
    plt.subplot(1, 2, 2)
    plt.scatter(frequencies, eta,s = 1, c="blue", label="OLF_loss Factor")
    #plt.scatter(frequencies, evaluated_values2, c="red",s = 1, label="OLF_Loss Factor, Default Phase")
    plt.axhline(0, color="black", linestyle="--", linewidth=0.8)
    plt.title("Loss factor vs Frequency ("+file_path+")")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Loss Factor [-]")
    plt.legend()


    # Scatter plot for the real part
    plt.subplot(1, 2, 1)
    plt.scatter(frequencies, nopi,s = 1, c="blue", label="Total Phase")
    #plt.scatter(frequencies, fake_angle, c="red",s = 1, label="Default Total Phase")
    plt.axhline(0, color="black", linestyle="--", linewidth=0.8)
    plt.title("Total Phase vs Frequency")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Total Phase [-]/pi]")
    plt.legend()
    return nopi, eta, evaluated_values
# file_path1 = "exl/10grad12.xls"
# file_path2 = "exl/30grad.xls"
# file_path3 = "exl/40grad.xls"
# file_path4 = "exl/50grad.xls"
# file_path5 = "exl/60grad.xls"
# file_path6 = "exl/70grad.xls"
# nopi1, eta1, evEta1 = evaluateOlfMeasurement(file_path1)
# nopi2, eta2, evEta2 = evaluateOlfMeasurement(file_path2)
# nopi3, eta3, evEta3 = evaluateOlfMeasurement(file_path3)
# nopi4, eta4, evEta4 = evaluateOlfMeasurement(file_path4)
# nopi5, eta5, evEta5 = evaluateOlfMeasurement(file_path5)
# nopi6, eta6, evEta6 = evaluateOlfMeasurement(file_path6)

# Folder containing the files
folder_path = "exl"
def evall(folder_path):
    # Initialize an empty list to store the vectors
    nopis = []
    etas = []
    evEtas = []
    # Loop through each file in the folder
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        if os.path.isfile(file_path):  # Ensure it's a file
            nopi, eta, evEta = evaluateOlfMeasurement(file_path)
            nopis.append(nopi)
            etas.append(eta)
            evEtas.append(evEta)

    # Convert the list of vectors into a matrix
    #matrix = np.column_stack(vectors)
    nopis = np.column_stack(nopis)
    etas = np.column_stack(etas)
    evEtas = np.column_stack(evEtas)
    return nopis, etas, evEtas

nopis1, etas1, evEtas1 = evall("exl")
# nopis2, etas2, evEtas2 = evall("SD660")

