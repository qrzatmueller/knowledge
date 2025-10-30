# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 08:34:01 2024
functions to reproduce the calculations done by Sika in their exel file. 
@author: prqrz
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import SolutionSurface as SS
import BWoptimization as BW

def read_measurement_data(file_path):
    # Load the Excel file
    excel_data = pd.ExcelFile(file_path)
    
    # Initialize lists for frequencies and measurement pairs
    frequencies = None
    complex_data = []
    # fake_angle = []
    # af = pd.read_excel(file_path, sheet_name="1", skiprows=5, header=None)
    # fake_angle = af.iloc[:,3]

    # Iterate through each sheet
    for sheet_name in excel_data.sheet_names:
        # Read the sheet, skipping the first rows with headers and units
        df = pd.read_excel(file_path, sheet_name=sheet_name, skiprows=5, header=None)
        
        # Extract columns by position: A = index 0, B = index 1, C = index 2
        freq_col = df.iloc[:, 0]  # First column (A) is frequencies
        real_col = df.iloc[:, 1]  # Second column (B) is the real part
        imag_col = df.iloc[:, 2]  # Third column (C) is the imaginary part
        
        # Extract data
        if frequencies is None:
            frequencies = freq_col.to_numpy()
        else:
            assert np.allclose(frequencies, freq_col.to_numpy()), "Frequencies do not match across sheets!"
        
        # Create complex numbers
        complex_numbers = real_col.to_numpy() + 1j * imag_col.to_numpy()
        
        # Append to the list
        complex_data.append(complex_numbers)
    
    # Combine all measurement pairs into a matrix
    complex_matrix = np.column_stack(complex_data)
    
    return frequencies, complex_matrix

def moving_average_complex_matrix(complex_matrix, window_size):
    """
    Applies a moving average to smooth the complex matrix.
    
    Parameters:
        complex_matrix (ndarray): Original matrix of complex values (frequencies x measurement pairs).
        window_size (int): Size of the moving average window (must be odd).
    
    Returns:
        ndarray: Smoothed complex matrix with the same shape as the input.
    """
    # Ensure the window size is odd
    assert window_size % 2 == 1, "Window size must be odd."
    
    # Number of rows (frequencies) and columns (measurement pairs)
    num_frequencies, num_pairs = complex_matrix.shape
    
    # Half window size for calculating bounds
    half_window = window_size // 2
    
    # Initialize the smoothed matrix
    smoothed_matrix = np.zeros_like(complex_matrix, dtype=complex)
    
    # Apply moving average for each frequency
    for i in range(num_frequencies):
        # Determine the range for the sliding window
        start = max(0, i - half_window)
        end = min(num_frequencies, i + half_window + 1)
        
        # Compute the average over the window
        smoothed_matrix[i, :] = np.mean(complex_matrix[start:end, :], axis=0)
    
    return smoothed_matrix
def process_smoothed_data(smoothed_matrix):
    """
    Processes the smoothed complex matrix step by step.

    Steps:
    1. Average over all measurement pairs (columns) for each frequency.
    2. Take the absolute value of each average.
    3. Compute the level (20log10 of the absolute value).
    
    Parameters:
        smoothed_matrix (ndarray): Smoothed complex matrix (frequencies x measurement pairs).
    
    Returns:
        Tuple of ndarrays: (average_vector, abs_average_vector, level_vector)
    """
    # Step 1: Average over all columns (measurement pairs) for each frequency
    average_vector = np.mean(smoothed_matrix, axis=1)  # Shape: (frequencies,)
    
    # Step 2: Take the absolute value of each average
    abs_average_vector = np.abs(average_vector)
    
    # Step 3: Compute the level (20log10 of the absolute value)
    level_vector = 20 * np.log10(abs_average_vector)
    
    return average_vector, abs_average_vector, level_vector

def compute_phase_lag(smoothed_matrix):
    """
    Computes the total phase lag between the second and first measurement points.
    
    Parameters:
        smoothed_matrix (ndarray): Smoothed complex matrix (frequencies x measurement pairs).
    
    Returns:
        ndarray: Vector of total phase lags for each frequency.
    """
    # Step 1: Calculate the phase (angle) for each measurement pair
    phase_matrix = np.angle(smoothed_matrix)  # Shape: (frequencies x measurement pairs)
    
    # Step 2: Unwrap the phase to avoid jumps at high frequencies
    unwrapped_phase_matrix = np.unwrap(phase_matrix, axis=0)  # Unwrap along the frequency axis (axis=0)
    
    return unwrapped_phase_matrix

def compute_eta(level_vector, unwrapped_phase):
    eta = 0.0732*2*level_vector/unwrapped_phase
    return eta

def fit_polynomial(frequencies, values, degree, freq_range=None):
    """
    Fits a polynomial to the data using linear regression.

    Parameters:
        frequencies (ndarray): Vector of frequency values.
        values (ndarray): Vector of values (e.g., phase or magnitude) to fit.
        degree (int): Degree of the polynomial.
        freq_range (tuple, optional): Range of frequencies to consider (min_freq, max_freq).
                                      If None, the full range is used.

    Returns:
        ndarray: Coefficients of the fitted polynomial (highest degree first).
    """
    # Step 1: Apply range filter if specified
    if freq_range is not None:
        min_freq, max_freq = freq_range
        mask = (frequencies >= min_freq) & (frequencies <= max_freq)
        freq_filtered = frequencies[mask]
        values_filtered = values[mask]
    else:
        freq_filtered = frequencies
        values_filtered = values

    # Step 2: Perform polynomial regression
    coefficients = np.polyfit(freq_filtered, values_filtered, degree)
    
    return coefficients

def olf_improved(b, h, rho, file_path, pointDistance, windowSize = 21, noise_coef = 0, noise_type = "phase" ):
    """
    Olf procedure without the assumption of small damping. b,h,rho only necessary
    to compute the youngs modulus.

    Parameters:
        b: width of a beam under OLF conditions
        h: thickness of beam under OLF conditions
        rho: density of the whole beam
        file_path: path to data from FEM simulation that contains velocity measurements in pairs
        it is assumed that the first half of all columns are paired with the second half of all columns
        values (ndarray): Vector of values (e.g., phase or magnitude) to fit.
        pointDIstance: separation between two measurements points of a pair.
        windowSize: moving average that is implemented in OLF for smoothing
                    reasons.The window size is in their process 21.
        

    Returns:
        eta: value of vector of values of loss factor for each frequency
        cb: value or vector of values of phase speed for each frequency
        E: value or vector of values of Youngs modulus for each frequency. 
        frequencies
    """
    frequencies, responses = SS.openFEMresults(file_path)
    if abs(noise_coef - 0)>1e-10:
        responses = BW.addNoise(responses, noise_coef, noise_type)
    n = int(responses.shape[1]/2)
    complex_matrix = responses[:,-1*n:]/responses[:,:n]
    smoothed_matrix = moving_average_complex_matrix(complex_matrix, windowSize)
    average_vector, abs_average_vector, level_vector = process_smoothed_data(smoothed_matrix)
    
    phase_lag = compute_phase_lag(average_vector)
    phaseLag_pi = phase_lag/np.pi
    kReal =-1*phase_lag/(1*pointDistance)
    kIm = np.log(abs_average_vector)/(-1*pointDistance)
    k_ = kReal -(1j*kIm)
    
    omega = 2* np.pi * frequencies
    k_4 = np.power(k_,4)
    C4_ = np.power(omega,4)/(k_4)
    C4 = np.real(C4_)
    eta = np.imag(C4_)/C4
    cb = np.power(C4,0.25)
    
    mu = b*h*rho
    I = b*h*h*h/12
    E = C4 * mu/I/omega/omega
    
    return eta, cb, E, frequencies

def save_scatter_plot(file_path,remark, x_vector, y_vector, title):
    """
    Create a scatter plot from x_vector and y_vector and save it as a PNG file.

    Parameters:
        file_path (str): The file path where the PNG should be saved.
        x_vector (list or array-like): The x values for the scatter plot.
        y_vector (list or array-like): The y values for the scatter plot.
    """
    # Check if x_vector and y_vector have the same length
    if len(x_vector) != len(y_vector):
        raise ValueError("x_vector and y_vector must have the same length.")
    
    # Create the scatter plot
    plt.figure(figsize=(8, 6))
    plt.plot(x_vector, y_vector, alpha=0.7, color='black') #changed it from plt.scatter to not have circles. color used to be edgecolors
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Loss Factot [-]")
    plt.title(title)
    plt.grid(True)
    #plt.ylim(0.15, 0.25)
    # plt.xlim(0,2000)
    # plt.gca().yaxis.set_major_locator(MultipleLocator(0.1))  # Major ticks every 0.1
    # plt.gca().yaxis.set_minor_locator(MultipleLocator(0.05))  # Minor ticks every 0.05
    # Generate the output file path with .png extension
    output_file = os.path.splitext(file_path)[0] +remark+ ".png"
    output_file_csv = os.path.splitext(file_path)[0] +remark+ ".csv"
    output_file_dat = os.path.splitext(file_path)[0] +remark+ ".dat"
    # Save the plot
    plt.savefig(output_file, dpi=300)
    plt.close()
    #plt.show()
    print(f"Scatter plot saved to: {output_file}")
    # Save data to CSV
    np.savetxt(output_file_dat, np.column_stack([x_vector, y_vector]), header='freq eta')


if __name__ == "__main__":
    # Example usage
    file_path = "exl/20grad.xls"

    
    if os.path.exists(file_path):
        print("File found!")
    else:
        print("File not found. Please check the file path.")
        
    excel_data = pd.ExcelFile(file_path)

    fake_angle = []
    af = pd.read_excel(file_path, sheet_name="1", skiprows=5, header=None)
    fake_angle = af.iloc[:,3].to_numpy()
    fake_angle = fake_angle[:1591]

    frequencies, complex_matrix = read_measurement_data(file_path)
    frequencies = frequencies[-1591:]

    # Display the results
    print("Frequencies:", frequencies)
    print("Complex Data Matrix Shape:", complex_matrix.shape)

    # Example usage
    window_size = 21  # Adjust the window size as needed (e.g., 17th to 37th frequency)
    smoothed_matrix = moving_average_complex_matrix(complex_matrix, window_size)
    smoothed_matrix = smoothed_matrix[-1591:,:]
    # Example usage
    average_vector, abs_average_vector, level_vector = process_smoothed_data(smoothed_matrix)

    # Example usage
    phase_lag = compute_phase_lag(average_vector)
    nopi = phase_lag/np.pi

    eta = compute_eta(level_vector, nopi)
    eta2 = compute_eta(level_vector, fake_angle)

    coefficients = fit_polynomial(frequencies, eta, 3, freq_range=(1074,6000))
    coefficients2 = fit_polynomial(frequencies, eta2, 3, freq_range=(1074,6000))


    evaluated_values = np.polyval(coefficients, frequencies)
    evaluated_values2 = np.polyval(coefficients2, frequencies)

    plt.figure(figsize=(12, 6))

    # Scatter plot for the real part
    plt.subplot(1, 2, 2)
    plt.scatter(frequencies, evaluated_values,s = 1, c="blue", label="OLF_loss Factor")
    plt.scatter(frequencies, evaluated_values2, c="red",s = 1, label="OLF_Loss Factor, Default Phase")
    plt.axhline(0, color="black", linestyle="--", linewidth=0.8)
    plt.title("Loss factor vs Frequency")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Loss Factor [-]")
    plt.legend()


    # Scatter plot for the real part
    plt.subplot(1, 2, 1)
    plt.scatter(frequencies, nopi,s = 1, c="blue", label="Total Phase")
    plt.scatter(frequencies, fake_angle, c="red",s = 1, label="Default Total Phase")
    plt.axhline(0, color="black", linestyle="--", linewidth=0.8)
    plt.title("Total Phase vs Frequency")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Total Phase [-]/pi]")
    plt.legend()