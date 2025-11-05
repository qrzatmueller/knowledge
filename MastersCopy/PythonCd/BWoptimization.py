# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 11:30:59 2025
In this file an optimization procedure is used to identify
the complex k from the measurements/ simulation values either from OLF or BW. 
@author: prqrz
"""
import numpy as np
from scipy.optimize import minimize
import SolutionSurface as SS
import os
import BeamWaves as BW
import matplotlib.pyplot as plt
import olf
# Define the matrix A(k) based on k
def compute_A(k, Xpos):
    k_complex = k[0] + 1j * k[1]  # k[0] = real(k), k[1] = imag(k)
    A = np.column_stack([
        np.exp(-1j*k_complex * Xpos),
        np.exp(-k_complex * Xpos),
        np.exp(1j*k_complex * Xpos),
        np.exp(k_complex * Xpos),
    ])
    return A

# Define the objective function for k
def Objective(k, b, Xpos):
    A = compute_A(k, Xpos)
    # Solve least squares for x
    x, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    # Compute the residual norm
    residual = b-A@x
    # norm_real = np.linalg.norm(residual.real)
    # norm_imag = np.linalg.norm(residual.imag)
    # weightedNorm = norm_real + 10*norm_imag
    # norm = weightedNorm
    norm = np.linalg.norm(residual**2)
    return norm
def ObjectiveNormal(k,b,Xpos):
    A = compute_A(k, Xpos)
    
    column_norms = np.linalg.norm(A,axis=0)#this normalization is new because from olf measurements it wasnt working without it, maybe because of big x?
    Anorm  = A / np.maximum(column_norms, 1e-10)
    # Solve least squares for x
    x, residualNorm, _, _ = np.linalg.lstsq(Anorm, b, rcond=None)
    # Compute the residual norm
    residual = b-A@x
    norm = np.linalg.norm(residual**2)
    return residualNorm
def Optimize(freq, Xpos, meas, objective):
    k0 = [0.0,0.0]#use 0.2 for eta and the
    beam = BW.EBBeam(b = 0.03) # for initial guess use steel beam
    starterEta = 0.5 #given the algorithm of Nelder-Mead (uses 5% change for vertices of Simplex), using a higher eta is better than close to 0, to give it a big area to look
    k_ = SS.GetKcomplex(starterEta, beam.E, beam.I, freq, beam.mu)
    k0 = [k_.real,k_.imag]
    result = minimize(objective, k0, args =(meas, Xpos),method = 'Nelder-Mead', tol = 1e-18)
    return result
def IdentifyEta_YACS(result_path, case, plot = "", remark = ""):
    """
    Using the BW procedure to identify the loss factor, using an optimization technique and
    no approximation of eta w.r.t. k' and k''.

    Parameters:
        result_path: file from simulation to get data from
        case: either BW for bending wave method or OLF for olf
        
    Returns:
        eta: value of vector of values of loss factor for each frequency
        cb: value or vector of values of phase speed for each frequency
        E: value or vector of values of Youngs modulus for each frequency. 
    """
    file_path = result_path
    frequencies, responses = SS.openFEMresults(file_path)
    if (case == "BW"):
        #Xpos for usual BW, I have not retested this
        start = 0.070 #it used to be 0.17, but that is completely wrong, although no real effect.  
        step = 0.020
        x_size = 21
        stop = start + ((x_size-1) * step)
        Xpos = np.linspace(start, stop, x_size)
        beam = BW.EBBeam(b=0.05, E=100e9) # for initial guess use steel beam
    if (case == "BWsplit"):
        #Xpos for usual BW, I have not retested this
        start = 0.070 #it used to be 0.17, but that is completely wrong, although no real effect.  
        step = 0.020
        x_size = 21
        stop = start + ((x_size-1) * step)
        Xpos = np.linspace(start, stop, x_size)
        beam = BW.EBBeam(b=0.05, E=100e9) # for initial guess use steel beam
    if (case == "OLF"):    
        #Xpos for OLF, and guessing beam: 
        start = 0.170
        step = 0.01 #it cannot work like this because the steps are not regelmäßig
        x_size = 12 #meaning one step but 2 datapoints
        stop = start + ((x_size-1) * step)
        Xpos = np.linspace(start, stop, x_size)
        Xpos2 = np.linspace(start+0.2, stop+0.2, x_size)
        Xpos = np.concatenate((Xpos, Xpos2), axis = 0)
        beam = BW.EBBeam(b = 0.03) # for initial guess use steel beam
        
    etas = np.zeros(frequencies.shape)
    cbs = np.zeros(frequencies.shape)
    k_s = np.zeros(frequencies.shape,dtype=complex)
    for i, freq in enumerate(frequencies):
        meas = responses[i,:] 
        print("Xposshape", Xpos.shape)
        print("meas shape", meas.shape)
        result = Optimize(freq, Xpos, meas, ObjectiveNormal)
        k_ = result.x[0]+1j*result.x[1]
        if (result.x[1]>0):
            print("imaginary part is positive at frequency "+str(freq)+ ", check formulation")
        omega = 2*np.pi*freq
        E_ = beam.mu * omega*omega / beam.I /(pow(k_,4))
        etas[i] = E_.imag / E_.real
        etas[i] = -(k_**4).imag / (k_**4).real #there is something funny here that I had to add the minus sign, but the line above delivers same result but with positive sign...
        k_s[i] = k_
        cbs[i] = (pow(E_*beam.I*omega*omega/beam.mu,0.25)).real
        cbs[i] = (omega/k_).real 
    if (plot != ""):
        olf.save_scatter_plot(file_path,"BWIDn"+remark, frequencies, etas, plot)
    return etas, cbs, E_, k_s 
def IdentifyEta_2Parts(result_path, case, plot = "", remark = ""):
    """
    Using the BW procedure to identify the loss factor, using an optimization technique and
    no approximation of eta w.r.t. k' and k''. Used for the beam with the small inhomogeneity
    where values are identified separately for the beam before and for the beam after the 
    inhomo

    Parameters:
        result_path: file from simulation to get data from
        case: either BW for bending wave method or OLF for olf
        
    Returns:
        eta: value of vector of values of loss factor for each frequency
        cb: value or vector of values of phase speed for each frequency
        E: value or vector of values of Youngs modulus for each frequency. 
    """
    file_path = result_path
    frequencies, responses = SS.openFEMresults(file_path)
    if (case == "BW"):
        #Xpos for usual BW, I have not retested this
        start = 0.070 #it used to be 0.17, but that is completely wrong, although no real effect.  
        step = 0.020
        x_size = 21
        stop = start + ((x_size-1) * step)
        Xpos = np.linspace(start, stop, x_size)
        beam = BW.EBBeam(b=0.05, E=100e9) # for initial guess use steel beam
    if (case == "BWsplit"):
        #Xpos for usual BW, I have not retested this
        start = 0.070 #it used to be 0.17, but that is completely wrong, although no real effect.  
        step = 0.020
        x_size = 7
        stop = start + ((x_size-1) * step)
        Xpos = np.linspace(start, stop, x_size)
        beam = BW.EBBeam(b=0.05, E=100e9) # for initial guess use steel beam

        nInhomo = 1 #number of nodes in the inhomogeneity or in general nodes in the middle that should not be part of the first or second beam
        start2 = stop + (nInhomo+1)*step
        x_size2 = 13
        stop2 = start2 +((x_size2-1)*step)
        Xpos2 = np.linspace(start2, stop2, x_size2)
        
    etas = np.zeros(frequencies.shape)
    print("etasshape", etas.shape)
    print("freqshape",frequencies.shape)
    cbs = np.zeros(frequencies.shape)
    k_s = np.zeros(frequencies.shape,dtype=complex)
    A_s = np.zeros((frequencies.size,4), dtype=complex)

    etas2 = np.zeros_like(etas)
    cbs2 = np.zeros_like(cbs)
    k_s2 = np.zeros_like(k_s)
    A_s2 = np.zeros_like(A_s)

    
    for i, freq in enumerate(frequencies):
        meas = responses[i,0:x_size] #testing change here it was responses[i,:]
        # print("Xposshape", Xpos.shape)
        # print("meas shape", meas.shape)
        result = Optimize(freq, Xpos, meas, ObjectiveNormal)
        k_ = result.x[0]+1j*result.x[1]
        if (result.x[1]>0):
            print("imaginary part is positive at frequency "+str(freq)+ ", check formulation")
        ehochk = compute_A(result.x, Xpos)
        # Solve least squares for x
        x, _, _, _ = np.linalg.lstsq(ehochk, meas, rcond=None)
        # print("xshape", x.shape)
        A_s[i] = x
        omega = 2*np.pi*freq
        E_ = beam.mu * omega*omega / beam.I /(pow(k_,4))
        etas[i] = E_.imag / E_.real
        etas[i] = -(k_**4).imag / (k_**4).real #there is something funny here that I had to add the minus sign, but the line above delivers same result but with positive sign...
        k_s[i] = k_
        cbs[i] = (pow(E_*beam.I*omega*omega/beam.mu,0.25)).real
        cbs[i] = (omega/k_).real 
    if (plot != ""):
        olf.save_scatter_plot(file_path,"BWIDn1"+remark, frequencies, etas, plot)

    for i, freq in enumerate(frequencies):
        meas2 = responses[i,x_size+1::] #testing change here it was responses[i,:]
        # print("Xposshape", Xpos2.shape)
        # print("meas shape", meas2.shape)
        result = Optimize(freq, Xpos2, meas2, ObjectiveNormal)
        k_ = result.x[0]+1j*result.x[1]
        if (result.x[1]>0):
            print("imaginary part is positive at frequency "+str(freq)+ ", check formulation")
        ehochk = compute_A(result.x, Xpos2)
        # Solve least squares for x
        x, _, _, _ = np.linalg.lstsq(ehochk, meas2, rcond=None)
        # print("xshape", x.shape)
        A_s2[i] = x
        omega = 2*np.pi*freq
        E_2= beam.mu * omega*omega / beam.I /(pow(k_,4))
        etas2[i] = E_.imag / E_.real
        etas2[i] = -(k_**4).imag / (k_**4).real #there is something funny here that I had to add the minus sign, but the line above delivers same result but with positive sign...
        k_s2[i] = k_
        cbs2[i] = (pow(E_*beam.I*omega*omega/beam.mu,0.25)).real
        cbs2[i] = (omega/k_).real 
    if (plot != ""):
        olf.save_scatter_plot(file_path,"BWIDn2"+remark, frequencies, etas2, plot)
    return etas, cbs, E_, k_s,A_s, etas2, cbs2, E_2, k_s2, A_s2
def addNoise(signal, noise_coef, noise_type):
    if noise_type == "phase":
        signal_phase = np.angle(signal)
        signal_magnitude = np.abs(signal)
        noise_std = abs(signal_phase*noise_coef)
        noise = np.random.normal(0,noise_std)
        noisy_phase = signal_phase+noise
        noisy_signal = signal_magnitude*np.exp(1j*noisy_phase)
    else:
        signal_magnitude = np.abs(signal)
        noise_std= signal_magnitude*noise_coef
    
        real_noise = np.random.normal(0,noise_std)
        imag_noise = np.random.normal(0,noise_std)
    
        noisy_signal = signal + real_noise + 1j*imag_noise
    return noisy_signal
def IdentifyEta_YACS_Noise(result_path,noise_coef, case, plot = "", remark = "", noise_type =""):
    """
    Using the BW procedure to identify the loss factor, using an optimization technique and
    no approximation of eta w.r.t. k' and k''.
    This version of the funciton incorporates noise into the measurement. It can be a general
    noise for the imaginary and real part, or a noise only on the phase. 

    Parameters:
        result_path: file from simulation to get data from
        case: either BW for bending wave method or OLF for olf
        
    Returns:
        eta: value of vector of values of loss factor for each frequency
        cb: value or vector of values of phase speed for each frequency
        E: value or vector of values of Youngs modulus for each frequency. 
    """
    file_path = result_path
    frequencies, responses = SS.openFEMresults(file_path)
    
    #noise_type = "phase"
    if noise_type == "phase": 
        signal_phase = np.angle(responses)
        signal_magnitude = np.abs(responses)
        noise_std = abs(signal_phase*noise_coef)
        noise = np.random.normal(0,noise_std)
        noisy_phase = signal_phase+noise
        noisyresponses = signal_magnitude*np.exp(1j*noisy_phase)
    elif noise_type == "max_noise":
        signal_magnitude = np.abs(responses)
        noise_std= signal_magnitude*noise_coef
        
    
        noisyresponses = responses + noise_std + 1j*noise_std
    elif noise_type == "max_noise_sign":
        signal_magnitude = np.abs(responses)
        noise_std= signal_magnitude*noise_coef

        real_response = responses.real + noise_std * np.sign(responses.real)
        imag_response = responses.imag + noise_std * np.sign(responses.imag)
    
        #noisyresponses = responses + noise_std + 1j*noise_std
        noisyresponses = real_response + 1j*imag_response
    else:
        signal_magnitude = np.abs(responses)
        noise_std= signal_magnitude*noise_coef
    
        real_noise = np.random.normal(0,noise_std)
        imag_noise = np.random.normal(0,noise_std)
    
        noisyresponses = responses + real_noise + 1j*imag_noise
    if (case == "BW"):
        #Xpos for usual BW, I have not retested this
        start = 0.070
        step = 0.020
        x_size = 21
        stop = start + ((x_size-1) * step)
        Xpos = np.linspace(start, stop, x_size)
        beam = BW.EBBeam(b=0.05) # for initial guess use steel beam
    if (case == "OLF"):    
        #Xpos for OLF, and guessing beam: 
        start = 0.170
        step = 0.01 #it cannot work like this because the steps are not regelmäßig
        x_size = 12 #meaning one step but 2 datapoints
        stop = start + ((x_size-1) * step)
        Xpos = np.linspace(start, stop, x_size)
        Xpos2 = np.linspace(start+0.2, stop+0.2, x_size)
        Xpos = np.concatenate((Xpos, Xpos2), axis = 0)
        beam = BW.EBBeam(b = 0.03) # for initial guess use steel beam
        
    etas = np.zeros(frequencies.shape)
    cbs = np.zeros(frequencies.shape)
    for i, freq in enumerate(frequencies):
        meas = noisyresponses[i,:]
        result = Optimize(freq, Xpos, meas, ObjectiveNormal)
        k_ = result.x[0]+1j*result.x[1]
        #print(k_) #todo remove
        if (result.x[1]>0):
            print("imaginary part is positive at frequency "+str(freq)+ ", check formulation")
        omega = 2*np.pi*freq
        E_ = beam.mu * omega*omega / beam.I /(pow(k_,4))
        etas[i] = E_.imag / E_.real
        etas[i] = -(k_**4).imag / (k_**4).real #there is something funny here that I had to add the minus sign, but the line above delivers same result but with positive sign...
        
        cbs[i] = pow(E_.real*beam.I*omega*omega/beam.mu,0.25)
        cbs[i] = (omega/k_).real   
    if (plot != ""):
        olf.save_scatter_plot(file_path,"BWID"+remark, frequencies, etas, plot)
    return etas, cbs, E_ , frequencies #todo I added thhis frequencies but is not maybe affects other functions that call this one
if __name__ == "__main__":
    current_file_path = os.path.abspath(__file__)
    mt_eq_path = os.path.abspath(os.path.join(current_file_path, '..', '..', '..'))

    #file_path = os.path.join(mt_eq_path, 'CodeAsterModels/YACSauto/YACS_BW/autoRes/res0.txt')#get freqiemcies and measurements
    file_path = os.path.join(mt_eq_path, 'CodeAsterModels/OLF/autoRes/03OlfId/res2.txt')
    #do a loop in which I iterate through the frequencies
    frequencies, responses = SS.openFEMresults(file_path)
    
    '''
    #example of using Identify Eta YACS, so using the optimization prcedure to 
    identify olf loss factor of simulated data through yacs
    '''

    #etas, cbs,_ = IdentifyEta_YACS(file_path, "OLF", plot="Loss Factor vs Frequency (BW)")
    
    '''
    #example of noise added to simulation
    '''

    # etas,cbs,_ = IdentifyEta_YACS_Noise(file_path,0.01, case = "BW", plot = "identified Eta from 1% Noisy Simulation", remark= "Noise01pct")
    #etas,cbs,_ = IdentifyEta_YACS_Noise(file_path,0.00, case = "BW", plot = "identified Eta from 1% Noisy Simulation", remark= "Noise01pct")
    # etas,cbs,_ = IdentifyEta_YACS_Noise(file_path,0.01, case = "BW", plot = "identified Eta from 1% PH-Noisy Simulation", remark= "Ph_Noise01pct")
    # etas,cbs,_ = IdentifyEta_YACS_Noise(file_path,0.05, case = "BW", plot = "identified Eta from 5% Noisy Simulation", remark= "Noise05pct")
    # etas,cbs,_ = IdentifyEta_YACS_Noise(file_path,0.05, case = "BW", plot = "identified Eta from 5% PH-Noisy Simulation", remark= "Ph_Noise05pct")
    '''
    
    #example for a single frequency:
    '''
    num= 145 #index of the single frequency to check, for all frequencies see example below below. 
    freq = frequencies[num]
    meas = responses[num,:]
    print("freq = ", freq)

    
    #Define X position vector, this should coincide with the amount of measurements in the .txt
    
    #Xpos for OLF, and guessing beam: 
    start = 0.170
    step = 0.01 #it cannot work like this because the steps are not regelmäßig
    x_size = 12 #meaning one step but 2 datapoints
    stop = start + ((x_size-1) * step)
    Xpos = np.linspace(start, stop, x_size)
    Xpos2 = np.linspace(start+0.2, stop+0.2, x_size)
    Xpos = np.concatenate((Xpos, Xpos2), axis = 0)
    #print(Xpos)
    beam = BW.EBBeam(b = 0.03) # for initial guess use steel beam
    
    #Xpos for usual BW, I have not retested this
    # start = 0.170
    # step = 0.020
    # x_size = 21
    # stop = start + ((x_size-1) * step)
    # Xpos = np.linspace(start, stop, x_size)
    # beam = BW.EBBeam(b=0.05) # for initial guess use steel beam
    
    # k0 = [0.0,0.0]#use 0.2 for eta and the
    
    starterEta = 0.01
    k_ = SS.GetKcomplex(starterEta, beam.E, beam.I, freq, beam.mu)
    print("starter eta "+ str(starterEta))
    print("steel beam k real = "+str(k_.real))
    k0 = [k_.real*0.9,k_.imag]
    # def constraint_ineq(k):
    #     return k[1] + 0.1  # Imaginary part should be greater than or equal to 0.1
    # con = {'type':'ineq', 'fun': constraint_ineq}
    result = minimize(ObjectiveNormal, k0, args =(meas, Xpos),method = 'Nelder-Mead', tol = 1e-18)
    print("Optimal k (real, imag):", result.x)
    print("Minimum residual:", result.fun)
    if result.fun >0.1:
        print("Alaaaarm, residual too high, check starting point.")
    k_ = result.x[0]+1j*result.x[1] #check here if the imaginary part is negative or not
    if (result.x[1]>0):
        print("imaginary part is positive!!, check formulation")
    omega = 2*np.pi*freq
    E_ = beam.mu * omega*omega / beam.I /(pow(k_,4))
    eta = E_.imag / E_.real
    eta2 = -k_.imag * 4 / k_.real 
    print("eta = "+str(eta)+" aprox eta = "+ str(eta2))
    #     #within loop calculate the ks
    #     #minimize
    #     #profit
    '''
    Attempt to calculate the coeficients of simulated OLF measurements to check
    if indeed there are only a few waves included. at least not the returning 
    wave. I will use the single frequency used in the example1 above 
    '''
    A = compute_A(result.x, Xpos)
    
    column_norms = np.linalg.norm(A,axis=0)#this normalization is new because from olf measurements it wasnt working without it, maybe because of big x?
    Anorm  = A / np.maximum(column_norms, 1e-10)
    coefsNorm, residuals, _, _ = np.linalg.lstsq(Anorm, meas, rcond=None)
    # Compute the coeficients
    coefs = coefsNorm / np.maximum(column_norms, 1e-10)
    # for i, coef in enumerate(coefs):
    #     if (abs(coef.real) < 0.0000001):
    #         coefs[i] = 0
    x = np.linspace(0, 0.7,201)
    #x=Xpos
    wave1 = coefs[0]*np.exp(-1j*k_*x)
    wave2 = coefs[1]*np.exp(-1*k_*x)
    wave3 = coefs[2]*np.exp(1j*k_*x)
    wave4 = coefs[3]*np.exp(1*k_*x)
    waveTot = wave1 + wave2 + wave3 + wave4
    output_file_dat = os.path.join(mt_eq_path, 'CodeAsterModels/OLF/autoRes/03OlfId/waveTerm2.dat')
    np.savetxt(output_file_dat, np.column_stack([x, wave1.real, wave2.real, wave3.real, wave4.real, waveTot.real]), header='x w1 w2 w3 w4 wt')

    aaa= plt.figure(figsize=(8, 6))
    # plt.subplot(1, 2, 1)
    plt.scatter(x, wave1.real,s = 5, c="blue", label="Forward ")
    plt.scatter(x, wave2.real,s = 5, c="red", label="Diminishing")
    plt.scatter(x, wave3.real,s = 5, c="black", label="Backward ")
    plt.scatter(x, wave4.real,s = 5, c="pink", label="Increasing ")
    plt.scatter(x, waveTot.real,s = 8, c="green", label="TOTAL ")
    for xi in Xpos:
        plt.axvline(xi, color="red", linestyle="--", alpha=0.7, label=f"x = {xi:.2f}")
    plt.axhline(0, color="black", linestyle="--", linewidth=0.8)
    plt.title("Terms of OLF Simulation at 30 Hz")
    plt.xlabel("X Position")
    plt.ylabel("Amplitude (Real)")
    plt.legend()
    plt.show()

    '''
    example for a whole measurement/simulation file: 
        using the file CodeAsterModels/OLF/EBlong_Files/resF/res0p2.txt
        here I learned that without normalization in the A matrix, it is not 
        possible to get convergence for higher frequencies. it becomes unstable
        numerically due to big numbers. 
        The coefficients although, cannot be taken from thes linalg.lstqs, since 
        they also have some coefficient inside. I do not know why trying to deno
        denormalized the coeficients did not work. but once k is found. the coefficients
        should be findable. 
        note that the X position here is that for OLF and not BW
    '''
    # etas = np.zeros(frequencies.shape)
    # cbs = np.zeros(frequencies.shape)
    # for i, freq in enumerate(frequencies):
    #     meas = responses[i,:]
    #     result = Optimize(freq, Xpos, meas, ObjectiveNormal)
    #     k_ = result.x[0]+1j*result.x[1]
    #     if (result.x[1]>0):
    #         print("imaginary part is positive at frequency "+str(freq)+ ", check formulation")
    #     E_ = beam.mu * freq*freq / beam.I /(pow(k_,4))
    #     etas[i] = E_.imag / E_.real
    #     omega = 2*np.pi*freq
    #     cbs[i] = pow(E_.real*beam.I*omega*omega/beam.mu,0.25)
    
    
    
        
        
        