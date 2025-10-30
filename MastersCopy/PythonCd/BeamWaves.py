# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 11:03:26 2024
complitaion of usefull functions regarding beam waves
some analytical solutions for semi infinite beams. They still need to be 
verified to make sure the coefficients are analytically correct. 
@author: prqrz
"""
import numpy as np
import matplotlib.pyplot as plt
class EBBeam:
    def __init__(self, h = 0.001, b = 0.05, E = 210000000000, rho = 7850):
        self.h = h
        self.b = b
        self.E = E
        self.rho = rho
        self.I = b*h*h*h/12
        self.A = b*h
        self.mu = self.A * self.rho
        
    def BendingWaveSpeed(self, frequency, plot=""):
        omega = frequency*2*np.pi
        cb = np.power(self.E*self.I*omega*omega/self.mu, 0.25)
        if plot != "":
            plt.figure(figsize=(12, 6))
            #plt.subplot(1, 2, 2)
            plt.scatter(frequency, cb ,s = 1, c="blue", label="Phase Speed E.B.beam")
            plt.axhline(0, color="black", linestyle="--", linewidth=0.8)
            plt.title(plot)
            plt.xlabel("Frequency [Hz]")
            plt.ylabel("Phase Speed (m/s)")
            plt.legend()
            plt.show()
        return cb
    
    def WaveNumber(self, frequency, plot= ""):
        omega = frequency*2*np.pi
        k = omega/ self.BendingWaveSpeed(frequency)
        if plot != "":
           plt.figure(figsize=(12, 6))
           #plt.subplot(1, 2, 2)
           plt.scatter(frequency, k ,s = 1, c="blue", label="WaveNumber E.B.beam")
           plt.axhline(0, color="black", linestyle="--", linewidth=0.8)
           plt.title(plot)
           plt.xlabel("Frequency [Hz]")
           plt.ylabel("K [1/m]")
           plt.legend()
           plt.show()
        return k
    def WaveLength(self, frequency, plot = ""):
        wavelength = 2*np.pi / self.WaveNumber(frequency)
        if plot != "":
           plt.figure(figsize=(12, 6))
           #plt.subplot(1, 2, 2)
           plt.scatter(frequency, wavelength ,s = 1, c="blue", label="WaveNumber E.B.beam")
           plt.axhline(0.2, color="black", linestyle="--", linewidth=0.8)
           plt.title(plot)
           plt.xlabel("Frequency [Hz]")
           plt.ylabel("Lambda [m]")
           plt.legend()
           plt.show()
        return wavelength
def semiInfBeam(EBBeam, case, excit, Xpos, frequency,eta = 0, plot = ""):
    
    k_ = EBBeam.WaveNumber(frequency)
    kim = eta*k_/4 #assuming small eta < 0.2
    k_ = k_ - 1j*kim
    #define matrix A from k complex
    kivect = np.zeros((1,4),dtype = complex)
    kivect[0,0] = -1j*k_
    kivect[0,1] = -k_
    kivect[0,2] = 1j*k_
    kivect[0,3] = k_
    kivect = np.outer(Xpos,kivect)

    A = np.ones((Xpos.size,4),dtype=complex)
    A = A * np.e
    A = A**kivect
    
    Coeffs = np.zeros(4, dtype = complex)
    if case == "displacement":
        Coeffs[0] = excit/2 #travelling wave out
        Coeffs[1] = excit/2 #nearfield left
    if case == "force":
        Coeffs[0] = (excit/(2*beam.E*beam.I))*(1+1j)
        Coeffs[1] = (excit/(2*beam.E*beam.I))*(1+1j)
    if case == "clamped-displacement":
        Coeffs[0] = (excit/2)*(1+1j)
        Coeffs[1] = (excit/2)*(1-1j)
    displacement = A@Coeffs
    term0 = A[:,0]*Coeffs[0]
    term1 = A[:,1]*Coeffs[1]
    term2 = A[:,2]*Coeffs[2] #not relevant because semi-infinite beam has only 2 terms
    term3 = A[:,3]*Coeffs[3] #not relevant because semi-infinite beam has only 2 terms
    
    if plot != "":
        plt.figure(figsize=(12, 6))
        plt.scatter(Xpos, np.real(displacement) ,s = 1, c="blue", label="Displacement")
        plt.scatter(Xpos, np.real(term0) ,s = 0.5, c="red", label="Forward travelling")
        plt.scatter(Xpos, np.real(term1) ,s = 0.5, c="black", label="Nearfield")
        plt.scatter(Xpos, np.real(term1)/np.abs(term0) ,s = 1, c="green", label="Nearfiel/Farfield")
        plt.axhline(0, color="black", linestyle="--", linewidth=0.8)
        plt.title(plot+", eta = "+str(eta)+", freq = "+str(frequency)+" Hz, excit = "+str(case))
        plt.xlabel("Position [m]")
        plt.ylabel("|w| [m]")
        plt.legend()
        plt.show()
    return displacement 
def fabricateVibration(eta=0.001, cb=65, freq=500, Vp=1, Vpj=0, Vm=0, Vmj=0, plot=""):
    # Vp = 1
    # Vpj = 0
    # Vm = 0
    # Vmj = 0
    
    Coefs = [Vp, Vpj, Vm, Vmj]
    Coefs = np.array(Coefs)
    kivect = np.zeros((1,4), dtype = complex)
    
    # eta = 0
    # cb = 765
    # freq = 300
    
    k_ = 2*np.pi*freq/cb * (1-1j*eta/4)
    
    kivect[0,0] = -1j*k_
    kivect[0,1] = -k_
    kivect[0,2] = 1j*k_
    kivect[0,3] = k_
    
    #Define X position vector
    start = 0.070
    step = 0.020
    x_size = 21
    stop = start + (x_size * step)
    Xpos = np.arange(start,stop,step)
    
    exp = np.outer(Xpos,kivect)
    
    A2 = np.ones((Xpos.size,4),dtype=complex)
    A2= A2 * np.e
    A2 = A2**exp
    
    faveVs = A2@Coefs
    
    if plot != "":
        plt.figure()
        plt.plot(Xpos, np.real(faveVs), marker='o')  # Add marker for clarity if desired
        plt.xlabel('X-axis Label')
        plt.ylabel('Y-axis Label')
        plt.title('X vs Y Plot')
        plt.title(plot+ " real")
        plt.show()
        
        plt.figure()
        plt.plot(Xpos, np.imag(faveVs), marker='o')  # Add marker for clarity if desired
        plt.xlabel('X-axis Label')
        plt.ylabel('Y-axis Label')
        plt.title('X vs Y Plot')
        plt.title(plot+" imag")
        plt.show()
    return faveVs
if __name__ == "__main__":  
    beam = EBBeam(b = 0.05, E = 100e9, h = 0.001 )
    frequencies = np.linspace(0,2000,201)
    cb = beam.BendingWaveSpeed(frequencies, plot = "h=1mm, b=30mm, Steel")  
    #k = beam.WaveNumber(frequencies, plot = "h=1mm, b=50mm, Steel")
    lam = beam.WaveLength(frequencies, plot = "h=1mm, b=50mm, Steel")
    Xpos = np.linspace(0,10,2001)
    beam2 = EBBeam(b= 0.03, E=250000000000)
    #cb2 = beam2.BendingWaveSpeed(frequencies, plot = "h=1mm, b= 30mm, ToughSteel")
    #lam2 = beam2.WaveLength(frequencies, plot = "h=1mm, b= 30mm, ToughSteel")
    #disp = semiInfBeam(beam, "displacement", 1, Xpos, 30,0.01, "Displacement vs Position x, disp excitation")
    #disp2 = semiInfBeam(beam, "force", 1, Xpos, 30,0.01, "Displacement vs Position x, force excitation")
    #disp3 = semiInfBeam(beam, "clamped-displacement", 1, Xpos, 30,0.0, "Displacement vs Position x")
    #disp4 = semiInfBeam(beam, "clamped-displacement", 1, Xpos, 30,0.5, "Displacement vs Position x")
    #fabricateVibration(eta=0.001, cb=60, freq=300, Vp=1, Vm=0, Vpj=0, Vmj=0, plot="HI")
        