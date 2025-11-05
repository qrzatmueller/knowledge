'''
In this file I will attempt to investigate the effects of having a thickness of one and a half of the value
for 10 mm starting at 0.27 of the OLF beam. What happens to the wavefield in terms of real and imaginary part.
We see that the inhomogenity causes a reflection, and for after the inhomogenity there is a delay on the phase
because of the small section(inhomo) in which the phaseVelocity and the wavelength were larger (for a fixed frequency)

THe files res480_th10_4000.txt are using a softer material E = 1,000,000,000 Pa, which in turn makes 
the wavelength inappropriately small for high frequencies at 2000 hz around 20mm. for our elements of 5mm
we should not exceed the 30mm wavelength. prefeably have larger that 50 mm wavelength so there are ten
elements per wavelength. This beam reaches 50 mm at 250 Hz. 

at the end of this file, loss factor identification is done, where it can be seen that the identification still
workds in a general sense, when using so many data points.

a second case is introduced for steel, whose wavelength is appropriate for the element length. 
Here I observe that the loss factor identified relatively good, with a maximum value at 0.216 for 97 point and 0.225 for 21 points
the loss factor identified is always bigger that the nominal, what is the relationship with a bigger thickness?. more thickness
more stiffness, bigger lamda for a given frequency. bigger lambda smaller k. but this is all the real part of k. could 
a higher stiffness magnitude be compensated by adding more on the imaginary side of E?
Thinking that this is the OLF bar with BW procedure, shouldnt the amplitude of the wave play a role?
The greatest difference is for lower frequencies reaching a peak at around 100 hz and then lowering. where the value is 0.204
Can it be that bigger thickness with the same eta is also a higher loss of energy that is not taken into account in the equations
since it is expecting a smaller thickness? so higher energy loss for the "same" thickness yields a higher loss factor. But how
is energy seen in the results. The waves have  slightly smaller amplitude which means a higher loss factor. this sounds more promising

At the end of the code I am identifying Eta with OLF. Obs1: It decreases the eta nstead of increasing it. Obs2, as always shitier
identification at lower frequencies, to a point where its just not valid. Obs3: between 150 and 500Hz there are observable oscillations
on the identified eta, not sure what they mean or why they disappear. They are present with and without inh.
Knowing that olf depends on 2 things, the Levelreduction and the phase difference, it makes sense to have a lower eta, 
since we observed that in the mayority of the cases, after the inhomogenity there is a delay in the phase. this makes the phase 
difference larger and reduces the loss factor. A way to test this is to put the inhomogenity at the beginning so it affects fewer 
measurement points. 

open questions in this example: If we have a material that always has a wavelength appropriate for the 
5mm element size,would the downwards trend with frequency subside?
If the inhomogenity was larger in width, how would this affect. 
What happens if the inhomogenity also contains a lower eta, does this affect?
what is the mechanism with which the inhomogenities affect the identification?
if less points are used for the identification, what happens?


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

current_file_path = os.path.abspath(__file__)
mt_eq_path = os.path.abspath(os.path.join(current_file_path, '..', '..', '..', '..', '..'))
'''
case1, softer material E = 1Gpa, with eta=0.1., inh is 1.5 mm thick 10 mm long yielded a very bad identification when using 21 points (responsesP4)
'''
file_path = os.path.join(mt_eq_path, 'CodeAsterModels/OLF/OLF_Working_Files/resF/res480_th15_4000.txt')
file_path2 = os.path.join(mt_eq_path, 'CodeAsterModels/OLF/OLF_Working_Files/resF/res480_th10_4000.txt')
file_path3 = os.path.join(mt_eq_path,'CodeAsterModels/OLF/OLF_Working_Files/resF/res_E210_eta0p2_Inh.txt')
file_path4 = os.path.join(mt_eq_path,'CodeAsterModels/OLF/OLF_Working_Files/resF/res_E210_eta0p2_NoInho.txt')

label1 = "No inhomogenity d = 0mm, h=1 mm"
label2 = "Inhomogenity d = 10mm, h=1.5mm"
label3 = "Inhomogenity d = 10mm, h=1.5mm"
'''
case2, steel material E = 210Gpa, with eta=0.2. same shape of inhomogenity
'''
# file_path = os.path.join(mt_eq_path, 'CodeAsterModels/OLF/OLF_Working_Files/resF/res480_E210_eta0p2_Inh.txt') #with inhomogenity
# file_path2 = os.path.join(mt_eq_path, 'CodeAsterModels/OLF/OLF_Working_Files/resF/res480_E210_eta0p2_NoInho.txt') #without
# file_path3 = os.path.join(mt_eq_path,'CodeAsterModels/OLF/OLF_Working_Files/resF/res_E210_eta0p2_Inh.txt')
# file_path4 = os.path.join(mt_eq_path,'CodeAsterModels/OLF/OLF_Working_Files/resF/res_E210_eta0p2_NoInho.txt')
# label1 = "No inhomogenity d = 0mm, h=1 mm"
# label2 = "Inhomogenity d = 10mm, h=1.5mm"
# label3 = "Inhomogenity d = 10mm, h=1.5mm"
'''
case3, steel material, nominal eta 0.2 and inh eta 0.3 on top of thickness inhomogenity
'''
# file_path = os.path.join(mt_eq_path, 'CodeAsterModels/OLF/OLF_Working_Files/resF/res480_E210_eta0p2_InhEta0p3.txt') #with inhomogenity
# file_path2 = os.path.join(mt_eq_path, 'CodeAsterModels/OLF/OLF_Working_Files/resF/res480_E210_eta0p2_Inh.txt') #without
# file_path3 = os.path.join(mt_eq_path,'CodeAsterModels/OLF/OLF_Working_Files/resF/res_E210_eta0p2_InhEta0p3.txt')
# file_path4 = os.path.join(mt_eq_path,'CodeAsterModels/OLF/OLF_Working_Files/resF/res_E210_eta0p2_Inh.txt')
# label1 = "Geometric Inhomogenity d = 10mm, h=1.5mm"
# label2 = "Geometric + Loss factor Inh. 97 Ppints d = 10mm, h=1.5mm eta_i=0.3"
# label3 = "Geometric + Loss factor Inh. 21Points d = 10mm, h=1.5mm eta_i=0.3"

'''
begins calculation
'''

case = "fullOLF"
plot = "a"
remark = "test"
frequencies, responses = SS.openFEMresults(file_path)
frequencies2, responses2 = SS.openFEMresults(file_path2)
#frequencies3, responses3 = SS.openFEMresults(file_path3)
frequecies4 = frequencies 
responsesP4 = responses[:,4::4] #take every fourth response starting from the fourth frequency. This should accomodate the BW points
start = 0.070
step = 0.020
x_size = 21
stop = start + ((x_size-1) * step)
Xpos4 = np.linspace(start, stop, x_size)
responsesP4 = responsesP4[:,:21]

#make Xpos for the first half, and the second half separate.
#is the Xpos for BW completely wrong??
if (case == "firstHalf"):
    start = 0.07
    step = 0.02
    x_size = 8
    stop = start + ((x_size-1) * step)
    Xpos = np.linspace(start, stop, x_size)
    beam = beamW.EBBeam(b=0.05) # for initial guess use steel beam
    #cut data to only the first 8 values
    responsesP = responses[:,:8]
if (case == "secondHalf"):
    start = 0.33
    step = 0.02
    x_size = 8
    stop = start + ((x_size-1) * step)
    Xpos = np.linspace(start, stop, x_size)
    beam = beamW.EBBeam(b=0.05) # for initial guess use steel beam

    responsesP = responses[:,-8:]
if (case == "BW"):
    #Xpos for usual BW, I have not retested this
    start = 0.070
    step = 0.020
    x_size = 21
    stop = start + ((x_size-1) * step)
    Xpos = np.linspace(start, stop, x_size)
    responsesP = responses
    
    beam = beamW.EBBeam(b=0.05) # for initial guess use steel beam
if (case == "OLF"):
    #Xpos for OLF, and guessing beam:
    start = 0.170
    step = 0.01 #it cannot work like this because the steps are not regelmäßig
    x_size = 12 #meaning one step but 2 datapoints
    stop = start + ((x_size-1) * step)
    Xpos = np.linspace(start, stop, x_size)
    Xpos2 = np.linspace(start+0.2, stop+0.2, x_size)
    Xpos = np.concatenate((Xpos, Xpos2), axis = 0)
    beam = beamW.EBBeam(b = 0.03) # for initial guess use steel beam
if (case == "fullBW"):
    start = 0.0
    step = 0.005
    stop = 0.48
    x_size = (stop-start)/step
    Xpos = np.linspace(start, stop, int(x_size)+1)
    responsesP = responses
    responsesP2 = responses2
    #responsesP3 = responses3
    beam = beamW.EBBeam(b=0.05) # for initial guess use steel beam
if (case == "fullOLF"):
    start = 0.0
    step = 0.005
    stop = 0.48
    x_size = (stop-start)/step
    Xpos = np.linspace(start, stop, int(x_size)+1)
    responsesP = responses
    responsesP2 = responses2
    #responsesP3 = responses3
    beam = beamW.EBBeam(b=0.05) # for initial guess use steel beam
    
n = 6 #n=8 for 100hz, n=159 for 1620
frequency = frequencies[n]
fig = plt.figure(figsize=(8,6))

ax0 = fig.add_subplot(221)

ax0.plot(Xpos, responsesP[n,:].real, label = "d = 10mm, h=1.5mm")
ax0.plot(Xpos, responsesP2[n,:].real, label = "d = 0, h = 1mm")
plt.axvline(0.27, color="black", linestyle="--", linewidth=0.8)
plt.axvline(0.28, color="black", linestyle="--", linewidth=0.8)
#ax.plot(Xpos, responsesP3[n,:].real, label = "thck1")
# Create reference plane at x = 0.27

ax0.set_xlabel("Position (x)")
ax0.set_ylabel("Real Part")
plt.title(str(frequency)+ "Hz, eta=0.1")
plt.legend()

ax1 = fig.add_subplot(222)
ax1.plot(Xpos, responsesP[n,:].imag, label = "d = 10mm, h=1.5mm")
ax1.plot(Xpos, responsesP2[n,:].imag, label = "d = 0, h = 1mm")
#ax.plot(Xpos, responsesP3[n,:].imag, label = "thck1")
plt.axvline(0.27, color="black", linestyle="--", linewidth=0.8)
plt.axvline(0.28, color="black", linestyle="--", linewidth=0.8)
# Create reference plane at x = 0.27
ax1.set_xlabel("Position (x)")
ax1.set_ylabel("Imag Part")
plt.legend()


ax2 = fig.add_subplot(223)
ax2.plot(Xpos, np.abs(responsesP[n,:]), label = "d = 10mm, h=1.5mm")
ax2.plot(Xpos, np.abs(responsesP2[n,:]), label = "d = 0, h = 1mm")
#ax2.plot(Xpos, np.abs(responsesP3[n,:]), label = "thck1")
plt.axvline(0.27, color="black", linestyle="--", linewidth=0.8)
plt.axvline(0.28, color="black", linestyle="--", linewidth=0.8)
ax2.set_xlabel("Position (x)")
ax2.set_ylabel("Magnitude")
plt.legend()

ax3 = fig.add_subplot(224)
ax3.plot(Xpos, np.angle(responsesP[n,:]), label = "d = 10mm, h=1.5mm")
ax3.plot(Xpos, np.angle(responsesP2[n,:]), label = "d = 0, h = 1mm")
#ax3.plot(Xpos, np.angle(responsesP3[n,:]), label = "thck1")
plt.axvline(0.27, color="black", linestyle="--", linewidth=0.8)
plt.axvline(0.28, color="black", linestyle="--", linewidth=0.8)
ax3.set_xlabel("Position (x)")
ax3.set_ylabel("Phase")
plt.legend()
plt.show()

fig2 = plt.figure(figsize=(8,6))
ax = fig2.add_subplot(111,projection = '3d')

line, = ax.plot(Xpos, responsesP[n,:].real, responsesP[n,:].imag, label = "d = 10mm, h=1.5mm")
#line2, =ax.plot(Xpos, responsesP2[n,:].real, responsesP2[n,:].imag, label = "d = 0, h = 1mm")
#ax.plot(Xpos, responsesP3[n,:].real, responsesP3[n,:].imag, label = "thck1")
proj_line, = ax.plot([], [], [], 'k--', alpha=0.5, label="Projection on Real Plane")  # Projection

# Create reference plane at x = 0.27
x_plane = np.full((10, 10), 0.27)  # A constant x-plane
y_plane = np.linspace(-0.0000001, 0.0000001, 10)  # Y-axis (real part range)
z_plane = np.linspace(-0.0000001, 0.0000001, 10)  # Z-axis (imaginary part range)
Y, Z = np.meshgrid(y_plane, z_plane)  # Create a grid
ax.plot_surface(x_plane, Y, Z, color='gray', alpha=0.3)  # Add semi-transparent plane
ax.set_xlabel("Position (x)")
ax.set_ylabel("Real Part")
ax.set_zlabel("Imag Part")
plt.legend()
'''
3D animation of complex pointer 
'''
# fps = 30
# omega = 0.005*frequency*2*np.pi
# t_max = 30
# num_frames = fps*t_max
# def update(frame):
#     t = frame / fps  # Convert frame to time
#     phase_factor = np.exp(1j * omega * t)  # e^(iωt) term
#     new_displacement = responsesP[n,:] * phase_factor  # Apply time evolution
    
#     # Extract real and imaginary parts
#     y_real = new_displacement.real
#     z_imag = new_displacement.imag
#     # Update the 3D line
#     line.set_data(Xpos, y_real)
#     line.set_3d_properties(z_imag)
    
#     proj_line.set_data(Xpos, y_real)
#     proj_line.set_3d_properties(np.zeros_like(y_real))
#     return line, proj_line
# ani = animation.FuncAnimation(fig2,update,frames = num_frames,interval =1000/fps, blit = False)
# plt.show()

'''
Identification procedure with BW procedure
'''
#calculate eta for the first file (no inhomogenity)
etas = np.zeros(frequencies.shape)
etasSp = np.zeros(frequencies.shape)
cbs = np.zeros(frequencies.shape)
for i, freq in enumerate(frequencies):
    meas = responsesP2[i,:]
    result = BW.Optimize(freq, Xpos, meas, BW.ObjectiveNormal)
    k_ = result.x[0]+1j*result.x[1]
    if (result.x[1]>0):
        print("imaginary part is positive at frequency "+str(freq)+ ", check formulation")
    E_ = beam.mu * freq*freq / beam.I /(pow(k_,4))
    etas[i] = E_.imag / E_.real
    k4 = k_**4
    etasSp[i] = k4.imag / k4.real
    omega = 2*np.pi*freq
    cbs[i] = pow(E_.real*beam.I*omega*omega/beam.mu,0.25)
    
# #calculate eta for the second file (inhomogenity)   
# etasInh = np.zeros(frequencies.shape)
# cbsInh = np.zeros(frequencies.shape)
# for i, freq in enumerate(frequencies):
#     meas = responsesP[i,:]
#     result = BW.Optimize(freq, Xpos, meas, BW.ObjectiveNormal)
#     k_ = result.x[0]+1j*result.x[1]
#     if (result.x[1]>0):
#         print("imaginary part is positive at frequency "+str(freq)+ ", check formulation")
#     E_ = beam.mu * freq*freq / beam.I /(pow(k_,4))
#     etasInh[i] = E_.imag / E_.real
#     omega = 2*np.pi*freq
#     cbsInh[i] = pow(E_.real*beam.I*omega*omega/beam.mu,0.25)

# etasInhFewPoints = np.zeros(frequencies.shape)
# cbsInhFewPoints = np.zeros(frequencies.shape)
# for i, freq in enumerate(frequencies):
#     meas = responsesP4[i,:]
#     result = BW.Optimize(freq, Xpos4, meas, BW.ObjectiveNormal)
#     k_ = result.x[0]+1j*result.x[1]
#     if (result.x[1]>0):
#         print("imaginary part is positive at frequency "+str(freq)+ ", check formulation")
#     E_ = beam.mu * freq*freq / beam.I /(pow(k_,4))
#     etasInhFewPoints[i] = E_.imag / E_.real
#     omega = 2*np.pi*freq
#     cbsInhFewPoints[i] = pow(E_.real*beam.I*omega*omega/beam.mu,0.25)


# fig4 = plt.figure(figsize=(8,6))
# ax4 = fig4.add_subplot(111)
# ax4.plot(frequencies,etas, label = label1)
# ax4.plot(frequencies,etasInh, label = label2 )
# ax4.plot(frequencies,etasInhFewPoints, label = label3)
# ax4.set_xlabel("Frequency [Hz]")
# ax4.set_ylabel("Loss Factor eta")
# plt.legend()
# plt.show()

# '''
# Identification procedure with OLF procedure
# '''
# h = 0.001
# b = 0.03
# rho = 7850
# pointDistance = 0.2
# eta_olf, cb_olf, E_olf, frequencies = olf.olf_improved(b, h, rho, file_path4, pointDistance, windowSize=1)
# eta_olf_inh, cb_olf_inh, E_olf_inh, frequenciesInh = olf.olf_improved(b, h, rho, file_path3, pointDistance, windowSize=1)
# fig5 = plt.figure(figsize=(8,6))
# ax5 = fig5.add_subplot(111)
# ax5.plot(frequencies,eta_olf, label = label1)
# ax5.plot(frequencies,eta_olf_inh, label = label2 )
# ax5.set_xlabel("Frequency [Hz]")
# ax5.set_ylabel("Loss Factor eta")
# plt.legend()
# plt.show()
