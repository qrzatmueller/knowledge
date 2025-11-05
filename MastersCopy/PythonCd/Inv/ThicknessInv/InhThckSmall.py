import sys
import BWoptimization as BW 
import SolutionSurface as SS
import os
import pandas as pd
import numpy as np
import olf
import matplotlib.pyplot as plt
'''
For thesis writing
first idea, do for small inhomogenities of increasing size 10 mm 20 mm 30 mm identification
see if there is a visible effect depending on the size
keep everything else constant, you wrote it on your notebook as a table. 
'''


# file_path2 = "Inv/06BWInh/RunCase_8/Result-Stage_3/res00.txt"
# frequencies2, responses2 = SS.openFEMresults(file_path2)
# etas, cbs, E = BW.IdentifyEta_YACS(file_path2, "BW", plot="Inhomo")


# file_path2 = "Inv/06BWInh/RunCase_9/Result-Stage_3/res00.txt"
# frequencies2, responses2 = SS.openFEMresults(file_path2)
# etas, cbs, E = BW.IdentifyEta_YACS(file_path2, "BW", plot="Inhomo")


# file_path2 = "Inv/06BWInh/RunCase_10/Result-Stage_3/res00.txt"
# frequencies2, responses2 = SS.openFEMresults(file_path2)
# etas, cbs, E = BW.IdentifyEta_YACS(file_path2, "BW", plot="Inhomo")

# file_path2 = "Inv/06BWInh/RunCase_13/Result-Stage_3/res00.txt"
# frequencies2, responses2 = SS.openFEMresults(file_path2)
# etas, cbs, E, k_ = BW.IdentifyEta_YACS(file_path2, "BW", plot="Inhomo")

# #no inhomogenity file
# file_path2 = "Inv/06BWInh/RunCase_11/Result-Stage_3/res00.txt"
# frequencies2, responses2 = SS.openFEMresults(file_path2)
# etas, cbs, E, k_ = BW.IdentifyEta_YACS(file_path2, "BW", plot="no homo")

'''
Second run with 1Gpa for lower stiffness
1Gpa turned out badly, but not because of inhomogenities, simply the wavelength is so small that the vibrations decay quicly and then 
there is nothing to identify
'''

# file_path2 = "Inv/06BWInh/RunCase_26/Result-Stage_3/res00.txt"
# frequencies2, responses2 = SS.openFEMresults(file_path2)
# etas, cbs, E, k_ = BW.IdentifyEta_YACS(file_path2, "BW", plot="Inhomo")


'''
remember what you changeeeed, the beam at identifyEta_yacs, I changed the E to match the simulations
'''

'''
to postprocess many quickly
'''
ns = np.array([54])
ks = np.zeros((577,ns.size),dtype=complex)
frequencies = np.zeros((1,1))
for i,  n in enumerate(ns): 
    file_path2 = "Inv/06BWInh/RunCase_" + str(n) + "/Result-Stage_3/res00.txt"
    frequencies2, responses2 = SS.openFEMresults(file_path2)
    frequencies = frequencies2
    etas, cbs, E, k_ = BW.IdentifyEta_YACS(file_path2, "BW", plot=("Inhomo RC" + str(i)))
    ks[:,i] = k_ 
# output_path = "Inv/06BWInh/ks51374849comp.dat"
# np.savetxt(output_path, np.column_stack([frequencies2, ks]), header='freq eta')

'''
to load the ones you want to see from .dat files and 
plot them 
'''

data = np.loadtxt("Inv/06BWInh/RunCase_40/Result-Stage_3/res00BWIDn.dat") #any file with the correct frequencies
frequencies = data[:,0] 
ns = np.arange(23,35)
ns = np.array([27,28,29,30,36,37,38,39]) #length comp
ns = np.arange(36,44) #position comparisson
ns = np.arange(40,48) #0.8 vs 1.2 at pos 6
ns = np.array([48,49,41,50,37])
# ns = np.array([34,35]) #control vs  homo at 0.1 eta
ns = np.array([50,54])

etas = np.zeros((frequencies.size, ns.size))
plt.figure()
plt.subplot(2,2,1)
for i, n in enumerate(ns):     
    # ------ Retrieve previously calculated data data ------
    file_path2 = "Inv/06BWInh/RunCase_" + str(n) + "/Result-Stage_3/res00BWIDn.dat"
    data = np.loadtxt(file_path2)
    etas[:,i] = data[:, 1]       # dimensionless
    
    plt.plot(frequencies,etas[:,i], label=str(n))
    plt.legend()
    
plt.show()

'''
to compare the ks directly for 34 and 35 (0.1 eta, without and with inh.)
I see that k'' increases while k' decreases a bit.
'''
#file_path2 = "Inv/06BWInh/ks3435comp.dat" #0.1 eta control vs inh
# file_path2 = "Inv/06BWInh/ks51374849comp.dat"
# data = np.loadtxt(file_path2,dtype=complex)
# frequencies= data[:,0]
# k1 = data[:,1]
# k2 = data[:,2]
# k3 = data[:,3]
# k4 = data[:,4]


# plt.subplot(2,2,2)
# plt.plot(frequencies.real,k1.real,label="normql")
# plt.plot(frequencies.real,k2.real,label="37")
# plt.plot(frequencies.real,k3.real,label="48")
# plt.plot(frequencies.real,k4.real,label="49")
# plt.legend()
# plt.subplot(2,2,3)
# plt.plot(frequencies.real,k1.imag,label="normql")
# plt.plot(frequencies.real,k2.imag,label="37")
# plt.plot(frequencies.real,k3.imag,label="48")
# plt.plot(frequencies.real,k4.imag,label="49")
# plt.legend()
# plt.show()

'''
Compare the actual signals
'''

