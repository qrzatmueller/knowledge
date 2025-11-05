import BWoptimization as BW 
import SolutionSurface as SS
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

'''In this file I will identify the cb and eta values of the 
beam before and after the inhomogenity separately'''

ns = np.array([54])
ks = np.zeros((577,ns.size),dtype=complex)
frequencies = np.zeros((1,1))
for i,  n in enumerate(ns): 
    file_path2 = "Inv/06BWInh/RunCase_" + str(n) + "/Result-Stage_3/res00.txt"
    frequencies2, responses2 = SS.openFEMresults(file_path2)
    frequencies = frequencies2
    etas, cbs, E, k_ = BW.IdentifyEta_YACS(file_path2, "BW", plot=("Inhomo RC" + str(i)))
    ks[:,i] = k_ 