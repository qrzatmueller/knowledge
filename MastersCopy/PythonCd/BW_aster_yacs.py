'''
This file is used by the simulation automation in yacs_inHomo_Th
'''

import sys
import BWoptimization as BW 
import SolutionSurface as SS

if len(sys.argv)>1:
    input_path = sys.argv[1]
    print(f"Received path : {input_path}")
else:
    print("no path provided")
    input_path = "/home/pqrz/Documents/MasterThesis/mt_eq/CodeAsterModels/OLF/autoRes/res0.txt"
file_path2 = input_path
frequencies2, responses2 = SS.openFEMresults(file_path2)
etas, cbs,_ = BW.IdentifyEta_YACS(file_path2, "BW", plot="OLF with BWID")
