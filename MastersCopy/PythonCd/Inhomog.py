'''
In this file I will attempt to investigate the effects of having a thickness of one half of the value
for the second half of the BW beam. What happens to the wavefield, what happens to the identification. 
What happens if I only use the first half of measurement points, what about the second half, all?
'''

import os
import sys
import pandas as pd
import numpy as np
import olf
import BWoptimization as BW

current_file_path = os.path.abspath(__file__)
mt_eq_path = os.path.abspath(os.path.join(current_file_path, '..', '..', '..', '..', '..'))
file_path = os.path.join(mt_eq_path, 'CodeAsterModels/BW_Woking/BW_Inhomo_Th_Files/resF/RC2_HalfBeam.txt')

etas,cbs,Es = BW.IdentifyEta_YACS(file_path, "BW", "Plot", "HalfHalf")


