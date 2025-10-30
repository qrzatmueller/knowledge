# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 12:08:38 2025

@author: prqrz
"""
import numpy as np

x = np.linspace(0, 10, 100)
y = np.sin(x)

np.savetxt("data.dat", np.column_stack((x, y)), header="x y", comments="")


