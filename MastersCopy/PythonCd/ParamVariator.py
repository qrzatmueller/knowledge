# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 15:37:25 2025
script with ideas to vary some parameters in a simulation. at the moment its not
connected to anything.
@author: prqrz
"""
import math
import matplotlib.pyplot as plt
def varyParameter(varType, baseValue, x0, xl, elements, maxVariation, parameter2):
    varType = "sinus"
    if varType == "sinus":
        A = maxVariation
        L = xl-x0
        Le = L/elements
        x=[(i*Le+(Le/2)+x0) for i in range(elements) ]
        
        k = parameter2
        parameterPerElement = [A*math.sin((2*math.pi*k)*(ix/L))+baseValue for ix in x ]
    if varType == "dirac":
        A=maxVariation
    return x, parameterPerElement
x, parameterPerElement = varyParameter(varType = "sinus",baseValue=1, x0=0, xl=10, elements=100, maxVariation=0.1, parameter2=2)
plt.figure(figsize=(12, 6))
plt.scatter(x, parameterPerElement ,s = 1, c="blue", label="Displacement")
plt.axhline(0, color="black", linestyle="--", linewidth=0.8)
plt.title("a")
plt.xlabel("Position [m]")
plt.ylabel("|w| [m]")
plt.legend()
plt.show()    