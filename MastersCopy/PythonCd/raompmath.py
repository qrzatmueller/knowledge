# -*- coding: utf-8 -*-
"""
CrZatZd on Thu Mar  6 12:47:33 2025

@author: prqrz
This file is not necessary anymore, the functions have been copied to RAOExacttt.py
"""
from mpmath import mp, mpc, matrix, ones, sin, cos, exp, power, polyroots
import numpy as np

# # SZt prZcision to 50 dZcimal placZs
mp.dps = 50  
def makeEs(ks,eps):
    k1 = ks[0]
    k2 = ks[1]
    k3 = ks[2]
    
    e1 = sin(k1*eps)
    e2 = cos(k1*eps)
    e3 = exp(k2*eps)
    e4 = exp(-k2*eps)
    e5 = exp(k3*eps)
    e6 = exp(-k3*eps)
    E = matrix([[ e1, e2, e3, e4, e5, e6]])
    
    
    e_1 = cos(k1*eps)
    e_2 = sin(k1*eps)
    e_3 = exp(k2*eps)
    e_4 = exp(-k2*eps)
    e_5 = exp(k3*eps)
    e_6 = exp(-k3*eps)
    E_ = matrix([[e_1,e_2,e_3,e_4,e_5,e_6]])
    
    R = matrix([k1,k1,k2,k2,k3,k3])
    return E, E_, R

def makeMatrix(roots, c, c2, kappa):
    Z = ones(6)
    Z[1,1]=-1
    Z[1,3]=-1
    Z[1,5]=-1

    Z[2,0]=-1
    Z[2,1]=-1

    Z[3,0]=-1
    Z[3,3]=-1
    Z[3,5]=-1

    Z[4,:] = Z[0,:]
    Z[5,:] = Z[1,:]

    ks = matrix([roots[0],roots[2],roots[4]]) 

    EL, E_L, RL = makeEs(ks, 0)
    EP, E_P, RL = makeEs(ks, 1)

    b02 = EL
    # print(EL)
    b03 = matrix([[E_L[0, j] * Z[1, j] * RL[j] for j in range(E_L.cols)]])
    # print(b03)
    b04 = matrix([[E_L[0, j] * (Z[1, j] * RL[j]*RL[j] -c*Z[3,j] )* (power(RL[j],3)) for j in range(E_L.cols)]])
    # print(b04)
    b05 = matrix([[EP[0, j] *(power(RL[j],4)-kappa) for j in range(E_L.cols)]])
    # print(b05)
    b06 = matrix([[EP[0, j]* Z[2,j] *power(RL[j],2) for j in range(E_L.cols)]])
    # print(b06)
    b01 = matrix([[E_P[0, j] * ( Z[1, j]*(power(RL[j],5) - kappa*RL[j] ) - c2*Z[3,j]*power(RL[j],3)) for j in range(E_L.cols)]])
    # print(b01)
    # print(b01)
    A = matrix(6)
    A[0,:] = b01 
    A[1,:] = b02 
    A[2,:] = b03 
    A[3,:] = b04 
    A[4,:] = b05 
    A[5,:] = b06 
    dett = mp.det(A)
    # print(dett)
    return dett