# -*- coding: utf-8 -*-
"""
Created on Sat Mar 22 15:04:26 2025

@author: prqrz
This file uses the formulas from the book of cremer heckl to do forward and inverse calculation of properties
The results were not completely satisfactory

Nos quedamos pensando en como verificar que programe bien este ultima parte. 
Puedes probar con una simulacion de elementos finitos, o puedes tratar de 
replicar la grafica 3.23, ojo, chance tiene que tener una capa superior delgada
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
class EinfachBeläge:
    def __init__(self, d1, E1, rho1, d2, E2, rho2, eta2, a = -1):
        self.d1 = d1
        self.E1 = E1
        self.rho1 = rho1
        self.d2 = d2
        self.E2 = E2
        
        self.rho2 = rho2
        self.eta2 = eta2
        self.E2_ = E2*(1+1j*self.eta2)
        
        if (a == -1):
            self.a = (d1 + d2 )/2
        else:
            self.a = a
        
            
    def CalcEtaTot(self):
        etaTot = self.eta2*self.E2 * self.d2*self.a**2 / self.CalcBTot()
        etaTot2 = self.eta2*self.E2 * self.d2*self.a**2 / self.CalcBTotEx() 
        chatStiff,_ = self.ChatStiffness()
        etaTot3 = self.eta2*self.E2 * self.d2*self.a**2 / chatStiff
        return etaTot2, etaTot, etaTot3
    def CalcBTot(self):
        '''
        Uses an aproximation given by 3.76a, there is a longer formula at 3.103

        Returns
        -------
        BTot : TYPE
            Total stiffness of the 2 layered material

        '''
        BTot = self.E1 * self.d1**3 / 12 + self.E2 * self.d2 * self.a**2
        return BTot
    def CalcBTotEx(self):
        '''
        More precise than the one above, from formula 3.103

        Returns
        -------
        BTot : float
            Total Stiffness of the 2 layered material.

        '''
        BTot = self.E1*self.d1**3/12 + self.E2*self.d2**3/12 + ((self.E2*self.d2)/(1+ self.E2*self.d2/self.E1/self.d1))*((self.d1/2+self.d2/2)**2 - self.d1*self.d2/3)
        return BTot
    def CalcBTot_Ex(self):
        '''
        More precise than the one above, from formula 3.103

        Returns
        -------
        BTot : float
            Total Stiffness of the 2 layered material.

        '''
        BTot_ = self.E1*self.d1**3/12 + self.E2_*self.d2**3/12 + ((self.E2_*self.d2)/(1+ self.E2_*self.d2/self.E1/self.d1))*((self.d1/2+self.d2/2)**2 - self.d1*self.d2/3)
        return BTot_        
    def IdentifyEta2(self, etaTotMeas, Bmeas, E2):
        Eta2 = (etaTotMeas * Bmeas)/(E2*self.d2*self.a**2)
        return Eta2
    def AproxE2(self, Bmeas):
        E2 = (Bmeas - (self.E1*self.d1**3/12) )/(self.d2 * self.a**2)
        return E2
    def ObjectiveE2(self, k0, Bmeas):
        self.E2 = k0[0]
        residum = self.CalcBTotEx() - Bmeas
        # print(residum)
        return np.abs(residum)
    def IdentifyParameters(self, BMeas, etaMeas):
        E2Guess = self.AproxE2(BMeas)
        eta2Guess =  self.IdentifyEta2(etaMeas, BMeas, E2Guess)
        result = minimize(self.ObjectiveE2, E2Guess, args=(BMeas), method="Nelder-Mead", tol= 1e-5)
        E2 = result.x[0]
        eta2 =  self.IdentifyEta2(etaMeas, BMeas, E2)
        return E2, eta2, E2Guess, eta2Guess 
    def ChatStiffness(self):
        y = (self.E1*self.d1**2/2 + self.E2*self.d2*(self.d1 + self.d2/2))/(self.E1*self.d1 + self.E2*self.d2)
        y1 = self.d1/2 
        y2 = self.d1 + self.d2/2
        I1 = self.d1**3/12 + self.d1*(y1-y)**2
        I2 = self.d2**3/12 + self.d2*(y2-y)**2
        EIeff = self.E1*I1 + self.E2*I2
        return EIeff, y

class CremerSandwich():
    def __init__(self, d1, E1, rho1, d2, G2, rho2, eta2, d3, E3, rho3):
        self.B1 = (d1**3)/12 * E1
        self.B3 = (d3**3)/12 * E3
        self.d1 = d1
        self.E1 = E1
        self.rho1 = rho1
        self.d2 = d2
        self.G2 = G2
        self.rho2 = rho2
        self.eta2 = eta2
        self.d3 = d3
        self.E3 = E3
        self.rho3 = rho3
        self.a = self.d2 + (self.d1 + self.d3)/2
        self.m_ = self.rho1*self.d1 + self.rho2*self.d2 + self.rho3*self.d3
        self.h = self.Calc_Y()
        self.Btot = -1
        self.etaTot = -1
        self.g = -1
        self.Btot_ = -1
        self.g_ = -1
    def Calc_Y(self):
        h = 1 /( ((self.B1+self.B3)/self.a**2) * (1/self.E1/self.d1 + 1/self.E3/self.d3) )
        return h
    def Calc_g(self, k):
        g = (self.G2/self.d2/k**2) * ((1/self.E1/self.d1)+(1/self.E3/self.d3))
        return g
    def Calc_k(self, freq, Btot):
        om = freq * 2 * np.pi
        k = (om**2 * self.m_ / Btot )**0.25
        return k
    def Calc_B_(self, g_):
        p1 = (self.B1+self.B3)
        p2 = 1+ g_*self.h/(1+g_)
        B_ = p1*p2
        return B_
    def CalcParameters(self, freq):
        BGuess = self.B1+self.B3
        result = minimize(self.ObjectiveBs, BGuess, args= (freq), method = "Nelder-Mead" ,tol=1e-8)
        Btot = result.x[0]
        k = self.Calc_k(freq, Btot)
        g = self.Calc_g(k)
        self.Btot = Btot
        self.g = g
        self.g_ = g*(1+1j*self.eta2)
        self.Btot_ = self.Calc_B_(self.g_)
        return 
    def ObjectiveBs(self, BGuess, freq):
        # print("BGuess: ", BGuess)
        k = self.Calc_k(freq, BGuess)
        # print(k)
        g = self.Calc_g(k)
        g_ = g*(1+1j*self.eta2)
        B_ = self.Calc_B_(g_)
        dif = BGuess - B_.real
        # print("diff", np.abs(dif))
        return np.abs(dif)
    def Calc_EtaTot(self, freq):
        if (self.g == -1 or self.Btot == -1):
            self.CalcParameters(freq)
        g_ = self.g * (1+1j*self.eta2)
        etaTot = self.eta2* ((self.h*self.g)/(np.abs(1+g_)**2 + (self.g*self.h*(1+ self.g*(1+self.eta2**2)))))
        return etaTot
    def Identifyg_(self, measB_):
        A = (measB_) / (self.B1 + self.B3) -1
        g_ = -A / (A - self.h)
        return g_
    def IdentifyEta2(self, ident_g_):
        g_ = ident_g_
        eta2 = g_.imag / g_.real
        return eta2
    def IdentifyG2(self, k, g_):
        A = (1/(self.E1*self.d1)+1/(self.E3*self.d3))
        g = g_.real
        G2 = g / A * self.d2 * k**2
        return G2
    def IdentifyParameters(self, freq, measB_):
        g_ = self.Identifyg_(measB_)
        eta2 = self.IdentifyEta2(g_)
        om = freq*2*np.pi
        B = measB_.real
        k = (om**2 * self.m_ /B)**0.25
        G2 = self.IdentifyG2(k, g_)
        return G2, eta2, g_, k      
        
        
        
if __name__ == "__main__": 
    '''
    plot Bild 3.19
    '''
    
    # E1 = 206000000000
    # d1 = 0.001
    # E2 = E1*1e-1#1000000000
    # rho1 = 7850
    # rho2 = 2600
    # eta2 = 1
    # start = d1/10
    # stop = d1 * 20
    # n = 10000
    # d2 = np.linspace(start, stop, n)
    # EtaTot = np.zeros([n,2])
    # for i, d2i in enumerate(d2):
    #     system = EinfachBeläge(d1, E1, rho1, d2i, E2, rho2, eta2)
    #     EtaTot[i,:] = system.CalcEtaTot()
    # dRatio = d2/d1    
    # etaRatio = EtaTot / system.eta2 
    # fig = plt.figure(figsize=(8,6))
    # ax0 = fig.add_subplot()
    # plt.xscale("log")
    # plt.yscale("log")
    # ax0.plot(dRatio, etaRatio[:,0], label = "exact")
    # ax0.plot(dRatio, etaRatio[:,1], label = "approx")
    # ax0.set_xlabel("thickness ratio")
    # ax0.set_ylabel("eta ratio")
    # plt.title("bild 3.19")
    # plt.legend()
    # plt.show()
    '''
    Identification test
    '''
    
    # BMeas = system.CalcBTot()*1.2
    # etaMeas = 0.4
    # E2, eta2, E2guess, eta2guess = system.IdentifyParameters(BMeas, etaMeas)
    '''
    Forward calculation using eta2 and E2 as a given, to test the identification
    of eta2 and E2 with the calculated BMeas and etaMeas
    '''
    # E1 = 210000000000
    # d1 = 0.001
    # E2 = E1*1e-3
    # rho1 = 7850
    # rho2 = 2600
    # eta2 = 0.8
    # n = 10000
    # d2 = 0.005
    # system2 = EinfachBeläge(d1, E1, rho1, d2, E2, rho2, eta2)
    # BMeas = np.zeros([2,1])
    # BMeas[0] = system2.CalcBTotEx()
    # BMeas[1] = system2.CalcBTotEx()+12
    # BMeasGuess = system2.CalcBTot()
    # etaMeas, etaMeas2 = system2.CalcEtaTot()
    # etaMeasArr =  np.zeros([2,1])
    # etaMeasArr[0] = etaMeas2
    # etaMeasArr[1] = etaMeas2 + 0.12
    # system2.eta2 = eta2 * 0.75
    # system2.E2 = E2 * 1.13
    
    # system = EinfachBeläge(d1, E1, rho1, d2, E2 = 0, rho2 = rho2, eta2 = 0)
    # E2_iden, eta2_iden, E2guess_iden, eta2guess_iden = system2.IdentifyParameters(BMeas, etaMeasArr)
    '''
    Forward calculation of loss factor a three layered beam, which requires iterations
    uses formulas of CremerBook 3.6.2.2
    '''
    # E1 = 210000000000
    # d1 = 0.008
    # rho1 = 7850
    # E3 = E1#70000000000
    # d3 = 0.008
    # rho3 = rho1
    # G2 = 1000000000
    # d2 = 0.00005
    # rho2 = 2600
    # eta2 = 1.0 
    # system3 = CremerSandwich(d1, E1, rho1, d2, G2, rho2, eta2, d3, E3, rho3)
    # # freq = 1000
    # # etaTot = system3.Calc_EtaTot(freq)
    
    E1 = 210000000000
    d1 = 0.001
    rho1 = 7800
    E3 = 70000000000
    d3 = 0.00015
    rho3 = 2700
    G2 = 0.1e8
    d2 = 0.00145
    rho2 = 1650
    eta2 = 1.5
    system3 = CremerSandwich(d1, E1, rho1, d2, G2, rho2, eta2, d3, E3, rho3)
    freq = 613
    etaTot = system3.Calc_EtaTot(freq)
    print(etaTot)
    '''
    IDentification of eta2 and G2 upon previous knowledge of total loss factor
    and complex total stiffness. Using the values from example above, no simu-
    lation or measurement data yet
    '''
    # measB_ = system3.Btot_
    # G2iden, eta2iden, g_iden, kiden = system3.IdentifyParameters(freq, measB_)
    
    