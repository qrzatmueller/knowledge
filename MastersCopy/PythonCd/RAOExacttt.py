# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 13:34:28 2025

@author: prqrz
I did many functions until it actually worked, in the end only a few are used
The first two matrices that I tried getting the determinant did not work, 
maybe not because of the matrices demselves but because I was not yet using 
mpmaths arbitrary precision. The first used exponential functions, the second
sinh and cosh functions. The third one, trigo functions from the paper of 
reference(6) in Raos paper with the example. All where ill conditioned and using
preconditionind did not help enough. The example 7.1 from Rao is replicated here
and the first, second and third mode are confirmed against the formulae he 
presented for the cantilever position. The ill conditioning of the matrices
that came from the boundary conditions proved render the calculation of the 
determinant useless (or not). Later, while using arbitrary precision, I learned
that the 0 of the determinant is really never found, but going through sign 
changes in the determinan means that you passed over the 0. I have not tested 
to see if with this approach and arbitrary precision I can find the same roots
with the first and second matrices. If they are equivalent I should!
Another important thing was the use of only 3 roots. Weirdly the reference 
(6) suggests using sqrt(z2,2j * rootsj), where z22 (for root1) is a -1. Why use
minus 1 here? idk. This could be answered if I understand his equivalent matrix
A, or test my sinh and cosh matrices with this new way of doing things. In my 
mind the 3 used roots should be positive wave numbers, right?
Additionally, one of the boundary conditions equations varies between that of 
RAO and his 6th reference..., the 6th reference is the one working now.
"""

import numpy as np
from sympy import symbols, solve, sqrt
import RKU
import sympy as sp
from scipy.optimize import minimize
from scipy.linalg import svd
from mpmath import mpf, mpc, diag, matrix, mp, svd, polyroots, sin, cos, exp, ones, power, lu_solve, residual, sqrt, norm, nprint, fabs, lu
import SolutionSurface as SS
import olf
import os
import matplotlib.pyplot as plt

iteration = 0
class SandwichBeam: 
    def __init__(self, h1,E1,rho1,h3,E3,rho3,h2,G2,eta2,rho2,b,L):
        self.h1 = h1
        self.E1 = E1
        self.rho1 = rho1
        self.h3 = h3
        self.E3 = E3
        self.rho3 = rho3
        self.h2 = h2
        self.G2 = G2
        self.rho2 = rho2 
        self.eta2 = eta2 
        self.b = b
        I1 = b*(h1*2)**3 / 12
        I3 = b*(h3*2)**3 / 12
        A1 = b*h1*2
        A2 = b*h2*2
        A3 = b*h3*2
        
        D = E1*I1+E3*I3
        m = rho1*A1 + rho2*A2 + rho3*A3 # mass per unit length of beam
        t0 = np.sqrt(m*L**4/D)
        G2_ = G2*(1+1j*eta2)
        g_ = (G2_*A2*L**2 / 4/h2**2) * ((E1*A1+E3*A3)/(E1*A1*E3*A3))  #todo L is also involved in g_ 
        c = 2*h2 + h1+ h3 
        Y = c**2/D * ((E1*A1*E3*A3)/(E1*A1+E3*A3))
        self.D = D
        self.t0 = t0
        self.G2_ = G2_
        self.g_ = g_ 
        self.Y = Y
        self.A2 = A2
        self.A1 = A1
        self.A3 = A3
        self.L = L
    def updateg_(self, eta2, G2):
        G2_ = G2*(1+1j*eta2)
        A2 = self.A2
        A1 = self.A1
        A3 = self.A3
        E1 = self.E1
        E3 = self.E3
        L = self.L
        h2 = self.h2
        g_ = (G2_*A2*L**2 / 4/h2**2) * ((E1*A1+E3*A3)/(E1*A1*E3*A3))
        self.G2 = G2
        self.eta2 = eta2 
        self.g_ = g_
        return 
    def Gfromg(self, g_):
        h2 = self.h2
        A2 = self.A2
        L = self.L
        E1 = self.E1
        A1 = self.A1
        E3 = self.E3
        A3 = self.A3
        G_ = g_*4*h2**2 /(A2*L**2)/(((E1*A1+E3*A3)/(E1*A1*E3*A3))) #todo check that L is the Xend-Xstart L
        
        return G_
    def CharEqn(self, Om_s):
        # Define symbols
        A, Omega, g, x = sp.symbols('A Omega g x', complex=True)
        
        # Define the cubic equation
        cubic_eq = x**3 - A*x**2 - Omega**2*x + Omega**2*g
        
        # Solve for x symbolically
        x_roots = sp.solve(cubic_eq, x)
        
        # Compute k values (square roots of x)
        k_roots = []
        for x_i in x_roots:
            sqrt_xi = sp.sqrt(x_i)  # Principal square root
            k_roots.append(sqrt_xi)
            k_roots.append(-sqrt_xi)  # Include negative sqrt
        
        # Substitute values for A, Omega, and g
        A_value = self.g_*(1+ self.Y) 
        g_value = self.g_
        
        # Create substitutions dictionary
        Om_ = Om_s[0] + Om_s[1]*1j
        subs_dict = {A: A_value, g: g_value, Omega: Om_}
        
        # Compute numerical values for k
        numerical_k = [k.subs(subs_dict).evalf() for k in k_roots]
        numerical_ks = [complex(ks.evalf()) for ks in numerical_k]
        nprint(polyroots([1,-1*A_value,-1*Om_*Om_,Om_*Om_*g_value],maxsteps=1000))
        return np.array(numerical_ks)
    
    def getDet(self, Om_):
        roots = self.CharEqn(Om_)
        det = self.DetFunc(Om_, roots)
        absdet = np.abs(det)
        return absdet
    
    def DetFunc(self, Om_s, k_s):
        '''
        Uses first developed matrix and symbolic notation. Didnt work.

        Parameters
        ----------
        Om_s : TYPE
            DESCRIPTION.
        k_s : TYPE
            DESCRIPTION.

        Returns
        -------
        determinant : TYPE
            DESCRIPTION.

        '''
        # Step 1: Define symbolic unknowns
        k1, k2, k3, k4, k5, k6, Om_, g_ = sp.symbols('k1 k2 k3 k4 k5 k6 Om_, g_', complex = True)
        t0, Y = sp.symbols('t0, Y') 
        k = sp.symbols('k')
        f_x2 = k**2 * sp.exp(k)
        row2_values = [f_x2.subs(sp.symbols('k'), ks) for ks in [k1, k2, k3, k4, k5, k6]]
        f_x3 = sp.exp(k)*(k**4  - (Om_**2 / t0 **2))
        row3_values = [f_x3.subs(sp.symbols('k'), ks) for ks in [k1, k2, k3, k4, k5, k6]]
        f_x4 = k**4  - g_*Y*k**3
        row4_values = [f_x4.subs(sp.symbols('k'), ks) for ks in [k1, k2, k3, k4, k5, k6]]
        f_x5 = sp.exp(k)*(k**5  - g_*(1+Y)*k**3 -(Om_**2 / t0 **2)*k)
        row5_values = [f_x5.subs(sp.symbols('k'), ks) for ks in [k1, k2, k3, k4, k5, k6]]
    
        # Step 2: Define a 6x6 matrix with unknowns
        M = sp.Matrix([
            [1, 1, 1, 1, 1, 1],
            [k1, k2, k3, k4, k5, k6],
            row2_values,
            row3_values,
            row4_values,
            row5_values
        ])
        
        # Step 3: Compute determinant symbolically
        #det_expr = M.det().simplify()
        
        M_func = sp.lambdify((k1, k2, k3, k4, k5, k6, Om_, g_, t0, Y), M, "numpy")
        #Normalizing rows before calculating determinant
        # Call it with specific values
        Om_ev = Om_s[0] + Om_s[1]*1j
        M_numpy_fast = np.array(M_func(k_s[0],k_s[1],k_s[2],k_s[3],k_s[4],k_s[5], Om_ev, self.g_, self.t0, self.Y), dtype=np.complex_)
        mp.dps = 100
        mpmath_matrix = mp.matrix([[mpc(m.real, m.imag) for m in row] for row in M_numpy_fast])
        det_A = mp.det(mpmath_matrix)
        print(det_A)
        singular_values = svd(M_numpy_fast)
        
        scaling_factors = [np.max(np.abs(M_numpy_fast[i, :])) for i in [-5,-4,-3,-2, -1]]
        A_normalized = normalize_last_rows(M_numpy_fast)
        A_norm2 = normalize_rows(M_numpy_fast)
        A_norm3 = normalize_columns(M_numpy_fast)
        A_norm4 = normalize_columns(A_norm2)
        cond = np.linalg.cond(M_numpy_fast)
        cond2 = np.linalg.cond(A_normalized)
        cond3 = np.linalg.cond(A_norm2)
        cond4 = np.linalg.cond(A_norm3)
        cond5 = np.linalg.cond(A_norm4)
        det_scaled = np.linalg.det(A_normalized)
        det_original = det_scaled * np.prod(scaling_factors)
        det_scaled2 = np.linalg.det(A_norm2)

        determinant = np.linalg.det(M_numpy_fast)
        
        return determinant
    def DetFuncTrig(self, Om_s, k_s):
        # Step 1: Define symbolic unknowns
        k1, k2, k3, k4, k5, k6, Om_, g_ = sp.symbols('k1 k2 k3 k4 k5 k6 Om_, g_', complex = True)
        t0, Y, k = sp.symbols('t0 Y k') 
        c1, c2, c3 = sp.symbols('c1 c2 c3', complex = True)
        #where c1 is -g_Y, c2 = Om_**2/t0**2, and c3 is g_*(Y+1)
    
        f_x21 = k**4 
        f_x22 = c1*k**3
        row22_values = [f_x21.subs(sp.symbols('k'), ks) if i%2 == 0 else f_x22.subs(sp.symbols('k'),ks) for i, ks in enumerate([k1, k1, k3, k3, k5, k5])]
        f_x31 = k**2*sp.cosh(k)
        f_x32 = k**2*sp.sinh(k)
        row33_values = [f_x31.subs(sp.symbols('k'), ks) if i%2 == 0 else f_x32.subs(sp.symbols('k'),ks) for i, ks in enumerate([k1, k1, k3, k3, k5, k5])]
        f_x41 = sp.cosh(k)*(k**4-(c2))
        f_x42 = sp.sinh(k)*(k**4-(c2))
        row44_values =  [f_x41.subs(sp.symbols('k'), ks) if i%2 == 0 else f_x42.subs(sp.symbols('k'),ks) for i, ks in enumerate([k1, k1, k3, k3, k5, k5])]
        f_x51 = sp.sinh(k)*(k**5- c3*k**3 -(c2)*k)
        f_x52 = sp.cosh(k)*(k**5- c3*k**3 -(c2)*k)
        row55_values = [f_x51.subs(sp.symbols('k'), ks) if i%2 == 0 else f_x52.subs(sp.symbols('k'),ks) for i, ks in enumerate([k1, k1, k3, k3, k5, k5])]
        
        M3 = sp.Matrix([
            [1, 0, 1, 0, 1, 0],
            [0, k1, 0, k3, 0, k5],
            row22_values,
            row33_values,
            row44_values,
            row55_values
        ])

        
        # Step 3: Compute determinant symbolically
        #det_expr = M.det().simplify()
        
        M_func = sp.lambdify((k1, k2, k3, k4, k5, k6, Om_, c1, c2, c3), M3, "numpy")

        #Normalizing rows before calculating determinant
            
        # Call it with specific values
        Om_ev = Om_s[0] + Om_s[1]*1j
        c1 = -1 * self.g_ * self.Y
        c2 = Om_ev**2 / (self.t0)**2
        c3 = self.g_ * (self.Y + 1)
        M_numpy_fast = np.array(M_func(k_s[0],k_s[0],k_s[2],k_s[2],k_s[4],k_s[4], Om_ev, c1, c2, c3), dtype=np.complex_)
        singular_values = svd(M_numpy_fast)
        
        #with crazy precision
        mp.dps = 100
        mpmath_matrix = mp.matrix([[mpc(m.real, m.imag) for m in row] for row in M_numpy_fast])
        det_A = mp.det(mpmath_matrix)
        print(det_A)
        scaling_factors = [np.max(np.abs(M_numpy_fast[i, :])) for i in [-3,-2, -1]]
        A_normalized = normalize_last_rows(M_numpy_fast)
        A_norm2 = normalize_rows(M_numpy_fast)
        A_norm3 = normalize_columns(M_numpy_fast)
        A_norm4 = normalize_columns(A_norm2)
        cond = np.linalg.cond(M_numpy_fast)
        cond2 = np.linalg.cond(A_normalized)
        cond3 = np.linalg.cond(A_norm2)
        cond4 = np.linalg.cond(A_norm3)
        cond5 = np.linalg.cond(A_norm4)
        det_scaled = np.linalg.det(A_normalized)
        det_original = det_scaled * np.prod(scaling_factors)
        det_scaled2 = np.linalg.det(A_norm2)
        det_scaled3 = np.linalg.det(A_norm4)
        determinant = np.linalg.det(M_numpy_fast)
        
        return determinant

    def checking(self, OmL):
        roots = self.CharEqn(OmL)
        Om_ = OmL[0]+ 1j*OmL[1]
        zeroes= np.zeros(6,dtype = complex)
        for i, root in enumerate(roots): 
            result = -root**6 + self.g_*(1+self.Y) * root**4 + Om_**2*(root**2-self.g_)
            zeroes[i] = result
        return zeroes
    def findModes(self, guessFreq):
        fSol = guessFreq
        EtaSol = 0.01
        pSol = fSol*2*np.pi 
        OmSol = pSol * self.t0 
        Om_Sol = OmSol*(1+1j*EtaSol)**0.5 
        Om_ = [Om_Sol.real, Om_Sol.imag]
        
        c1 = self.g_*self.Y
        c2 = self.g_ *(1+self.Y)
        kappa = Om_Sol**2
        
        # make precision variables
        mp.dps = 50 #desired precision
        # mc1 = mpc(c1)
        # mc2 = mpc(c2)
        # mg = mpc(self.g_)
        mc1 = c1
        mc2 = c2
        mg = self.g_
        x = Optimize(fSol, EtaSol, self.t0, mc1, mc2, mg, Objectiveplz)
        
        optFreq = x.x[0]
        optEta = x.x[1] 
        print("Optimal:", x.x)
        return x
    def ForwardCalc(self, result_path, case, plot = "", remark = "", Xposi = [], resp = [0,1,-1]):
        file_path = result_path
        frequencies, responses = SS.openFEMresults(file_path)
        if (case == "BW"):
            #Xpos for usual BW, I have not retested this
            start = 0.170
            step = 0.020
            x_size = 21
            stop = start + ((x_size-1) * step)
            Xpos = np.linspace(start, stop, x_size)
             # for initial guess use steel beam
        if (case == "OLF"):    
            #Xpos for OLF, and guessing beam: 
            start = 0.170
            step = 0.01 #it cannot work like this because the steps are not regelmäßig
            x_size = 12 #meaning one step but 2 datapoints
            stop = start + ((x_size-1) * step)
            Xpos = np.linspace(start, stop, x_size)
            Xpos2 = np.linspace(start+0.2, stop+0.2, x_size)
            Xpos = np.concatenate((Xpos, Xpos2), axis = 0)
             # for initial guess use steel beam
        if(case == "Custom"):
            Xpos = Xposi
            
        betas = np.zeros(frequencies.shape)
        G_s = np.zeros(frequencies.shape, dtype="complex")    
        for i, freq in enumerate(frequencies):
            i = 19 #todo delete this with the break at the end
            freq = frequencies[i]
            meas = responses[i,resp[0]:resp[1]:resp[2]]
            femETA = 0.114 #0.048 #global eta, then try with the exact one
            # make precision variables
            mp.dps = 15 #desired precision
            t0 = self.t0 
            fSol = freq 
            pSol = fSol*2*np.pi 
            OmSol = pSol * t0 
            Om_Sol = OmSol*(1+1j*femETA)**0.5
            mc2 = self.g_ *(1+self.Y)
            mg = self.g_
            coeffs = [1, -mc2, -Om_Sol**2, Om_Sol**2 * mg]
            # print("Coeffs1", coeffs)
            x_roots = np.roots(coeffs)
            k_roots = np.column_stack((np.sqrt(x_roots), -np.sqrt(x_roots))).ravel()
            
            # k_roots = k_roots * 1j #todo delete this just added to test if an M-mj configuration is better
            # k_roots[2] = k_roots[2].real - k_roots[2].imag*1j #todo this whole manipulation of
            # k_roots[3] = k_roots[3].real - k_roots[3].imag*1j #the wavenumbers did help, but why?? and why was it wrong???
            # k_roots[4] = k_roots[4].real - k_roots[4].imag*1j
            # k_roots[5] = k_roots[5].real - k_roots[5].imag*1j
            
            # k_roots = k_roots * 1j #todo delete this just added to test if an M-mj configuration is better
            k_roots[2] = k_roots[2].real - k_roots[2].imag*1j #todo this whole manipulation of
            k_roots[3] = k_roots[3].real - k_roots[3].imag*1j #the wavenumbers did help, but why?? and why was it wrong???
            k_roots[4] = k_roots[4].real - k_roots[4].imag*1j
            k_roots[5] = k_roots[5].real - k_roots[5].imag*1j
            
            print("forward roots", k_roots)
            if (abs(k_roots[0].real) > 50): #todo this 10 is completely arbitrary
                k_roots = k_roots[2:6]
                
                # print(k_roots)
                A = getSmallNumpyAMatrix(k_roots, Xpos)
                A2 = getSmallAMatrix(k_roots, Xpos)
            elif (k_roots[2].real > 50):
                k_roots = k_roots[0:2]
                print("Super Reduced roots", k_roots)
                A = getTinyAMatrix(k_roots,Xpos)
            else: 
                k_roots[0] = k_roots[0]*1j
                A = getAMatrix(k_roots, Xpos) #matrix from paper
            x2, residualNorm, _, _ = np.linalg.lstsq(A, meas, rcond=None)
            if residualNorm.size == 0:
                U, S, Vt = np.linalg.svd(A)
                dependent_col_index = np.argmax(np.abs(Vt[-1]))
                print("removing dependent column at index:", dependent_col_index)
                A_reduced = np.delete(A, dependent_col_index, axis=1)
                x2, residualNorm, rank, s = np.linalg.lstsq(A_reduced, meas, rcond=None)
                A = A_reduced
            
            # x = lu_solve(A, meas) #not working because of really large numbers (numerically singular matrix)
            # res = residual(A, x, meas)
            # l2norm = norm(res, p=2)
            print("residualNorm", residualNorm)
            curveMp = A@x2
            curve = np.array(curveMp.tolist(), dtype=complex)
            fig = plt.figure(figsize=(8,6))
            ax0 = fig.add_subplot(221)
            ax0.plot(Xpos, curve.real, label = "curve real")
            ax0.plot(Xpos, meas.real, label = "meas real")
            ax0.set_xlabel("Position (x)")
            ax0.set_ylabel("Real Part")
            plt.title(str(freq)+ "Hz")
            plt.legend()
            
            ax0 = fig.add_subplot(222)
            
            ax0.plot(Xpos, curve.imag, label = "curve imag")
            ax0.plot(Xpos, meas.imag, label = "meas imag")
            ax0.set_xlabel("Position (x)")
            ax0.set_ylabel("Imag Part")
            plt.title(str(freq)+ "Hz")
            plt.legend()
            
            ax0 = fig.add_subplot(223)
            
            ax0.plot(Xpos, np.angle(curve) , label = "curve phase")
            ax0.plot(Xpos, np.angle(meas), label = "meas phase")
            ax0.set_xlabel("Position (x)")
            ax0.set_ylabel("Phase")
            plt.title(str(freq)+ "Hz")
            plt.legend()
            
            ax0 = fig.add_subplot(224)
            
            ax0.plot(Xpos, np.abs(curve), label = "curve mag")
            ax0.plot(Xpos, np.abs(meas), label = "meas mag")
            ax0.set_xlabel("Position (x)")
            ax0.set_ylabel("Magnitude")
            plt.title(str(freq)+ "Hz")
            plt.legend()
            plt.show()
            
            break

        return
    def ForwardWave(self, result_path, case, plot = "", remark = "", Xposi = [], resp = [0,1,-1]):
        file_path = result_path
        frequencies, responses = SS.openFEMresults(file_path)
        if (case == "BW"):
            #Xpos for usual BW, I have not retested this
            start = 0.170
            step = 0.020
            x_size = 21
            stop = start + ((x_size-1) * step)
            Xpos = np.linspace(start, stop, x_size)
             # for initial guess use steel beam
        if (case == "OLF"):    
            #Xpos for OLF, and guessing beam: 
            start = 0.170
            step = 0.01 #it cannot work like this because the steps are not regelmäßig
            x_size = 12 #meaning one step but 2 datapoints
            stop = start + ((x_size-1) * step)
            Xpos = np.linspace(start, stop, x_size)
            Xpos2 = np.linspace(start+0.2, stop+0.2, x_size)
            Xpos = np.concatenate((Xpos, Xpos2), axis = 0)
             # for initial guess use steel beam
        if(case == "Custom"):
            Xpos = Xposi
            
        betas = np.zeros(frequencies.shape)
        G_s = np.zeros(frequencies.shape, dtype="complex")    
        for i, freq in enumerate(frequencies):
            i = 19 #todo delete this with the break at the end
            freq = frequencies[i]
            meas = responses[i,resp[0]:resp[1]:resp[2]]
            femETA = 0.228 #0.114 #0.048 #global eta, then try with the exact one
            # make precision variables
            mp.dps = 15 #desired precision
            t0 = self.t0 
            fSol = freq 
            pSol = fSol*2*np.pi 
            OmSol = pSol * t0 
            Om_Sol = OmSol*(1+1j*femETA)**0.5
            mc2 = self.g_ *(1+self.Y)
            mg = self.g_
            coeffs = [1, -mc2, -Om_Sol**2, Om_Sol**2 * mg]
            # print("Coeffs1", coeffs)
            x_roots = np.roots(coeffs)
            k_roots = np.column_stack((np.sqrt(x_roots), -np.sqrt(x_roots))).ravel()
            
            # k_roots = k_roots * 1j #todo delete this just added to test if an M-mj configuration is better
            # k_roots[2] = k_roots[2].real - k_roots[2].imag*1j #todo this whole manipulation of
            # k_roots[3] = k_roots[3].real - k_roots[3].imag*1j #the wavenumbers did help, but why?? and why was it wrong???
            # k_roots[4] = k_roots[4].real - k_roots[4].imag*1j
            # k_roots[5] = k_roots[5].real - k_roots[5].imag*1j
            
            # k_roots = k_roots * 1j #todo delete this just added to test if an M-mj configuration is better
            k_roots[2] = k_roots[2].real - k_roots[2].imag*1j #todo this whole manipulation of
            k_roots[3] = k_roots[3].real - k_roots[3].imag*1j #the wavenumbers did help, but why?? and why was it wrong???
            k_roots[4] = k_roots[4].real - k_roots[4].imag*1j
            k_roots[5] = k_roots[5].real - k_roots[5].imag*1j
            
            print("forward roots", k_roots)
            
            if (abs(k_roots[0].real) > 50): #todo this 10 is completely arbitrary
                k_roots = k_roots[2:6]
                
                # print(k_roots)
                A = getSmallNumpyAMatrix(k_roots, Xpos)
                A2 = getSmallAMatrix(k_roots, Xpos)
            elif (k_roots[2].real > 50):
                k_roots = k_roots[0:2]
                print("Super Reduced roots", k_roots)
                A = getTinyAMatrix(k_roots,Xpos)
            else: 
                k_roots[0] = k_roots[0]*1j
                A = getAMatrix(k_roots, Xpos) #matrix from paper
            curve = 4e-05*np.exp(k_roots[1]*Xpos)
            fig = plt.figure(figsize=(8,6))
            ax0 = fig.add_subplot(221)
            ax0.plot(Xpos, curve.real, label = "curve real")
            ax0.plot(Xpos, meas.real, label = "meas real")
            ax0.set_xlabel("Position (x)")
            ax0.set_ylabel("Real Part")
            plt.title(str(freq)+ "Hz")
            plt.legend()
            
            ax0 = fig.add_subplot(222)
            
            ax0.plot(Xpos, curve.imag, label = "curve imag")
            ax0.plot(Xpos, meas.imag, label = "meas imag")
            ax0.set_xlabel("Position (x)")
            ax0.set_ylabel("Imag Part")
            plt.title(str(freq)+ "Hz")
            plt.legend()
            
            ax0 = fig.add_subplot(223)
            
            ax0.plot(Xpos, np.angle(curve) , label = "curve phase")
            ax0.plot(Xpos, np.angle(meas), label = "meas phase")
            ax0.set_xlabel("Position (x)")
            ax0.set_ylabel("Phase")
            plt.title(str(freq)+ "Hz")
            plt.legend()
            
            ax0 = fig.add_subplot(224)
            
            ax0.plot(Xpos, np.abs(curve), label = "curve mag")
            ax0.plot(Xpos, np.abs(meas), label = "meas mag")
            ax0.set_xlabel("Position (x)")
            ax0.set_ylabel("Magnitude")
            plt.title(str(freq)+ "Hz")
            plt.legend()
            plt.show()
            
            break

        return
    def IdentifyG_(self, result_path, case, plot = "", remark = "", Xposi = [], resp = [0,1,-1]):
        file_path = result_path
        frequencies, responses = SS.openFEMresults(file_path)
        if (case == "BW"):
            #Xpos for usual BW, I have not retested this
            start = 0.170
            step = 0.020
            x_size = 21
            stop = start + ((x_size-1) * step)
            Xpos = np.linspace(start, stop, x_size)
             # for initial guess use steel beam
        if (case == "OLF"):    
            #Xpos for OLF, and guessing beam: 
            start = 0.170
            step = 0.01 #it cannot work like this because the steps are not regelmäßig
            x_size = 12 #meaning one step but 2 datapoints
            stop = start + ((x_size-1) * step)
            Xpos = np.linspace(start, stop, x_size)
            Xpos2 = np.linspace(start+0.2, stop+0.2, x_size)
            Xpos = np.concatenate((Xpos, Xpos2), axis = 0)
             # for initial guess use steel beam
        if(case == "Custom"):
            Xpos = Xposi
            
        betas = np.zeros(frequencies.shape)
        etas = np.zeros(frequencies.shape)
        G_s = np.zeros(frequencies.shape, dtype="complex")    
        for i, freq in enumerate(frequencies):
            i = 19 #todo delete this with the break at the end
            freq = frequencies[i]
            meas = responses[i,resp[0]:resp[1]:resp[2]]
            # result3 = self.Optimize4(freq, Xpos, meas, self.Objective4)
            result = self.Optimize5(freq, Xpos, meas, self.Objective5)
            # result2 = self.Optimize3(freq, Xpos, meas, self.Objective3)
            # result = self.Optimize2(freq, Xpos, meas, self.Objective2)
            
            self.updateg_(eta2=result.x[2], G2=result.x[1])
            if (result.x[2]<0):
                print("imaginary part is negative at frequency "+str(freq)+ ", check formulation")
            G_ = self.Gfromg(self.g_) #insert formula for getting G from g
            betas[i] = G_.imag / G_.real
            etas[i] = result.x[0]
            G_s[i] = G_ #todo insert formula that makes sense 
            print("freq: ", freq)
            print("Success: ", result.success)
            print("G_", G_)
            print("Eta2, ", betas[i])
            break

        if (plot != ""):
            olf.save_scatter_plot(file_path,"BWID"+remark, frequencies, betas, plot)
        return etas, G_s ,betas
    
    def Optimize2(self,freq, Xpos, meas, objective):
        guessEta = 0.5 #0.048 #global eta, then try with the exact one
        k0 = [guessEta, self.g_.real,self.g_.imag]#use the guessed values
             #given the algorithm of Nelder-Mead (uses 5% change for vertices of Simplex), using a higher eta is better than close to 0, to give it a big area to look
        c1 = self.g_*self.Y
        c2 = self.g_ *(1+self.Y)
        # make precision variables
        mp.dps = 15 #desired precision
        mc1 = (c1) #I removed all the mpc variables.
        mc2 = (c2)
        mg = (self.g_)
        t0 = self.t0 #todo think about this, Do I want to specify a length? can I calculate without normalizing time? what happens to the boundary conditions
        result = minimize(objective, k0, args =(freq, meas, Xpos, t0, mc1, mc2, mg),method = 'Nelder-Mead')
        return result
    
    def Objective2(self, k0, freq, meas, Xpos, t0, mc1, mc2, mg):
        #print(k0)
        fSol = freq
        EtaSol = (k0[0])
        g_real = (k0[1])
        g_imag = (k0[2])
        penalty = 0.0
        penalty = penalty + 1e4 if k0[0] < 0 else penalty +0.0
        penalty = penalty + 1e4 if k0[1] < 0 else penalty +0.0
        penalty = penalty + 1e4 if k0[2] < 0 else penalty +0.0
        g_ = g_real + 1j*g_imag
        fSol = freq 
        pSol = fSol*2*np.pi 
        OmSol = pSol * t0 
        Om_Sol = OmSol*(1+1j*EtaSol)**0.5

        ''' I need to calculate the mcs and mgs here!!!!'''
        mc2 = g_ *(1+self.Y)
        mg = g_

        # kappa = Om_Sol**2
        
        # mkappa = (kappa)
        # mroots = (polyroots([1,-mc2,-mkappa, mkappa*mg],maxsteps=100))
        # mrootsd = matrix([mroots[0],mroots[0], mroots[1], mroots[1], mroots[2],mroots[2]])
        # mrootsd[0] = mrootsd[0]*-1 #according to the paper this needs this change of sign...
        # for i,r in enumerate(mrootsd): 
        #     mrootsd[i] = sqrt((mrootsd[i]))
        coeffs = [1, -mc2, -Om_Sol**2, Om_Sol**2 * mg]
        # print("Coeffs2", coeffs)
        x_roots = np.roots(coeffs)
        k_roots = np.column_stack((np.sqrt(x_roots), -np.sqrt(x_roots))).ravel()
        if (abs(k_roots[0].real) > 10): #todo this 10 is completely arbitrary
            k_roots = k_roots[2:6]
            # print(k_roots)
            A = getSmallAMatrix(k_roots, Xpos)
        else: 
            k_roots[0] = k_roots[0]*1j
            A = getAMatrix(k_roots, Xpos) #matrix from paper
        
        x = lu_solve(A, meas) #not working because of really large numbers (numerically singular matrix)
        res = residual(A, x, meas)
        curveMp = A*x
        curve = np.array(curveMp.tolist(), dtype=complex)
        global iteration
        if (iteration%100 == 0 ):
            updateGraph(Xpos, meas, curve)
        iteration = iteration +1
        # C = A[0:-40:1,:] #reduced A matrix until the last forty
        # x2 = lu_solve(C, meas[0:-40:1])
        # # res = C * x2 - meas[0:-40:1]
        # res = residual(C, x2, meas[0:-40:1])

        # B = np.array(A.tolist(), dtype=complex) #A matrix but as numpy array
        # column_norms = np.linalg.norm(B,axis=0)
        # Bnorm  = B / np.maximum(column_norms, 1e-10)
        # x, res, _, _ = np.linalg.lstsq(Bnorm, meas, rcond=None)
        # res2 = Bnorm@x - meas

        # B2 = B[0:-40:1, : ] #reducdd numpy array
        # x1, res2 = np.linalg.lstsq(B2, meas[0:-40:1], rcond = None)
        # # res = residual(B2, x1, meas)
        real_parts = [z.real for z in res]
        imaginary_parts = [z.imag for z in res]
        l2normReal = norm(real_parts,p=2)
        l2normImag = norm(imaginary_parts, p=2)
        l2norm = (norm(res, p=2))
        sepNorm = 10*l2normReal+l2normImag
        print("normal norm", l2norm)
        print("Separate norm", sepNorm)

        return sepNorm+penalty
    def Optimize3(self,freq, Xpos, meas, objective):
        guessEta = 0.048 #global eta, then try with the exact one
        k0 = [guessEta, self.G2,self.eta2]#use the guessed values
             #given the algorithm of Nelder-Mead (uses 5% change for vertices of Simplex), using a higher eta is better than close to 0, to give it a big area to look
        # make precision variables
        mp.dps = 20 #desired precision
        result = minimize(objective, k0, args =(freq, meas, Xpos),method = 'Nelder-Mead')
        return result
    def Objective3(self, k0, freq, meas, Xpos):
        #print(k0)
        fSol = freq
        EtaSol = (k0[0])
        G2 = (k0[1])
        eta2 = (k0[2])
        self.updateg_(eta2, G2)
        c2 = self.g_ *(1+self.Y)
        
        penalty = 0.0
        penalty = penalty + 1e4 if k0[0] < 0 else penalty +0.0
        penalty = penalty + 1e4 if k0[1] < 0 else penalty +0.0
        penalty = penalty + 1e4 if k0[2] < 0 else penalty +0.0
        penalty = penalty + 1e4 if k0[2] > 3 else penalty +0.0
        fSol = freq 
        pSol = fSol*2*np.pi 
        OmSol = pSol * self.t0 
        Om_Sol = OmSol*(1+1j*EtaSol)**0.5

        ''' I need to calculate the mcs and mgs here!!!!'''
        mc2 = self.g_ *(1+self.Y)
        mg = self.g_

        # kappa = Om_Sol**2
        
        # mkappa = (kappa)
        # mroots = (polyroots([1,-mc2,-mkappa, mkappa*mg],maxsteps=100))
        # mrootsd = matrix([mroots[0],mroots[0], mroots[1], mroots[1], mroots[2],mroots[2]])
        # mrootsd[0] = mrootsd[0]*-1 #according to the paper this needs this change of sign...
        # for i,r in enumerate(mrootsd): 
        #     mrootsd[i] = sqrt((mrootsd[i]))
        coeffs = [1, -mc2, -Om_Sol**2, Om_Sol**2 * mg]
        # print("Coeffs2", coeffs)
        x_roots = np.roots(coeffs)
        k_roots = np.column_stack((np.sqrt(x_roots), -np.sqrt(x_roots))).ravel()
        if (abs(k_roots[0].real) > 10): #todo this 10 is completely arbitrary
            k_roots = k_roots[2:6]
            # print(k_roots)
            A = getSmallAMatrix(k_roots, Xpos)
            A2 = getNumpyAMatrix(k_roots, Xpos)
        else: 
            k_roots[0] = k_roots[0]*1j
            A = getAMatrix(k_roots, Xpos) #matrix from paper
        column_norms = np.linalg.norm(A2,axis=0)#this normalization is new because from olf measurements it wasnt working without it, maybe because of big x?
        Anorm  = A2 / np.maximum(column_norms, 1e-10)
        # Solve least squares for x
        x2, residualNorm, _, _ = np.linalg.lstsq(Anorm, meas, rcond=None)
        
        # x = lu_solve(A, meas) #not working because of really large numbers (numerically singular matrix)
        # res = residual(A, x, meas)
        # curveMp = A*x
        # curve = np.array(curveMp.tolist(), dtype=complex)
        # global iteration
        # if (iteration%20 == 0 ):
        #     print(k0)
        #     updateGraph(Xpos, meas, curve)
        # iteration = iteration +1
        # # C = A[0:-40:1,:] #reduced A matrix until the last forty
        # # x2 = lu_solve(C, meas[0:-40:1])
        # # # res = C * x2 - meas[0:-40:1]
        # # res = residual(C, x2, meas[0:-40:1])

        # # B = np.array(A.tolist(), dtype=complex) #A matrix but as numpy array
        # # column_norms = np.linalg.norm(B,axis=0)
        # # Bnorm  = B / np.maximum(column_norms, 1e-10)
        # # x, res, _, _ = np.linalg.lstsq(Bnorm, meas, rcond=None)
        # # res2 = Bnorm@x - meas

        # # B2 = B[0:-40:1, : ] #reducdd numpy array
        # # x1, res2 = np.linalg.lstsq(B2, meas[0:-40:1], rcond = None)
        # # # res = residual(B2, x1, meas)
        # real_parts = [z.real for z in res]
        # imaginary_parts = [z.imag for z in res]
        # l2normReal = norm(real_parts,p=2)
        # l2normImag = norm(imaginary_parts, p=2)
        # l2norm = (norm(res, p=2))
        # sepNorm = l2normReal+l2normImag
        # print("normal norm", l2norm)
        # print("Separate norm", sepNorm )
        return residualNorm + penalty
    def Optimize4(self,freq, Xpos, meas, objective): #only optimizing for 2 quantities
        guessEta = 0.048 #global eta, fixed this time
        k0 = [self.G2,self.eta2]#use the guessed values
             #given the algorithm of Nelder-Mead (uses 5% change for vertices of Simplex), using a higher eta is better than close to 0, to give it a big area to look
        # make precision variables
        mp.dps = 20 #desired precision
        result = minimize(objective, k0, args =(freq, meas, Xpos, guessEta),method = 'Nelder-Mead')
        return result
    def Objective4(self, k0, freq, meas, Xpos, guessEta):
        #print(k0)
        fSol = freq
        EtaSol = guessEta
        G2 = (k0[0])
        eta2 = (k0[1])
        self.updateg_(eta2, G2)
        c2 = self.g_ *(1+self.Y)
        
        penalty = 0.0
        penalty = penalty + 1e4 if k0[0] < 0 else penalty +0.0
        penalty = penalty + 1e4 if k0[1] < 0 else penalty +0.0
        penalty = penalty + 1e4 if k0[1] > 3 else penalty +0.0
        fSol = freq 
        pSol = fSol*2*np.pi 
        OmSol = pSol * self.t0 
        Om_Sol = OmSol*(1+1j*EtaSol)**0.5

        ''' I need to calculate the mcs and mgs here!!!!'''
        mc2 = self.g_ *(1+self.Y)
        mg = self.g_

        # kappa = Om_Sol**2
        
        # mkappa = (kappa)
        # mroots = (polyroots([1,-mc2,-mkappa, mkappa*mg],maxsteps=100))
        # mrootsd = matrix([mroots[0],mroots[0], mroots[1], mroots[1], mroots[2],mroots[2]])
        # mrootsd[0] = mrootsd[0]*-1 #according to the paper this needs this change of sign...
        # for i,r in enumerate(mrootsd): 
        #     mrootsd[i] = sqrt((mrootsd[i]))
        coeffs = [1, -mc2, -Om_Sol**2, Om_Sol**2 * mg]
        # print("Coeffs2", coeffs)
        x_roots = np.roots(coeffs)
        k_roots = np.column_stack((np.sqrt(x_roots), -np.sqrt(x_roots))).ravel()
        if (abs(k_roots[0].real) > 10): #todo this 10 is completely arbitrary
            k_roots = k_roots[2:6]
            # print(k_roots)
            A = getSmallAMatrix(k_roots, Xpos)
        else: 
            k_roots[0] = k_roots[0]*1j
            A = getAMatrix(k_roots, Xpos) #matrix from paper
        
        x = lu_solve(A, meas) #not working because of really large numbers (numerically singular matrix)
        res = residual(A, x, meas)
        curveMp = A*x
        curve = np.array(curveMp.tolist(), dtype=complex)
        global iteration
        if (iteration%20 == 0 ):
            print(k0)
            updateGraph(Xpos, meas, curve)
        iteration = iteration +1
        # C = A[0:-40:1,:] #reduced A matrix until the last forty
        # x2 = lu_solve(C, meas[0:-40:1])
        # # res = C * x2 - meas[0:-40:1]
        # res = residual(C, x2, meas[0:-40:1])

        # B = np.array(A.tolist(), dtype=complex) #A matrix but as numpy array
        # column_norms = np.linalg.norm(B,axis=0)
        # Bnorm  = B / np.maximum(column_norms, 1e-10)
        # x, res, _, _ = np.linalg.lstsq(Bnorm, meas, rcond=None)
        # res2 = Bnorm@x - meas

        # B2 = B[0:-40:1, : ] #reducdd numpy array
        # x1, res2 = np.linalg.lstsq(B2, meas[0:-40:1], rcond = None)
        # # res = residual(B2, x1, meas)
        real_parts = [z.real for z in res]
        imaginary_parts = [z.imag for z in res]
        l2normReal = norm(real_parts,p=2)
        l2normImag = norm(imaginary_parts, p=2)
        l2norm = (norm(res, p=2))
        sepNorm = l2normReal+l2normImag
        print("normal norm", l2norm)
        print("Separate norm", sepNorm )
        return l2norm + penalty
    def Optimize5(self,freq, Xpos, meas, objective):
        guessEta = 0.114 #global eta, then try with the exact one
        k0 = [guessEta, self.G2,self.eta2]#use the guessed values
             #given the algorithm of Nelder-Mead (uses 5% change for vertices of Simplex), using a higher eta is better than close to 0, to give it a big area to look
        # make precision variables
        mp.dps = 20 #desired precision
        result = minimize(objective, k0, args =(freq, meas, Xpos),method = 'Nelder-Mead')
        print("result: ", result.x)
        return result
    def Objective5(self, k0, freq, meas, Xpos):
        #print(k0)
        fSol = freq
        EtaSol = (k0[0])
        G2 = (k0[1])
        eta2 = (k0[2])
        self.updateg_(eta2, G2)
        
        penalty = 0.0
        penalty = penalty + 1e4 if k0[0] < 0 else penalty +0.0
        penalty = penalty + 1e4 if k0[1] < 0 else penalty +0.0
        penalty = penalty + 1e4 if k0[2] < 0 else penalty +0.0
        penalty = penalty + 1e4 if k0[2] > 3 else penalty +0.0
        fSol = freq 
        pSol = fSol*2*np.pi 
        OmSol = pSol * self.t0 
        Om_Sol = OmSol*(1+1j*EtaSol)**0.5

        ''' I need to calculate the mcs and mgs here!!!!'''
        mc2 = self.g_ *(1+self.Y)
        mg = self.g_

        coeffs = [1, -mc2, -Om_Sol**2, Om_Sol**2 * mg]
        # print("Coeffs2", coeffs)
        x_roots = np.roots(coeffs)
        k_roots = np.column_stack((np.sqrt(x_roots), -np.sqrt(x_roots))).ravel()
        k_roots = k_roots * 1j #todo delete this just added to test if an M-mj configuration is better
        k_roots[2] = k_roots[2].real - k_roots[2].imag*1j #todo this whole manipulation of
        k_roots[3] = k_roots[3].real - k_roots[3].imag*1j #the wavenumbers did help, but why?? and why was it wrong???
        k_roots[4] = k_roots[4].real - k_roots[4].imag*1j
        k_roots[5] = k_roots[5].real - k_roots[5].imag*1j
        k_roots[0] = k_roots[0].real*1j - k_roots[0].imag
        k_roots[1] = k_roots[1].real*1j - k_roots[1].imag
        # print(k_roots)
        if (abs(k_roots[0].real) > 100): #todo this 10 is completely arbitrary
            k_roots = k_roots[2:6]
            # print(k_roots)
            A = getSmallAMatrix(k_roots, Xpos)
            A2 = getSmallNumpyAMatrix(k_roots, Xpos)
        else: 
            k_roots[0] = k_roots[0]*1j
            A2 = getNumpyAMatrix(k_roots, Xpos)
            print("k0included")
        # column_norms = np.linalg.norm(A2,axis=0)#this normalization is new because from olf measurements it wasnt working without it, maybe because of big x?
        # Anorm  = A2 / np.maximum(column_norms, 1e-10)
        
        # Solve least squares for x
        # x, residualNorm, _, _ = np.linalg.lstsq(Anorm, meas, rcond=None)
        Aoriginal = A2
        x2, residualNorm2, _, _ = np.linalg.lstsq(A2, meas, rcond=None)
        while residualNorm2.size == 0:
            U, S, Vt = np.linalg.svd(A2)
            dependent_col_index = np.argmax(np.abs(Vt[-1]))
            # print("removing dependent column at index:", dependent_col_index)
            A_reduced = np.delete(A2, dependent_col_index, axis=1)
            x2, residualNorm2, rank, s = np.linalg.lstsq(A_reduced, meas, rcond=None)
            A2 = A_reduced
        
        curve = A2@x2
        global iteration
        if (iteration<300 and iteration%20 == 0):
            
            print("k0 at"+str(iteration),k0)
            print("roots:", k_roots)
            # print("resnorm", residualNorm2 + penalty)
            updateGraph(Xpos, meas, curve)
        iteration = iteration +1
        # print("k0",k0)
        objective = residualNorm2 + penalty
        if objective.size != 1:
            a = "stop"
        return objective
def updateGraph(Xpos, meas, curve):
    global iteration
    fig = plt.figure(figsize=(8,6))
    ax0 = fig.add_subplot(221)
    ax0.plot(Xpos, curve.real, label = "curve real")
    ax0.plot(Xpos, meas.real, label = "meas real")
    ax0.set_xlabel("Position (x)")
    ax0.set_ylabel("Real Part")
    plt.title(str(iteration)+ "th iteration")
    plt.legend()
    
    ax0 = fig.add_subplot(222)
    
    ax0.plot(Xpos, curve.imag, label = "curve imag")
    ax0.plot(Xpos, meas.imag, label = "meas imag")
    ax0.set_xlabel("Position (x)")
    ax0.set_ylabel("Imag Part")
    plt.title(str(iteration)+ "th iteration")
    plt.legend()
    
    ax0 = fig.add_subplot(223)
    
    ax0.plot(Xpos, np.angle(curve) , label = "curve phase")
    ax0.plot(Xpos, np.angle(meas), label = "meas phase")
    ax0.set_xlabel("Position (x)")
    ax0.set_ylabel("Phase")
    plt.title(str(iteration)+ "th iteration")
    plt.legend()
    plt.show()
    
    ax0 = fig.add_subplot(224)
    
    ax0.plot(Xpos, np.abs(curve), label = "curve mag")
    ax0.plot(Xpos, np.abs(meas), label = "meas mag")
    ax0.set_xlabel("Position (x)")
    ax0.set_ylabel("Magnitude")
    plt.title(str(iteration)+ "th iteration")
    plt.legend()
    plt.show()
    
    return
def normalize_columns(A):
    col_norms = np.linalg.norm(A, axis=0)  # Compute column norms
    return A / col_norms  # Normalize each column
def normalize_last_rows(A):
    A = A.astype(np.complex_)  # Ensure it's complex if necessary
    for i in [-3,-2, -1]:  # Normalize only the last 4 rows
        row_max = np.max(np.abs(A[i, :]))  # Largest absolute value in row
        if row_max > 0:  # Avoid division by zero
            A[i, :] /= row_max
    return A
def normalize_rows(matrix):
    # Calculate the magnitude (Euclidean norm) of each row
    row_magnitudes = np.linalg.norm(matrix, axis=1)  # norm of each row

    # Normalize each row by dividing by its magnitude
    normalized_matrix = matrix / row_magnitudes[:, np.newaxis]
    
    return normalized_matrix
def Objectiveplz(k0,t0, mc1, mc2, mg):
    penalty = 1e100 if k0[0] < 1 else 0 
    fSol = k0[0]
    EtaSol = k0[1]
    pSol = fSol*2*np.pi 
    OmSol = pSol * t0 
    Om_Sol = OmSol*(1+1j*EtaSol)**0.5
    kappa = Om_Sol**2
    mkappa = kappa #todo adding this line to avoid such a big kappa
    #mkappa = mpc(kappa)
    
    #print(k0)
    # mroots = (polyroots([1,-mc2,-mkappa, mkappa*mg],maxsteps=2000))
    # mrootsd = matrix([mroots[0],mroots[0], mroots[1], mroots[1], mroots[2],mroots[2]])
    # mrootsd[0] = mrootsd[0]*-1 #according to the paper this needs this change of sign...
    # for i,r in enumerate(mrootsd): 
    #     mrootsd[i] = sqrt((mrootsd[i]))
    # print(mrootsd)
    # Create the function once
    

   
    #A_value = mc2  # Example values
    #g_value = mg
    #Omega_value = Om_Sol
    # Now, k_solver can be called multiple times without re-solving the cubic equation
    # mrootsd3 = k_solver(A_value, g_value, Omega_value)
    # print(mrootsd3)
    
    coeffs = [1, -mc2, -Om_Sol**2, Om_Sol**2 * mg]
    x_roots = np.roots(coeffs)
    k_roots = np.column_stack((np.sqrt(x_roots), -np.sqrt(x_roots))).ravel()
    # a = k_roots[0]
    # k_roots[0] = k_roots[2]
    # k_roots[2] = a
    k_roots[0] = k_roots[0]*1j

    # print(k_roots)
    #mrootsd2 = CharEqn(mc2, mg, Om_Sol)
    #print(mrootsd2)
    # mdett = makeMatrix(mrootsd3, mc1, mc2, mkappa)
    mdett2 = makeMatrix(k_roots, mc1, mc2, mkappa)
    # print(mdett)
    # print(mdett2)
    mdett = complex(mdett2)
    obj = abs(mdett) #using the absolute number of the (complex) determinant 
    obj = abs(mdett.real)+abs(mdett.imag)
    return obj + penalty
def Optimize(fSol, EtaSol, t0, mc1, mc2, mg, objective):

    k0 = [fSol,EtaSol]#use 0.2 for eta and the
    result = minimize(objective, k0, args =(t0,mc1,mc2,mg), method = 'Nelder-Mead')
    return result
def CharEqn(mc2, mg, Om_s): #todo, symbolic part could be done outside of the loop of the optimization, so there is already a solved cubic eqn as a function just to be used
    # Define symbols
    A, Omega, g, x = sp.symbols('A Omega g x', complex=True)
    
    # Define the cubic equation
    cubic_eq = x**3 - A*x**2 - Omega**2*x + Omega**2*g
    
    # Solve for x symbolically
    x_roots = sp.solve(cubic_eq, x)
    x_roots[0] = x_roots[0]* -1
    # Compute k values (square roots of x)
    k_roots = []
    for x_i in x_roots:
        sqrt_xi = sp.sqrt(x_i)  # Principal square root
        k_roots.append(sqrt_xi)
        k_roots.append(-sqrt_xi)  # Include negative sqrt
    
    # Substitute values for A, Omega, and g
    A_value = mc2
    g_value = mg
    
    # Create substitutions dictionary
    Om_ = Om_s#
    #Om_ = Om_s[0] + Om_s[1]*1j
    subs_dict = {A: A_value, g: g_value, Omega: Om_}
    
    # Compute numerical values for k
    numerical_k = [k.subs(subs_dict).evalf() for k in k_roots]
    numerical_ks = [complex(ks.evalf()) for ks in numerical_k]
    # nprint(polyroots([1,-1*A_value,-1*Om_*Om_,Om_*Om_*g_value],maxsteps=1000))
    return np.array(numerical_ks)

def generate_k_function():
    # Define symbols
    A, Omega, g, x = sp.symbols('A Omega g x', complex=True)
    
    # Define the cubic equation
    cubic_eq = x**3 - A*x**2 - Omega**2*x + Omega**2*g
    
    # Solve for x symbolically (only once)
    x_roots = sp.solve(cubic_eq, x)
    x_roots[0] = x_roots[0]* -1
    # Compute k values (square roots of x)
    k_roots = []
    for x_i in x_roots:
        sqrt_xi = sp.sqrt(x_i)  # Principal square root
        k_roots.append(sqrt_xi)
        k_roots.append(-sqrt_xi)  # Include negative sqrt

    # Convert k_roots into a lambda function
    def k_function(A_value, g_value, Omega_value):
        subs_dict = {A: A_value, g: g_value, Omega: Omega_value}
        numerical_k = [k.subs(subs_dict).evalf() for k in k_roots]
        return [complex(ks) for ks in numerical_k]
    
    return k_function



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
    c = mpc(c)
    c2 = mpc(c2)
    kappa  = mpc(kappa)
    
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
    
    
    # print(A)
    dett = mp.det(A)
    
    #is preconditioning here a good idea?
    # j=0
    # max_real = max(abs(A[i, j].real) for i in range(A.rows))
    # for i in range(A.rows):
    #     A[i, j] /= max_real
    # j=1
    # max_real = max(abs(A[i, j].real) for i in range(A.rows))
    # for i in range(A.rows):
    #     A[i, j] /= max_real
    # dett2 = mp.det(A)
    
    j = 0 

    # Compute mean absolute real part
    mean_real = sum(abs(A[i, j].real) for i in range(A.rows)) / A.rows
    
    # Normalize by mean real part
    for i in range(A.rows):
        A[i, j] /= mean_real
    
    j = 1  

    # Compute mean absolute real part
    mean_real = sum(abs(A[i, j].real) for i in range(A.rows)) / A.rows
    
    # Normalize by mean real part
    for i in range(A.rows):
        A[i, j] /= mean_real
    dett3 = mp.det(A)
    # dett2 = mp.det(A_precond)
    # print(dett)
    # print(dett3)
    #In the end all the preconditioning attemps are shite. 
    return dett3

def getAMatrix(roots, Xpos):
    k1 = roots[0]
    k2 = roots[2]
    k3 = roots[4]
    n = Xpos.size
    # L = Xpos[-1]-Xpos[0]
    # start = Xpos[0]
    # xis = (Xpos - start)/L #should be a nondimensional Xpos from 0 to 1
    A = matrix(n, 6)
    for i, xi in enumerate(Xpos):
        V = matrix([[sin(k1*xi), cos(k1*xi), exp(k2*xi), exp(-1*k2*xi) ,exp(k3*xi) ,exp(-1*k3*xi) ]])
        A[i,:]=V
    return A
def getSmallAMatrix(roots, Xpos):
    
    k2 = roots[0]
    k3 = roots[2]
    n = Xpos.size
    # L = Xpos[-1]-Xpos[0]
    # start = Xpos[0]
    # xis = (Xpos - start)/L #should be a nondimensional Xpos from 0 to 1
    A = matrix(n, 4)
    for i, xi in enumerate(Xpos):
        V = matrix([[exp(k2*xi), exp(-1*k2*xi) ,exp(k3*xi) ,exp(-1*k3*xi) ]])
        A[i,:]=V
    return A
def getSmallNumpyAMatrix(roots, Xpos):
    k2 = roots[0]
    k3 = roots[2]
    A = np.column_stack([
        np.exp(k2 * Xpos),
        np.exp(-k2 * Xpos),
        np.exp(k3 * Xpos),
        np.exp(-k3 * Xpos),
    ])
    return A
def getNumpyAMatrix(roots, Xpos):
    k1 = roots[0]
    k2 = roots[2]
    k3 = roots[4]
    A = np.column_stack([
        np.sin(k1*Xpos),
        np.cos(k1*Xpos),
        np.exp(k2 * Xpos),
        np.exp(-k2 * Xpos),
        np.exp(k3 * Xpos),
        np.exp(-k3 * Xpos),
    ])
    return A
def getTinyAMatrix(roots, Xpos):
    
    k2 = roots[0]
    n = Xpos.size
    # L = Xpos[-1]-Xpos[0]
    # start = Xpos[0]
    # xis = (Xpos - start)/L #should be a nondimensional Xpos from 0 to 1
    A = matrix(n, 2)
    for i, xi in enumerate(Xpos):
        V = matrix([[exp(k2*xi), exp(-1*k2*xi)]])
        A[i,:]=V
    return A
if __name__ == "__main__":  
    #Example Input
    L = 0.48 #to match the length of a bw beam
    rho1 = 7850
    rho3 = rho1
    E1 = 206e9
    E3 = E1
    
    rho2 = 2600
    eta2 = 0.0 #originally 0.1
    G2 = 0.98e10
    
    b = 0.05 #this is not given in the paper. check if it affects anything. 
    h1 = 0.003 / 2
    h2 = 0.005 /2
    h3 = 0.008 / 2
    
    #Create Sandwich object
    example = SandwichBeam(h1,E1, rho1, h3, E3, rho3, h2, G2, eta2, rho2, b,L)
    xx = example.findModes(1000)
    print("First Mode Freq and Eta:", xx.x)
    '''
    testing the identification procedure
    '''
    # example.L = 0.48 
    # current_file_path = os.path.abspath(__file__)
    # mt_eq_path = os.path.abspath(os.path.join(current_file_path, '..', '..', '..'))

    # file_path = os.path.join(mt_eq_path, 'CodeAsterModels/YACSauto/YACS_BW/autoRes/res0.txt')#get freqiemcies and measurements
     
    # example.IdentifyG_(file_path, "BW")
    '''
    procedure to manually get determinants instead of the above automatic line
    '''
    mp.dps = 200 #desired precision
    k_solver = generate_k_function()# generate the cubic function
    #Solution Data
    fSolCorrect1 =1309.6227567376307661106693 #100000050 #1309 #first mode
    fSolCorrect2 = 6984.5539 #second mode
    fSolCorrect3 = 16849.8283 #third mode
    EtaSol =0.01 #0.0061111813975 #0.00696499999 # 0.006965 
    fSol =800 #1308.9999999 #1504.633 this using complex values now, but correct roots?
    # the identification algorith is very sensitive to the initial guesses :(
    
    #Calculate Frequency Parameter
    pSol = fSol*2*np.pi 
    OmSol = pSol * example.t0 
    Om_Sol = OmSol*(1+1j*EtaSol)**0.5 
    Om_ = [Om_Sol.real, Om_Sol.imag]
    
    # roots = example.CharEqn(Om_)
    c1 = example.g_*example.Y
    c2 = example.g_ *(1+example.Y)
    kappa = Om_Sol**2
    # dett = makeMatrix(roots,c1 , c2, kappa)
    # print("Normal numerical determinant",dett)
    # make precision variables
    
    mc1 = mpc(c1)
    mc2 = mpc(c2)
    mkappa = mpc(kappa)
    mg = mpc(example.g_)
    
    #calculate 3 cubic roots, and make them repeat themselves. The definitoion
    #for the roots is different for the undamped case. See paper "new method bla bla"
    # mroots= (polyroots([1,-c2,-kappa, kappa*example.g_],maxsteps=1000)) #for some reason this formulation works and the one below doesnt, I will use this one then
    #mroots2 = (polyroots([1,-mc2,-mkappa, mkappa*mg],maxsteps=1000))
    # mrootsd = matrix([mroots[0],mroots[0], mroots[1], mroots[1], mroots[2],mroots[2]])
    # mrootsd[0] = mrootsd[0]*-1 #according to the paper this needs this change of sign...
    # for i,r in enumerate(mrootsd): 
    #     mrootsd[i] = sqrt((mrootsd[i])) #get the 6th root by sqrt of cubic root
    
    # mdett = makeMatrix(mrootsd, mc1, mc2, mkappa)
    print("Guess frequency:",fSol)
    print("Guess eta:", EtaSol)
    # print("Determinant:", mdett)
    
    #use optimization technique. not member function yet.
    x = Optimize(fSol, EtaSol, example.t0, c1, c2, example.g_, Objectiveplz)
    optFreq = x.x[0]
    optEta = x.x[1] 
    print("Optimal:", x.x)
    '''
    older code that I dont use anymore.
    '''
    # det = example.DetFuncTrig(Om_, roots)
    # det2 = example.getDet(Om_)
    #root_lambda = example.find_lambda_root(Om_)
    # result = minimize(example.getDet, Om_,method = 'Nelder-Mead', tol = 1e-30)
    # Om_test = result.x

