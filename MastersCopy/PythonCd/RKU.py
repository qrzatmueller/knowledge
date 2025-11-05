# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 12:27:35 2025
This file calculates the homogenized properties of a beam, starting from the 
properties of the layers of an unconstrained damped laminates. Basically 2
layers, one stiff (steel) and one viscoelastic. Afterwards an optimization 
algorithm identifies the loss factor of the layer based on the results

From the identification process, I can say that identifying only one parameter
like the loss factor of the viscoelastic layer, works, but trying to identify
both the youngs modulus and the loss factor does not work. There are more than
one combination that satisfy the equation. Because we have 1 equation and 2 
variables. But we do have the overall stiffness, that can probably make a second
equation that limits the stiffness. How can forward calculations give me also 
an overall youngs modulus?



@author: chatgpt for rku_method function, me with the RKU publication for the 
rku_plates function. 
"""
import numpy as np
from scipy.optimize import minimize
from sympy import symbols, solve, sqrt
def rku_method(E1, h1, E3, h3, G2, eta2, h2):
    """
    Computes the effective bending stiffness, complex modulus, and loss factor of a sandwich beam using RKU.
    
    Parameters:
    E1  : float  - Young's modulus of top face layer (Pa)
    h1  : float  - Thickness of top face layer (m)
    E3  : float  - Young's modulus of bottom face layer (Pa)
    h3  : float  - Thickness of bottom face layer (m)
    G2  : float  - Shear modulus of viscoelastic core (Pa)
    eta2: float  - Loss factor of viscoelastic core (dimensionless)
    h2  : float  - Thickness of viscoelastic core (m)
    
    Returns:
    D_eff    : float - Effective bending stiffness (Pa m^3)
    E_eff    : float - Effective complex Young's modulus (Pa)
    eta_eff  : float - Effective loss factor (dimensionless)
    """
    
    # Total thickness
    h = h1 + h2 + h3
    
    # Shear modulus with loss factor
    G2_complex = G2 * (1 + 1j * eta2)
    
    # Contribution from face layers
    D1 = (E1 * h1**3) / 12 + (E1 * h1) * (h/2)**2
    D3 = (E3 * h3**3) / 12 + (E3 * h3) * (h/2)**2
    
    # Contribution from viscoelastic core
    D2 = (G2_complex * h2) / 2
    
    # Effective bending stiffness
    D_eff = D1 + D3 + D2
    
    # Extract real and imaginary parts
    D_real = D_eff.real
    D_imag = D_eff.imag
    
    # Effective Young's modulus (approximate, assuming thin core)
    E_eff = D_real / (h**3 / 12)
    
    # Effective loss factor
    eta_eff = D_imag / D_real
    
    return D_eff, E_eff, eta_eff
def fromMIT(E1,h1,E2,h2,eta2):
    '''
    Returns
    -------
    eta : calculated loss factor
        DESCRIPTION.
        This formula I took from an MIT lecture github jupyter notebook.
        https://github.com/acoular/structure-borne-sound-lecture/blob/master/09_Damping.ipynb
        THe printed equation in the jupyter notebook is somehow faulty, or the
        assumption do not work with our case. Do not use this function or that
        formula. The later symbolic evaluation does work...
        The assumptions were the limit where E2 tends to 0, which is weird 
        becuase E2 is still part of the formula. The jupyter notebook one is 
        yielding the same result almoss as that from RKU_plates below :( 
        I saved the notebook under the jupyter folder anyway...
    '''
    O1 = E2*eta2*h2
    O2 = 3*E1**2*h1**4
    O3 = 6*E1**2*h1**3*h2 
    O4 = 4*E1**2*h1**2*h2**2
    O5 = 2*E1*E2*h1*h2**3
    O6 = E2**2*h2**4
    U1 = E1*h1
    U2 = E1**2*h1**4
    U3 = 2*E1*E2*h1**3*h2
    U4 = 4*E2**2*h1**2*h2**2 
    U5 = 6*E2**2*h1*h2**3 
    U6 = 3*E2**2*h2**4
    eta = O1 * (O2+O3+O4+O5+O6)/(U1*(U2+U3+U4+U5+U6))
    eta_ratio = eta/eta2
    return eta, eta_ratio 


def rku_plates(H1, H2, E1, E2, b, width=1):
    """
    Follows the special case of "Plates with Viscoelastic Coatings" given in 
    page 601 of Ungars Sctructural Damping.
    From the equations one can observe that the width is not important if both
    layers are of the same width and rectangular. 
    Parameters:
    H1    : float  - Thickness of carrying plate (m)
    E1    : float  - Young's modulus of carrying plate (Pa)
    H2    : float  - Thickness of viscoelastic layer (m)
    E2    : float  - Young's modulus of viscoelastic layer (Pa)
    b     : float  - beta, loss factor of viscoelastic layer (m)
    width : float  - width of the plates (m)
    
    Returns:
    eta      : float - Loss factor of the composite (Pa m^3)
    etaratio : float - ratio of Eta/Beta(Pa)
    """
    H12 = 0.5*(H1+H2)
    r1 = H1/np.sqrt(12)
    r2 = H2/np.sqrt(12)
    A1 = H1*width
    A2 = H2*width
    K1 = E1*A1
    K2 = E2*A2
    k = K2/K1
    a = (1+k)**2 + (b*k)**2
    etaratio = 1/(1+((k**2 * (1+b**2)+(r1/H12)**2 * a)/(k*(1+(r2/H12)**2 * a))))
    eta = etaratio * b
    return eta, etaratio

def rku(E1, H1, eta1, E2, H2, eta2, width):
    I1 = width * H1**3 /12
    E2_ = E2*(1+eta2*1j)
    E1_ = E1*(1 + eta1*1j)
    e_ = E2_/E1_
    h = H1/H1
    StiffRatio = 1+e_*h**3 + 3*(1+h)**3*(e_*h/(1+e_*h))
    E_I = StiffRatio * E1_*I1
    eta =(E_I/E_I.real)
    return eta.imag
def objective(guessb,H1,H2,E1,E2,width,MeasEta):
    calcEta, a = rku_plates(H1,H2,E1,E2,guessb,1)
    return np.abs(MeasEta-calcEta)

def objective2(guessParam,H1,H2,E1,width,MeasEta):
    guessb = guessParam[0]
    guessE2 = guessParam[1]
    calcEta, a = rku_plates(H1,H2,E1,guessE2,guessb,1)
    return np.abs(MeasEta-calcEta)
def RAO(L,h1, h2, h3, E1, G2, E3, rho1, rho2, rho3, eta2, b):
    '''
    from the paper by RAO, which takes the 6th order PDE by DiTaranto. 
    use the default parameters to do the 7.1 example
    It is only programed for teh first mode.3 modes could be available if I 
    transcribe the formulas... all modes are available if I solve the equation.
    I will test my CodeAster model to compare with this results
    Parameters
    ----------
    L : TYPE float
        DESCRIPTION. length
    h1 : TYPE float
        DESCRIPTION. half thickness first layer
    h2 : TYPE float
        DESCRIPTION. half thickness 
    h3 : TYPE float
        DESCRIPTION. half thickness bottom layer
    E1 : TYPE float
        DESCRIPTION. Youngs modulus top layer
    G2 : TYPE float
        DESCRIPTION. Shear modulus of mid-layer 
    E3 : TYPE float
        DESCRIPTION. youngs modulus bottom layer
    rho1 : TYPE float
        DESCRIPTION. density of top layer
    rho2 : TYPE float
        DESCRIPTION. density of middle layer (viscous)
    rho3 : TYPE float
        DESCRIPTION. density of bottom (structural) layer
    eta2 : TYPE float
        DESCRIPTION. loss factor of middle layer
    b : TYPE float
        DESCRIPTION. width

    Returns
    -------
    Y : TYPE float
        DESCRIPTION. geometrical parameter to use in diagrams from paper
    g : TYPE float
        DESCRIPTION. shear parameter, eqn 7
    eta : TYPE float
        DESCRIPTION. loss factor of all the sandwich
    f : TYPE float
        DESCRIPTION. first eigenmode frequency of the beam

    '''
    
    I1 = b*(h1*2)**3 / 12
    I3 = b*(h3*2)**3 / 12
    A1 = b*h1*2
    A2 = b*h2*2
    A3 = b*h3*2
    
    D = E1*I1+E3*I3
    m = rho1*A1 + rho2*A2 + rho3*A3 # mass per unit length of beam
    t0 = np.sqrt(m*L**4/D)
    G2_ = G2*(1+1j*eta2)
    g_ = (G2_*A2*L**2 / 4/h2**2) * ((E1*A1+E3*A3)/(E1*A1*E3*A3))
    g = g_.real
    c = 2*h2 + h1+ h3 
    Y = c**2/D * ((E1*A1*E3*A3)/(E1*A1+E3*A3))
    
    ''' the following 3 lines are not used, but optional diagram way'''
    Om0 = 3.516 # choose from table, first mode
    Om_ = 4.5 # find in diagram, with Y = 21.421 and g= 43.6
    eta_ = 0.09 # or 0.08, not easy to see from diagram...
    
    
    Om = Om0*Om_
    eta_ = np.zeros(3)
    Om = np.zeros_like(eta_)
    for i in range(3):
        eta_[i], Om[i] = get_eta_Om_(g,Y,i) #asking for second mode. 
    eta = eta_ * eta2
    f = Om / (2*np.pi*t0)
    return Y, g, eta, f

def get_eta_Om_(g, Y, mode):
    x = np.log10(g)
    y = np.log10(Y)
    if(mode == 0):
        
        eta_ = ( 0.142675 * pow(1.76507,x) * pow(0.381063,pow(x,2)) * 
                pow(0.958039,pow(x,3)) * pow(1.06549,pow(x,4)) * pow(5.78986,y) * 
                pow(0.539055,x*y) * pow(1.08780,pow(x,2)*y) * pow(1.08745,pow(x,3)*y) *
                pow(0.976757,pow(x,4)*y) * pow(0.586965, pow(y,2)) * 
                pow(1.02241 ,x*pow(y,2)) * pow(1.15050, pow(x,2)*pow(y,2)) * 
                pow(0.970210,pow(x,3)*pow(y,2)) * pow(0.983707, pow(x,4)*pow(y,2)))
        
        Om = ( 3.97720 + 0.593861*x + 0.151641*pow(x,2) - 0.061331*pow(x,3) - 
               0.022252*pow(x,4) + y*(0.154016 + 0.553244*x + 0.144153*pow(x,2) - 
                                      0.075502*pow(x,3) - 0.014758*pow(x,4)) +
               pow(y,2)*(2.42584 + 2.14246*x + 0.294401*pow(x,2) - 0.194389*pow(x,3)
                         - 0.048345*pow(x,4)))
    if(mode == 1):
        eta_= ( 0.045809 * pow(5.23955,x) * pow(0.660211,pow(x,2)) * 
               pow(0.897479, pow(x,3)) * pow(0.992201, pow(x,4)) *
               pow(9.96967, y) * pow(0.708600,x*y) * pow(0.737315, pow(x,2)*y) * 
                 pow(1.01858,pow(x,3)*y) * pow(1.04428, pow(x,4)*y) * 
                pow(0.632614,pow(y,2)) * pow(0.782958, x*pow(y,2)) * pow(1.16677, pow(x,2)*pow(y,2))
                * pow(1.03702, pow(x,3)*pow(y,2)) * pow(0.976208, pow(x,4)*pow(y,2)))
      
        Om = (22.6437 + 1.53915*x + 1.16909*x*x + 0.088024*pow(x,3) - 0.096522*pow(x,4)
              + y*(-1.37115 + 0.330011*x + 1.61566*x*x+ 0.233375*x*x*x
                   - 0.169197*pow(x,4)) +y*y*(6.24946 + 8.70502*x + 3.57741*x*x 
                                              - 0.22708*x*x*x - 0.306308*pow(x,4)) )
    if(mode == 2):
       eta_= ( 0.018430 * pow(8.13049,x) * pow(0.803410,pow(x,2)) * 
              pow(0.890123, pow(x,3)) * pow(0.976944, pow(x,4)) *
              pow(10.9134, y) * pow(0.90008,x*y) * pow(0.789523, pow(x,2)*y) * 
                pow(0.964562,pow(x,3)*y) * pow(1.02404, pow(x,4)*y) * 
               pow(0.763190,pow(y,2)) * pow(0.730160, x*pow(y,2)) * pow(0.974320, pow(x,2)*pow(y,2))
               * pow(1.05259, pow(x,3)*pow(y,2)) * pow(1.00884, pow(x,4)*pow(y,2)))
     
       Om = (62.2273 + 1.82075*x + 1.96034*x*x + 0.59842*pow(x,3) - 0.000547*pow(x,4)
             + y*(-3.00169 - 2.56321*x + 1.23125*x*x+ 1.22688*x*x*x
                  - 0.158933*pow(x,4)) +y*y*(8.08929 + 16.5705*x + 10.5650*x*x 
                                             + 0.634568*x*x*x - 0.768191*pow(x,4)) )                                           
                                              
                                              
    return eta_ , Om
# Example usage
if __name__ == "__main__":
    # Material properties
    E1 = 2.1e11  # Pa (steel)
    h1 = 0.001   # m
    E3 = 2.1e11  # Pa (steel)
    h3 = 0.001   # m
    G2 = 1e6     # Pa (viscoelastic core)
    eta2 = 0.5   # Loss factor of core
    h2 = 0.002   # m
    E2 = 210000000
    Eratio = E2/E1
    relH = h2/h1
    eta, etaratio = rku_plates(h1, h2, E1, E2, b=eta2, width = 0.05)
    eta_MIT, eta_ratio_MIT = fromMIT(E1, h1, E2, h2, eta2)

    etaother = rku(E1, h1, 0, E2, h2, eta2, 0.05)
    
    L = 0.1
    rho1 = 7850
    rho3 = rho1
    E1 = 206e9
    E3 = E1
    
    rho2 = 2600
    eta2 = 0.1
    G2 = 0.98e10
    
    b = 0.05 #this is not given in the paper. check if it affects anything. 
    h1 = 0.003 / 2
    h2 = 0.005 /2
    h3 = 0.008 / 2
    Y2, g, etaRao, f = RAO(L, h1, h2, h3, E1, G2, E3, rho1, rho2, rho3, eta2, b)
    
    L = 0.48
    rho1 = 7850
    rho3 = rho1
    E1 = 206e9
    E3 = E1
    
    rho2 = 2600
    eta2 = 0.1
    G2 = 0.98e10
    
    b = 0.05 #this is not given in the paper. check if it affects anything. 
    h1 = 0.003 / 2
    h2 = 0.005 /2
    h3 = 0.008 / 2
    Y2, g, etaRaoLong, fLong = RAO(L, h1, h2, h3, E1, G2, E3, rho1, rho2, rho3, eta2, b)
    
    
    L = 0.1
    rho1 = 7771
    rho3 = 2700
    E1 = 200.56e9
    E3 = 69e9
    
    rho2 = 2870
    eta2 = 2.52
    G2 = 1504699.80207024
    
    b = 0.05 #this is not given in the paper. check if it affects anything. 
    h1 = 0.001 / 2
    h2 = 0.00152 /2
    h3 = 0.00015 / 2
    Y2B, gB, etaRaoB, fB = RAO(L, h1, h2, h3, E1, G2, E3, rho1, rho2, rho3, eta2, b)
    print(etaRaoB)
    print(fB)
    #D_eff, E_eff, eta_eff = rku_method(E1, h1, 0, 0, E2/3, eta2, h2) #this one is not working for a 2 layered material, need testing for 3 layers...
    
    # starterEta = 0.3 #given the algorithm of Nelder-Mead (uses 5% change for vertices of Simplex), using a higher eta is better than close to 0, to give it a big area to look
    
    
    # result = minimize(objective, starterEta, args =(h1,h2,E1,E2,1,eta),method = 'Nelder-Mead', tol = 1e-18)
    # eta2ident = result.x[0]
    # print(f"Identified viscoelastic loss factor was (eta2ident): {eta2ident:.3e} [-]")
    # print(f"Predefined loss factor was (eta2): {eta2:.3e} [-]")
    # # print(f"Effective Loss Factor (eta_eff): {eta_eff:.3f}")
    # starterParam = [0.8,210000000]
    # result2 = minimize(objective2, starterParam, args =(h1,h2,E1,1,eta),method = 'Nelder-Mead')
    
    # etaB, f = rku_plates(h1,h2,E1,result2.x[1],result2.x[0],1)
    
    
    
