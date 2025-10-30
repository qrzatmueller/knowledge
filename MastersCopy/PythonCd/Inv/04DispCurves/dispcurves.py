'''
This file plots the adimensional dispercion curve from the graff book
it normalizes the wavenumber by multiplying with a characteristic length
(radius or thickness). The one in the book is for circular bars, 
here I also added for rectangular bars. kappa depends also on this geometry, so change accordingly
below there is also a dispercion graph, with units and dependent on frequency
if that is more helpfull for a specific case.
'''

import numpy as np
import matplotlib.pyplot as plt

# Material and geometry
E = 1e9           # Young's modulus [Pa]
rho = 7850           # Density [kg/m^3]
nu = 0.3             # Poisson's ratio
G = E / (2 * (1 + nu))   # Shear modulus [Pa]
kappa = 5/6          # Shear correction factor 9/8 for round?? checkkk

b = 0.05             # Width [m]
h = 0.001            # Thickness [m]
A = b * h            # Cross-sectional area [m^2]
I = (b * h**3) / 12  # Second moment of area [m^4]
r = 0.001

A_c = np.pi * r**2
I_c = (np.pi * (2*r)**4) /64
kappa_c=10/9 #for circular bar 
#using circular bar
# A=A_c
# I = I_c
# kappa=kappa_c

f = np.linspace(0,20000, 20000)
f = np.linspace(0,1000000, 20000)
om = f * 2 * np.pi

a = E*I/rho/A 
b = -I/A*(1+E/G/kappa)*om**2
c = -1*(om**2) + (rho*I/G/A/kappa)*om**4
ksquared1 =(-b - np.sqrt(b**2 - 4 * a * c)) / (2 * a)
ksquared2 =(-b + np.sqrt(b**2 - 4 * a * c)) / (2 * a)
k1 = ksquared1**0.5
k2 = ksquared2**0.5



k_eb = ((rho*A/E/I)*om**2)**0.25


c = om/k2
c_eb = om/k_eb

#normalization for circular bar
k2norm = r*k2 / 2 /np.pi 
k_eb_norm = k_eb *r /2/np.pi
c0 = np.sqrt(E/rho)
c_norm = c/c0
c_eb_norm = c_eb/c0

#normalization for rectangular 
k2norm = h * k2 /2/np.pi #dividing thickness over wavelength
k_eb_norm = k_eb *h /2/np.pi


# --- Plot ---
plt.figure(figsize=(8, 5))
plt.plot(k2norm, c_norm, '--', label='Timoshenko (bending mode)', lw=2)
plt.plot(k_eb_norm, c_eb_norm, label='Euler-Bernoulli', lw=2)
plt.xlabel("Norm. Wavenumber")
plt.ylabel("Norm. Phase Velocity")
plt.ylim(0.0, 1.0)
plt.title("Normalized Phase Velocity vs Normalized Wave Number (Graff)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# plt.figure(figsize=(8, 5))
# plt.plot(k2, c, '--', label='Timoshenko (bending mode)', lw=2)
# plt.plot(k_eb, c_eb, label='Euler-Bernoulli', lw=2)
# plt.xlabel("Wavenumber [1/m]")
# plt.ylabel("Phase Velocity [m/s]")
# plt.title("Phase Velocity vs Wave Number")
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.show()

plt.figure(figsize=(8, 5))
plt.plot(f, c, '--', label='Timoshenko (bending mode)', lw=2)
plt.plot(f, c_eb, label='Euler-Bernoulli', lw=2)
plt.xlabel("Frequency [Hz]")
plt.ylabel("Phase Velocity [m/s]")
plt.title("c vs f")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# The same again but for other properties
E = 210e6           # Young's modulus [Pa]
rho = 7850           # Density [kg/m^3]
nu = 0.3             # Poisson's ratio
G = E / (2 * (1 + nu))   # Shear modulus [Pa]
kappa = 5/6          # Shear correction factor

b = 0.05             # Width [m]
h = 0.01            # Thickness [m]
A = b * h            # Cross-sectional area [m^2]
I = (b * h**3) / 12  # Second moment of area [m^4]
r = 0.001

A_c = np.pi * r**2
I_c = (np.pi * (2*r)**4) /64
kappa_c=10/9 #for circular bar 
#using circular bar
# A=A_c
# I = I_c
# kappa=kappa_c

f = np.linspace(0,1000, 20000)
om = f * 2 * np.pi

a = E*I/rho/A 
b = -I/A*(1+E/G/kappa)*om**2
c = -1*(om**2) + (rho*I/G/A/kappa)*om**4
ksquared1 =(-b - np.sqrt(b**2 - 4 * a * c)) / (2 * a)
ksquared2 =(-b + np.sqrt(b**2 - 4 * a * c)) / (2 * a)
k1 = ksquared1**0.5
k2 = ksquared2**0.5



k_eb = ((rho*A/E/I)*om**2)**0.25


c = om/k2
c_eb = om/k_eb


plt.figure(figsize=(8, 5))
plt.plot(k2, c, '--', label='Timoshenko (bending mode)', lw=2)
plt.plot(k_eb, c_eb, label='Euler-Bernoulli', lw=2)
plt.xlabel("Frequency [Hz]")
plt.ylabel("Phase Velocity [m/s]")
plt.title("Plot from Graff2")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()