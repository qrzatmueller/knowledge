'''
in this file I am postprocessing data files created between the yacs BW2 and BW_aster_yacs.py
the objective was to show that for high wavelengths, the method starts failing, 
here the loss factor and frequency are converted to files with loss factor relative error
and a ratio of lamda and the measurement aperture. 
for this a lamda is calculated 
'''

import numpy as np
import matplotlib.pyplot as plt
import BeamWaves as BW

# ----- Given parameters -----
eta_ideal = 0.5   # Ideal loss factor (dimensionless)
E = 210e09*(1+eta_ideal*1j)        # Young's modulus [Pa]
rho = 7850        # Density [kg/m^3]
h = 0.010         # Thickness [m]
b = 0.05          # Width [m]
L_meas = 0.40     # Measurement length [m] 

# ----- Read the .dat file -----
# Format: freq [Hz]    loss_factor
data = np.loadtxt('Inv/02BWLowFreq/Data/res0BWIDn.dat')
frequencies = data[:, 0]       # Hz
loss_factors = data[:, 1]      # dimensionless

# ----- Beam geometry -----
I = (b * h**3) / 12            # Second moment of area [m^4]
A = b * h                      # Cross-sectional area [m^2]

# ----- Compute wavelength -----
omega = 2 * np.pi * frequencies
wavelengths = 2 * np.pi * (( (E * I) / (rho * A * omega**2) )**0.25).real
beam = BW.EBBeam(b = b, E = E.real, h = h, rho=rho )
#frequencies = np.linspace(0,1000,101)
#cb = beam.BendingWaveSpeed(frequencies, plot = "h=1mm, b=30mm, Steel")  
#k = beam.WaveNumber(frequencies, plot = "h=1mm, b=50mm, Steel")
#lam = beam.WaveLength(frequencies, plot = "h=1mm, b=50mm, Steel")

# ----- Non-dimensional length -----
non_dim_lengths = L_meas / wavelengths 

# ----- Relative loss factor error -----
loss_factor_error = np.abs((loss_factors - eta_ideal) / eta_ideal)

# ----- Optional plotting here
plt.figure(figsize = (8,5))
plt.plot(non_dim_lengths, loss_factor_error, 'o-', label='Relative Error')
plt.axhline(0, color='gray', linestyle='--', linewidth=1)
plt.xlabel('Non-dimensional Length (Lm / λ)')
plt.ylabel('Relative Loss Factor Error')
plt.title('Loss Factor Error vs. Non-dimensional Length')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


# ----- Save output -----
output = np.column_stack((non_dim_lengths, loss_factor_error))
np.savetxt('Inv/02BWLowFreq/Data/errVsLamdratio.dat', output, header='NonDimLength  RelLossFactorError', fmt='%.6e')
