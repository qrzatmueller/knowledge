'''
Simple file to compare one of the approximations comonly used for low loss 
factor.  which says that 1/(1+ja) =approx 1-ja
'''

import numpy as np
import matplotlib.pyplot as plt

# Define the range of zk_ values (imaginary component)
eta = np.linspace(0, 1.5, 100)

# Compute the actual and approximate expressions
actual = 1 / (1 + 1j*eta)
approximation = 1 - 1j*eta 
etaBecker = eta * (1-eta*eta/16)/ (1-3*eta*eta/8+eta*eta*eta*eta/256)
approximationB = 1- 1j*etaBecker
# Compute the real and imaginary differences
diff_real = np.real(actual - approximation)
diff_imag = np.imag(actual - approximation)

# Create scatter plots for the differences
plt.figure(figsize=(12, 6))

# Scatter plot for the real part
plt.subplot(1, 2, 1)
plt.scatter(eta, np.real(actual), c="blue", label="Real Part Exact")
plt.scatter(eta, np.real(approximation), c="red", label="Real Part Approximation")
plt.axhline(0, color="black", linestyle="--", linewidth=0.8)
plt.title("Difference in Real Part")
plt.xlabel("Eta")
plt.ylabel("Difference")
plt.legend()

# Scatter plot for the imaginary part
plt.subplot(1, 2, 2)

plt.scatter(eta, np.imag(actual), c="blue", label = "Imaginary Exact")
plt.scatter(eta, np.imag(approximation), c="red", label = "Imag Approximation")
plt.axhline(0, color="black", linestyle="--", linewidth=0.8)
plt.title("Difference in Imaginary Part")
plt.xlabel("Eta")
plt.ylabel("Difference")
plt.legend()

plt.tight_layout()
plt.show()