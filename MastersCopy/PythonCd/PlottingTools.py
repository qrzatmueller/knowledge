import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

'''
File that shows a way to print the vibration in the complex plane along a 3rd coordinate
in this case along space.
'''

# Define x values
x = np.linspace(0, 10, 100)  # Positions along x

# Define a complex function
complex_values = np.exp(1j * x)  # Example function: e^(ix)

# Extract real and imaginary parts
real_part = complex_values.real
imag_part = complex_values.imag

# Create 3D plot
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')

# Plot the complex function evolution
ax.plot(x, real_part, imag_part, label="Complex function")

# Labels
ax.set_xlabel("Position (x)")
ax.set_ylabel("Real Part")
ax.set_zlabel("Imaginary Part")
ax.set_title("Complex Function Evolution")

ax.legend()
plt.show()
