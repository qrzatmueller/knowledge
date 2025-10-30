# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 18:01:46 2025

@author: chatgpt
chatgpt example of how to annimate a plot
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

# Parameters
n_points = 50                      # Number of points along the beam
x = np.linspace(0, 1, n_points)    # Beam position (normalized from 0 to 1)
omega = 2 * np.pi                  # Frequency of oscillation
t_max = 2                           # Animation time duration
fps = 30                            # Frames per second
num_frames = fps * t_max            # Total number of frames

# Define the initial complex displacement (Example: Some wave pattern)
complex_displacement = np.exp(-5 * (x - 0.5)**2) * (1 + 1j)  # Example Gaussian shape

# Create figure and 3D axis
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Initial plot (empty)
line, = ax.plot([], [], [], 'o-', markersize=5, label=r"$\psi(x,t)$")

# Labels
ax.set_xlabel("Position $x$")
ax.set_ylabel("Real Part")
ax.set_zlabel("Imaginary Part")
ax.set_title("Complex Displacement Evolution")
ax.legend()

# Set axes limits
ax.set_xlim(0, 1)
ax.set_ylim(-1.5, 1.5)
ax.set_zlim(-1.5, 1.5)

# Animation update function
def update(frame):
    t = frame / fps  # Convert frame to time
    phase_factor = np.exp(1j * omega * t)  # e^(iωt) term
    new_displacement = complex_displacement * phase_factor  # Apply time evolution
    
    # Extract real and imaginary parts
    y_real = new_displacement.real
    z_imag = new_displacement.imag
    
    # Update the 3D line
    line.set_data(x, y_real)
    line.set_3d_properties(z_imag)
    return line,

# Create animation
ani = animation.FuncAnimation(fig, update, frames=num_frames, interval=1000/fps, blit=False)

plt.show()
