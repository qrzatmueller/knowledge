# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 14:57:10 2025

@author: prqrz
"""

import numpy as np
import numpy as np
from mpmath import mp, matrix, mpf, mpc

def cubic_roots(a, b, c, d):
    """Finds the roots of the cubic polynomial ax^3 + bx^2 + cx + d = 0"""
    if a == 0:
        raise ValueError("The coefficient 'a' must be nonzero for a cubic equation.")
    
    # Compute roots
    coefficients = [a, b, c, d]
    roots = np.roots(coefficients)
    
    return roots

def check_roots(a, b, c, d, roots):
    """Substitutes each root into the polynomial and checks if it evaluates to zero"""
    for root in roots:
        value = a * root**3 + b * root**2 + c * root + d
        print(f"Checking root: {root}, P(root) = {value}")
        if not np.isclose(value, 0, atol=1e-6):
            print("Warning: Root check failed!")

# Example usage
a, b, c, d = 1.264, 6.44, 11, 6.17  # Example coefficients for x^3 - 6x^2 + 11x - 6 = 0
roots = cubic_roots(a, b, c, d)
print("Roots:", roots)
check_roots(a, b, c, d, roots)


# Set high precision
mp.dps = 100  # 100 decimal places

# Define an ill-conditioned matrix with very large and very small values
numpy_matrix = np.array([[1e20, 1], [1, 1e-20]])

# Convert to mpmath matrix with arbitrary precision
mpmath_matrix = matrix([[mpf(str(m)) for m in row] for row in numpy_matrix])

# Compute the determinant
det_A = mpmath_matrix.det()

# Print the result
print(det_A)
