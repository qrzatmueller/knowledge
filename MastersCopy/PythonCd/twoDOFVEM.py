# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 10:59:56 2024

@author: prqrz
Example 8.2 from structural dynamics script, erweitert um die VEM mit 
frequenz abhängige eigenschaften zu benutzen
"""
import numpy as np
#bestimmen eigenschaften k1, k2, m1, m2 für frequenz unabhängie Fall 

k1 = 200 #[N/m]
k2 = 125
m1 = 80 #[kg] 
m2 = 8 

#eigenwert lösung durch determinant 
a= m1*m2
b= -(k1*m1 + k1*m2 + k2*m1)
c=  k1*k2
l1 = (-b + np.sqrt(b*b-4*a*c)) / (2*a)
l2 = (-b - np.sqrt(b*b-4*a*c)) / (2*a)

#eigenfrequenzen
om1 = np.sqrt(l1)
om2 = np.sqrt(l2)

#frequency dependent properties
f_array = [1,10,20,30,40]
k1_array = [180, 200 ,220, 200, 150]
k2_array = [100, 125, 145, 125, 80]

f__unt = np.linspace(0, 50,101)

Kc = np.array([[4+2j, 1-1j],[1+1j, 3-3j]])
M = np.array([[2,0],[0,1]])
M_inv = np.linalg.inv(M)
A = np.dot(M_inv, Kc)

eigenvalues, eigenvectors = np.linalg.eig(A)

print("Eigenvalues", eigenvalues)
print("eigenvectors", eigenvectors) 
