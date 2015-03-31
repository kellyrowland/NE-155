import numpy as np
import matplotlib.pyplot as plt
from math import *

x_max = 4.0
x_min = -4.0
D = 1.0
Sigma_a = 0.2
h = 0.1
L = sqrt(D/Sigma_a)
n = int((x_max - x_min)/h)
x_range = np.linspace(x_min,x_max,num=n+1)
case = str(raw_input("Please enter case (b or c): "))

A = np.matrix(np.zeros((n+1,n+1)))
a = np.repeat(0.0,n+1)
b = np.repeat(0.0,n+1)
c = np.repeat(0.0,n+1)
for i in range(0,n+1):
        A[i,i] = 2 + (h/L)**2
        b[i] = A[i,i]
        if i <= n-1:
                A[i+1,i] = -1.0; A[i,i+1] = -1.0
                a[i] = -A[i+1,i]; c[i] = -A[i,i+1]
phi = np.repeat(0.0,n+1)
S = np.repeat(0.0,n+1)
S_matrix = np.repeat(0.0,n+1)
if case == "b":
	for i in range(0,n+1): S[i] = 8.0
elif case == "c":
	for i in range(0,n+1): S[i] = cos(x_range[i])
for i in range(0,n+1): S_matrix[i] = (h**2)*(S[i]/D)

u = np.repeat(0.0,n+1)
v = np.repeat(0.0,n+1)

u[0] = b[0]; v[0] = S_matrix[0]
for i in range(1,n+1):
        u[i] = b[i] - (a[i]*c[i-1])/u[i-1]
        v[i] = S_matrix[i] + (a[i]*v[i-1])/u[i-1]

phi[n] = v[n]/u[n]
for i in range(n-1,0,-1):
        phi[i] = (v[i] + c[i]*phi[i+1])/u[i]

numer = np.repeat(0.0, n+1)
denom = np.repeat(0.0, n+1)
for i in range(0,n+1):
	if case == "b":
        	numer[i] = -(S[i]/Sigma_a)*(exp(sqrt(Sigma_a/D)*x_range[i]) + exp(-sqrt(Sigma_a/D)*x_range[i]))
        	denom[i] = exp(sqrt(Sigma_a/D)*x_max) + exp(-sqrt(Sigma_a/D)*x_max)
	elif case == "c":
        	numer[i] = -(cos(x_max))*(exp(sqrt(Sigma_a/D)*x_range[i]) + exp(-sqrt(Sigma_a/D)*x_range[i]))
        	denom[i] = (D+Sigma_a)*(exp(sqrt(Sigma_a/D)*x_max) + exp(-sqrt(Sigma_a/D)*x_max))
if case == "b": phi_x = (S/Sigma_a)+(numer/denom)
elif case == "c": phi_x = (S/(D+Sigma_a)) + (numer/denom)
plt.plot(x_range, phi_x, 'x')
plt.plot(x_range, phi, '.')
plt.legend(("Analytical","Numerical"), loc="upper right")
plt.xlabel('x [cm]'); plt.ylabel('Flux [n/cc/s]')
plt.title('Flux Solutions')
plt.savefig('hw6_2.pdf')
