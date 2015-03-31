#jacobi
import numpy as np
import matplotlib.pyplot as plt
from math import *

def first_sum(i,A,x):
	sum = 0.0
	for j in range(0,i): sum += A[i,j]*x[j]
	return sum

def second_sum(i,A,x):
	sum = 0.0
	for j in range(i+1,len(x)): sum += A[i,j]*x[j]
	return sum

def Jacobi(A,b,x_old,eps):
	delta = 1.0
	n_iter = 0
	n = len(b)
	x_new = np.repeat(0.0,n)
	while delta > eps:
		n_iter += 1
		for i in range(0,n):
			x_new[i] = (b[i] - first_sum(i,A,x_old) - second_sum(i,A,x_old))/A[i,i]
		delta = max(abs(np.linalg.solve(A,b) - x_new))
		x_old = x_new
	return n_iter, x_new

x_max = 4.0
x_min = -4.0
D = 1.0
Sigma_a = 0.2
h = float(raw_input("Please enter mesh size [cm]: "))
case = str(raw_input("Please enter case (b or c): "))
L = sqrt(D/Sigma_a)
n = int((x_max - x_min)/h)
x_range = np.linspace(x_min,x_max,num=n+1)
eps = float(raw_input("Please enter absolute error tolerance: "))

A = np.matrix(np.zeros((n+1,n+1)))
for i in range(0,n+1):
        A[i,i] = 2 + (h/L)**2
        if i <= n-1:
                A[i+1,i] = -1.0; A[i,i+1] = -1.0
phi_old = np.repeat(0.0,n+1)
S = np.repeat(0.0,n+1)
S_matrix = np.repeat(0.0,n+1)
if case == "b":
	for i in range(0,n+1): S[i] = 8.0
elif case == "c":
	for i in range(0,n+1): S[i] = cos(x_range[i])
for i in range(0,n+1): S_matrix[i] = (h**2)*(S[i]/D)

n_iter = Jacobi(A,S_matrix,phi_old,eps)[0]

print("absolute error tolerance: "+str(eps))
print("mesh size: "+str(h))
print("number of iterations: "+str(n_iter))

h_plot = np.array([1.0,0.5,0.1,0.05])
e3b = np.array([38,131,2969,11733])
e5b = np.array([53,189,4297,16990])
e3c = np.array([20,58,1701,6828])
e5c = np.array([36,116,3029,12085])
plt.plot(h_plot, e3b, 'x', color = 'g', linestyle = '--')
plt.plot(h_plot, e5b, '.', color = 'g', linestyle = '--')
plt.plot(h_plot, e3c, 'x', color = 'b', linestyle = '--')
plt.plot(h_plot, e5c, '.', color = 'b', linestyle = '--')
plt.yscale('log')
plt.legend(("Case (b), $\epsilon = 10^{-3}$","Case (b), $\epsilon = 10^{-5}$","Case (c), $\epsilon = 10^{-3}$","Case (c), $\epsilon = 10^{-5}$"), loc="upper right")
plt.xlabel('Mesh Size [cm]'); plt.ylabel('Number of Iterations (log scale)')
plt.title("Number of Iterations vs. Mesh Size for Jacobi Method")
plt.savefig('hw6_4a.pdf')
