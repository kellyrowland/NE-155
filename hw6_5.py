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

def SOR(A,b,x_old,omega):
	delta = 1.0
	n_iter = 0
	n = len(b)
	x_new = np.repeat(0.0,n)
	while delta > 1E-6:
		n_iter += 1
		for i in range(0,n):
			x_new[i] = (1 - omega)*x_old[i] + omega*(b[i] - first_sum(i,A,x_new) - second_sum(i,A,x_old))/A[i,i]
		delta = max(abs(np.linalg.solve(A,b) - x_new))
		x_old = x_new
	return n_iter, x_new

x_max = 4.0
x_min = -4.0
D = 1.0
Sigma_a = 0.7
nuSig_f = 0.6
h = 0.1
L = sqrt(D/Sigma_a)
n = int((x_max - x_min)/h)
omega = 1.2
eps = 0.0001
n_iter = 0
converged = False

A = np.matrix(np.zeros((n,n)))
for i in range(0,n):
        A[i,i] = 2*D/(h**2) + Sigma_a
        if i <= n-2:
                A[i+1,i] = -D/(h**2); A[i,i+1] = -D/(h**2)

k_old = 1.0
k_new = 1.0
phi_old = np.array(np.repeat(1.0, n)/np.linalg.norm(np.repeat(1.0, n)))
phi_new = np.array(np.repeat(1.0, n)/np.linalg.norm(np.repeat(1.0, n)))
Q_old = nuSig_f*phi_old
Q_new = nuSig_f*phi_new

while(not converged):
	n_iter += 1
	phi_old = phi_new
	Q_old = nuSig_f*phi_old
	k_old = k_new
	phi_new = SOR(A,Q_old/k_old,phi_old,omega)[1]
	Q_new = nuSig_f*phi_new
	k_new = k_old*((0.5*h*sum(Q_new))/(0.5*h*sum(Q_old)))
	if(abs(k_new-k_old) < eps and max(abs(phi_new-phi_old)) < eps):
		converged = True

print("k: " + str(k_new))
print("number of iterations: " + str(n_iter))
phi_plot = phi_new/np.linalg.norm(phi_new)
x_range = np.linspace(x_min,x_max,num=n)
plt.plot(x_range, phi_plot)
plt.xlabel('x [cm]'); plt.ylabel('Flux [n/cc/s]')
plt.title("Power Iteration Flux Solution")
plt.savefig('hw6_5.pdf')
