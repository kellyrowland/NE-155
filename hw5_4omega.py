import matplotlib.pyplot as plt
import numpy as np

n = 5

A = np.matrix(np.zeros((n,n)))
for i in range(0,n):
        A[i,i] = 4.0
        if i <= n-2: A[i+1,i] = -1.0; A[i,i+1] = -1.0

b = np.repeat(100.0,n)
omega = np.linspace(0.5,1.975,1475)
n_iter = np.repeat(0,len(omega))

for i in range(len(omega)):
	delta = 1.0
	x_old = np.repeat(0.0,n)
	x_new = np.repeat(0.0,n)
	while delta > 10**(-6):
		n_iter[i] += 1
		for k in range(0,n):
			sum1 = 0.0
			sum2 = 0.0
			for j in range(0,k):
				sum1 += A[k,j]*x_new[j]
			for j in range(k+1,n):
				sum2 += A[k,j]*x_old[j]
			x_new[k] = (1 - omega[i])*x_old[k] + omega[i]*(b[k] - sum1 - sum2)/A[k,k]
	        delta = max(abs(np.linalg.solve(A,b) - x_new))
	        x_old = x_new

plt.plot(omega,n_iter)
plt.xlabel('$\omega$')
plt.ylabel('Number of Iterations')
plt.title('$\omega$ vs. Number of Iterations')
plt.axis([0.5, 2.0, 0, 700])
plt.savefig('hw5_omega.pdf')
