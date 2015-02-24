import numpy as np

n = 5

A = np.matrix(np.zeros((n,n)))
for i in range(0,n):
        A[i,i] = 4.0
        if i <= n-2: A[i+1,i] = -1.0; A[i,i+1] = -1.0

b = np.repeat(100.0,n)
x_old = np.repeat(0.0,n)
x_new = np.repeat(0.0,n)
delta = 1.0
n_iter = 0

while delta > 10**(-6):
	n_iter += 1
	for i in range(0,n):
		sum1 = 0.0
		sum2 = 0.0
		for j in range(0,i):
			sum1 += A[i,j]*x_old[j]
		for j in range(i+1,n):
			sum2 += A[i,j]*x_old[j]
		x_new[i] = (b[i] - sum1 - sum2)/A[i,i]
	delta = max(abs(np.linalg.solve(A,b) - x_new))
	x_old = x_new

print n_iter
print x_new