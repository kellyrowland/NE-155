import numpy as np
import matplotlib.pyplot as plt

h = np.array([1.0,0.5,0.1,0.05,0.01])
b = np.array([3636363.63636,975609.756098,39960.03996,9997.50062484,399.99600004])
c = np.array([297110.736756,79712.6366907,3264.95315117,816.850313501,32.6818542246])

plt.plot(h, b, 'x', linestyle = '--')
plt.plot(h, c, '.', linestyle = '--')
plt.yscale('log')
plt.xscale('log')
plt.legend(("Case (b)","Case (c)"), loc="upper left")
plt.xlabel('Mesh Size [cm] (log scale)'); plt.ylabel('Maximum Relative Error (log scale)')
plt.title('Maximum Relative Error vs. Mesh Size')
plt.savefig('hw6_3.pdf')
