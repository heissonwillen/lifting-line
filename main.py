import numpy as np
import matplotlib.pyplot as plt

# Inputs

# Wing geometry / aerodynamic
wing_span = 10
theta_t = 0.0
alpha_L_0 = 0.0
a0 = 2 * np.pi
c = 1
wing_area = c * wing_span
wing_AR = wing_span**2.0 / wing_area

# Flight condition
U_inf = 40  # m/s
alpha = 3 * np.pi / 180

# Discretization
N = 5
eps = wing_span * 0.02
theta = np.linspace(0+eps, np.pi-eps, N)

# Assembly of M matrix and b vector
M = np.zeros((N, N))
b = np.zeros(N)

for i in range(N):
    for j in range(N):
        M[i, j] = 4 * (wing_span / (a0 * c)) * np.sin((j+1) * theta[i]) + \
            (j+1) * np.sin((j+1) * theta[i]) / np.sin(theta[i])
    b[i] = alpha + theta_t - alpha_L_0

A = np.linalg.solve(M, b)

plt.plot(A, 'ob')

# Lift coeficient
CL = np.pi * wing_AR * A[0]

# Induced drag coeficient
delta = np.sum(np.arange(2, N+1, 1) * (A[1:] / A[0])**2.0)
CDi = (1 + delta) * CL**2.0/(np.pi * wing_AR)

print('Lift coefficient: ', CL)
print('Induced drag coefficient: ', CDi)

plt.figure()
Nvis = 300
y = np.linspace(-wing_span/2, wing_span/2, Nvis)
theta = np.arccos(-2*y/wing_span)

Gamma = np.zeros(Nvis)
for i in range(N):
    Gamma += A[i] * np.sin((i+1) * theta)
Gamma *= 2 * U_inf * wing_span

plt.plot(y, Gamma, '-b')
plt.xlabel('Spanwise position [m]')
plt.ylabel('Circulation')
plt.show()
