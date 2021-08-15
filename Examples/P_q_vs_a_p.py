"""
Plots Mie and Rayleigh form factors at varying particle radii.
"""

from scattering import *
import numpy as np
import matplotlib.pyplot as plt

# Default Scattering Parameters
n_p = 2.0  # particle refractive index
n_s = 1.332  # medium refractive index (water)
lambda_vac = 685e-9  # wavelength of light in vacuum [meters]
phi = 0.03  # particle volume fraction
n_ang = 100  # number of sampled angles

a_p_range = [50e-9, 100e-9, 250e-9, 500e-9, 1000e-9]  # particle radii [meters]
p_q = np.empty(shape=(n_ang, np.size(a_p_range), 2))  # initialize result array

# Collect Scattering Data
for i in range(np.size(a_p_range)):
    _, _, _, _, [q, p_q[:, i, 0], _, _, _, _, _, _] = \
        mie_scattering(n_p, n_s, a_p_range[i], lambda_vac, phi, n_ang=n_ang, struct='PY')
    _, _, _, _, [_, p_q[:, i, 1], _, _] = \
        rayleigh_scattering(n_p, n_s, a_p_range[i], lambda_vac, phi, n_ang=n_ang)

# Generate Colors
colormap = plt.get_cmap('plasma')
c = np.empty(shape=(np.size(a_p_range), 3))
for i in range(np.size(a_p_range)):
    c[i, :] = colormap.colors[round(256 * i / np.size(a_p_range))]

# Plot
fig, ax1 = plt.subplots()
for i in range(np.size(a_p_range)):
    if i == 0:
        ax1.plot(q, p_q[:, i, 0] / p_q[0, i, 0], label=r'Mie, $a_p$ = %i nm' % (a_p_range[i] * 1e9), c=c[i])
        ax1.plot(q, p_q[:, i, 1] / p_q[0, i, 1], linestyle='--',
                 label=r'Rayleigh, $a_p$ = %i nm' % (a_p_range[i] * 1e9), c=c[i])
    else:
        ax1.plot(q, p_q[:, i, 0] / p_q[0, i, 0], label=r'$a_p$ = %i nm' % (a_p_range[i] * 1e9), c=c[i])
        ax1.plot(q, p_q[:, i, 1] / p_q[0, i, 1], linestyle='--', c=c[i])
ax1.set(xlabel='q', ylabel='Normalized P(q)', title=r'$n_p$ = %.3f, $n_s$ = %.3f, $\lambda$ = %i nm, $\phi$ = %.2f' % (n_p, n_s, lambda_vac*1e9, phi))
ax1.legend(loc='upper right')
ax1.set(ylim=[-0.05, 1.1])
plt.show()

