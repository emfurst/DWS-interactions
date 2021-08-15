"""
Plots l* for varying input parameters for both Mie and Rayleigh scattering.
"""

from scattering import *
import numpy as np
import matplotlib.pyplot as plt

# Default Scattering Parameters
n_p = 1.8  # particle refractive index
n_s = 1.332  # medium refractive index (water)
a_p = 250e-9  # particle radii [meters]
lambda_vac = 685e-9  # wavelength of light in vacuum [meters]
phi = 0.03  # particle volume fraction

# Generate Parameter Ranges
n_range = 100  # number of data points between values
phi_range = np.linspace(0.01, 0.15, n_range)
a_p_range = np.linspace(150e-9/2, 1000e-9/2, n_range)
n_p_range = np.linspace(1.6, 2, n_range)

# Initialize Result Array
l_star = np.empty(shape=(n_range, 3, 2))

# Collect Scattering Data
for i in range(n_range):
    _, _, _, l_star[i, 0, 0], _ \
        = mie_scattering(n_p, n_s, a_p, lambda_vac, phi_range[i])
    _, _, _, l_star[i, 1, 0], _ \
        = mie_scattering(n_p, n_s, a_p_range[i], lambda_vac, phi)
    _, _, _, l_star[i, 2, 0], _ \
        = mie_scattering(n_p_range[i], n_s, a_p, lambda_vac, phi)

    _, _, _, l_star[i, 0, 1], _ \
        = rayleigh_scattering(n_p, n_s, a_p, lambda_vac, phi_range[i])
    _, _, _, l_star[i, 1, 1], _ \
        = rayleigh_scattering(n_p, n_s, a_p_range[i], lambda_vac, phi)
    _, _, _, l_star[i, 2, 1], _ \
        = rayleigh_scattering(n_p_range[i], n_s, a_p, lambda_vac, phi)

# Plot Data
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
ax1.plot(phi_range, l_star[:, 0, 0]*1e6, 'k', label='Mie', linewidth=2)
ax1.plot(phi_range, l_star[:, 0, 1]*1e6, 'k--', label='Rayleigh', linewidth=2)
ax1.legend()
ax1.set(xlabel=r'$\phi$', ylabel='l* [µm]')

ax2.plot(a_p_range*1e9, l_star[:, 1, 0]*1e6, 'k', label='Mie', linewidth=2)
ax2.plot(a_p_range*1e9, l_star[:, 1, 1]*1e6, 'k--', label='Rayleigh', linewidth=2)
ax2.set(xlabel=r'$a_p$ [nm]')

ax3.plot(n_p_range, l_star[:, 2, 0]*1e6, 'k', label='Mie', linewidth=2)
ax3.plot(n_p_range, l_star[:, 2, 1]*1e6, 'k--', label='Rayleigh', linewidth=2)
ax3.set(xlabel=r'$n_p$')

# title with default parameters listed
plt.suptitle(r'$n_p$ = %.3f, $n_s$ = %.3f, $a_p$ = %i nm, $\lambda$ = %i nm, $\phi$ = %.2f' % (n_p, n_s, a_p*1e9, lambda_vac*1e9, phi))
plt.show()
