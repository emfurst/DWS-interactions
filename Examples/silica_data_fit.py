"""
Generates l* values versus particle volume concentration for silica particles with varying form factor models.
Plots model versus experimental data (courtesy of Qi Li).
"""

from scattering import *
import numpy as np
import matplotlib.pyplot as plt

# Scattering Parameters
n_p = 1.446  # particle refractive index (silica)
n_s = 1.332  # medium refractive index (water)
a_p = 345e-9/2  # particle radii [meters]
lambda_vac = 685e-9  # wavelength of light in vacuum [meters]

# Generate Parameter Range
n_range = 100  # number of data points between values
phi_range = np.linspace(0.03, 0.15, n_range)

# Initialize Result Array
l_star_phi = np.empty(shape=(n_range, 4))

# Silica Data (courtesy of Qi Li)
silica_data = np.array([[0.065, 0.067, 0.069, 0.072, 0.075, 0.135],
                        [0.000230, 0.000225, 0.000216, 0.000215, 0.000205, 0.000130]])

# Collect Scattering Data
for i in range(n_range):
    _, _, _,  l_star_phi[i, 0], _ \
        = mie_scattering(n_p, n_s, a_p, lambda_vac, phi_range[i])
    _, _, _, l_star_phi[i, 1], _ \
        = mie_scattering(n_p, n_s, a_p, lambda_vac, phi_range[i], struct='PY')
    _, _, _, l_star_phi[i, 2], _ \
        = mie_scattering(n_p, n_s, a_p, lambda_vac, phi_range[i], struct='SHS', optional_params=[0.2, 0.05])
    _, _, _, l_star_phi[i, 3], _ \
        = mie_scattering(n_p, n_s, a_p, lambda_vac, phi_range[i], struct='PYAnnulus', optional_params=[1.2*a_p])

# Generate Colors
n_colors = 4
colormap = plt.get_cmap('plasma')
c = np.empty(shape=(n_colors, 3))
for i in range(n_colors):
    c[i, :] = colormap.colors[round(256 * i / n_colors)]

# Plot Data
fig, ax1 = plt.subplots()
ax1.plot(phi_range, l_star_phi[:, 0]*1e6, label='No Interactions', linewidth=2, c=c[0])
ax1.plot(phi_range, l_star_phi[:, 1]*1e6, linestyle='--', label='Hard Sphere', linewidth=2, c=c[1])
ax1.plot(phi_range, l_star_phi[:, 2]*1e6, linestyle=':', label='Sticky Hard Sphere', linewidth=2, c=c[2])
ax1.plot(phi_range, l_star_phi[:, 3]*1e6, linestyle='-.', label='1.2$a_p$ Excluded Annulus', linewidth=2, c=c[3])
ax1.plot(silica_data[0, :], silica_data[1, :]*1e6, 'k.', label='Silica Data', markersize=10)
ax1.legend()
ax1.set(xlabel=r'$\phi$', ylabel='l* [µm]')

plt.show()
