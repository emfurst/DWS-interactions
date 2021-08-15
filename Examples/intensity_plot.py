"""
Plots Mie intensities and form factor versus scattering angle.
"""

from scattering import *
import matplotlib.pyplot as plt

# Default Scattering Parameters
a_p = 500e-9  # # particle radius [meters]
n_p = 1.5  # particle refractive index
n_s = 1.0  # medium refractive index
lambda_vac = 685e-9  # wavelength of light in vacuum [meters]
phi = 0.03  # particle volume fraction
n_ang = 100  # number of sampled angles

# Collect Data
theta, i1, i2, _, _ \
    = mie_scattering(n_p, n_s, a_p, lambda_vac, phi)
_, i1r, i2r, _, _ \
    = rayleigh_scattering(n_p, n_s, a_p, lambda_vac, phi)
theta = np.append(theta, np.pi + theta)
i1 = np.append(i1, i1[::-1])
i2 = np.append(i2, i2[::-1])
i1r = np.append(i1r, i1r[::-1])
i2r = np.append(i2r, i2r[::-1])

# Generate Colors
n_colors = 3
colormap = plt.get_cmap('plasma')
c = np.empty(shape=(n_colors, 3))
for i in range(n_colors):
    c[i, :] = colormap.colors[round(256 * i / n_colors)]

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw=dict(projection='polar'))
ax1.plot(theta, i1, label='Perpendicular Mie', c=c[0])
ax1.plot(theta, i2, label='Parallel Mie', c=c[2])
ax1.plot(theta, 0.5*(i1+i2), label='Mie P(q)', c=c[1])
# ax1.plot(theta, i1r, label='Perpendicular Rayleigh', c=c[0], linestyle='--')
# ax1.plot(theta, i2r, label='Parallel Rayleigh', c=c[2], linestyle='--')
# ax1.plot(theta, 0.5*(i1r+i2r), label='Rayleigh P(q)', c=c[1], linestyle='--')
ax1.set_title("Intensity")
ax1.legend()

ax2.plot(theta, np.log10(i1 + 1), label='Perpendicular Mie', c=c[0])
ax2.plot(theta, np.log10(i2 + 1), label='Parallel Mie', c=c[2])
ax2.plot(theta, np.log10(0.5*(i1+i2) + 1), label='Mie P(q)', c=c[1])
# ax2.plot(theta, np.log10(i1r + 1), label='Perpendicular Rayleigh', c=c[0], linestyle='--')
# ax2.plot(theta, np.log10(i2r + 1), label='Parallel Rayleigh', c=c[2], linestyle='--')
# ax2.plot(theta, np.log10(0.5*(i1r+i2r) + 1), label='Rayleigh P(q)', c=c[1], linestyle='--')
ax2.set_title("Log Intensity")
ax2.legend()

plt.suptitle(r'$n_p$ = %.3f, $n_s$ = %.3f, $a_p$ = %i nm, $\lambda$ = %i nm, $\phi$ = %.2f'
             % (n_p, n_s, a_p*1e9, lambda_vac*1e9, phi))
plt.show()
