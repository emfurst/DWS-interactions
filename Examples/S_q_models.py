"""
Plots various structure factor models.
"""

from scattering import *
import numpy as np
import matplotlib.pyplot as plt

# Default Scattering Parameters
a_p = 200e-9  # meters
n_p = 1.5  # particle refractive index
n_s = 1.332  # medium refractive index (water)
lambda_vac = 685e-9  # wavelength of light in vacuum [meters]
phi = 0.03  # particle volume fraction
n_ang = 100  # number of sampled angles

s_q = np.empty(shape=(n_ang, 3))  # initialize result array

# Collect Scattering Data
_, _, _, _, [q, _, s_q[:, 0], _, _, _, _, _] = mie_scattering(n_p, n_s, a_p, lambda_vac, phi, n_ang=n_ang,
                                                              struct='PY')
_, _, _, _, [_, _, s_q[:, 1], _, _, _, _, _] = mie_scattering(n_p, n_s, a_p, lambda_vac, phi, n_ang=n_ang,
                                                              struct='PYAnnulus', optional_params=[1.2*a_p])
_, _, _, _, [_, _, s_q[:, 2], _, _, _, _, _] = mie_scattering(n_p, n_s, a_p, lambda_vac, phi, n_ang=n_ang,
                                                              struct='SHS', optional_params=[0.2, 0.05])

# Plot
fig, ax1 = plt.subplots()
ax1.plot(q, np.ones_like(s_q[:, 0]), 'k', label='No Interaction', linewidth=2)
ax1.plot(q, s_q[:, 0], 'k-.', label='Hard Sphere', linewidth=2)
ax1.plot(q, s_q[:, 1], 'k:', label='1.2$a_p$ Excluded Annulus', linewidth=2)
ax1.plot(q, s_q[:, 2], 'k--', label='Sticky Hard Sphere', linewidth=2)
ax1.set(xlabel=r'$q$', ylabel='S(q)',
        title=r'$n_p$ = %.3f, $n_s$ = %.3f, $a_p$ = %i nm, $\lambda$ = %i nm, $\phi$ = %.2f'
              % (n_p, n_s, a_p*1e9, lambda_vac*1e9, phi))
ax1.legend()

plt.show()
