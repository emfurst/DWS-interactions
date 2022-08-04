# DWS
Function package for generating light transport properties in multiple scattering media  
Nicholas Sbalbi, Qi Li, and Eric M. Furst

<img src="https://github.com/emfurst/DWS-interactions/blob/main/Examples/Results/Intensity%20Plot.png" width="60%">

## Usage

The functions `mie_scattering` and `rayleigh_scattering` calculate light scattering properties for scatterers, including scattering intensities and photon mean-free path. The input parameters are the particle and medium refractive indices, particle radius, and incident light wavelength, along with the desired structure factor model. The default is to assume non-interacting particles by setting the structure factor S(q)=1. See the documentation within `scattering.py` for further details. Examples are given below and example scripts can be found in the `Examples` folder.

## Examples

In the simplest case, only the following system parameters meed to be input:
* particle refractive index, n_p
* medium refractive index, n_s
* particle radii, a_p
* wavelength of the incident light in a vacuum, lambda_vac
* particle volume fraction, phi

```python
theta, i1, i2, l_star, stats = mie_scattering(n_p, n_s, a_p, lambda_vac, phi)
```

Plotting i1 and i2 versus theta generated the polar plots seen above.

To specify the structure factor, the argument "struct" can be used. Currently implimented is a Percus-Yevick hard sphere structure factor ('PY'), a hard sphere with an excluded annulus ('PYAnnulus'), and a sticky hard sphere model with a square potential well ('SHS'). For some models, additional parameters must be specified using the "optional_params" arguement. Examples of each are seen below along with an plot of example parameters. See the documentation of scattering.py for parameter details and citations.

```python
_, _, _, _, [q, _, s_q, _, _, _, _, _] = mie_scattering(n_p, n_s, a_p, lambda_vac, phi, struct='PY')
...  = mie_scattering(n_p, n_s, a_p, lambda_vac, phi, sruct='PYAnnulus', optional_params=[1.2*a_p])
...  = mie_scattering(n_p, n_s, a_p, lambda_vac, phi, sruct='SHS', optional_params=[0.2, 0.05])
```

<img src="https://github.com/emfurst/DWS-interactions/blob/main/Examples/Results/S(q)%20Models.png" width="45%">

The code is sufficiently optimized so that 100+ calls to mie/rayleigh_scattering can occur in under a second. See mie_vs_rayleigh_l_star.py for an example; the results are shown below.

<img src="https://github.com/emfurst/DWS-interactions/blob/main/Examples/Results/Mie%20vs%20Rayleigh%20l_star.png" width="90%">

## Other Resources

Other resources for scattering calculations can be found here: 
* https://miepython.readthedocs.io/en/latest/index.html
* https://omlc.org/calc/mie_calc.html

## License

[GNU General Public License v3.0](https://choosealicense.com/licenses/gpl-3.0/)

## Acknoledgments

The development of this code was primarily supported by NSF through the University of Delaware Materials Research Science and Engineering Center DMR-2011824. Additional financial support was received from the Chemours Company, Wilmington, DE. 
