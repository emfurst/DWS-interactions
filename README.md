# Scattering
Function package for generating theoretical light scattering models. Composed during the 2021 CHARM REU Program for the Furst Group.

<img src="https://github.com/nsbalbi/Scattering/blob/main/Examples/Results/Intensity%20Plot.png" width="60%">

## Usage

The functions mie_scattering and rayleigh_scattering calculate theoretical light scattering properties, including scattering intensities and photon mean-free path. Inputs are the system parameters (particle and medium refractive indeces, particle radii, incident light wavelength) along with the desired structure factor model (default S(q)=1). See the documentation within scattering.py for further detail. Examples are given below and example scripts can be found in the "Examples" folder.

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

<img src="https://github.com/nsbalbi/Scattering/blob/main/Examples/Results/S(q)%20Models.png" width="45%">

Another great resourse and similar code can be found here: https://miepython.readthedocs.io/en/latest/index.html
