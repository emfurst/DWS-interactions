"""
Function package for generating theoretical light scattering models.

Originally composed by Nicholas Sbalbi (nsbalbi@umass.edu) during the 2021 CHARM REU Program for the Furst Group in the
University of Delaware Department of Chemical and Biomolecular Engineering. Last updated 08/13/2021
"""

import numpy as np
import scipy.special as sps
import sys
from scattering_test import *


def mie_pi_tau(theta, alpha):
    """
    Calculate the functions pi and tau which appear in the Mie scattering calculation.

    These functions match those described by van de Hulst in Light Scattering by Small Particles, 1981, and are
    calculated using the method outlined by Kerker in The Scattering of Light and Other Electromagnetic Radiation, 1969.

    Args:
        theta: the desired scattering angles
        alpha: the size parameter of the spheres (2 * pi * radius / wavelength)

    Returns:
        pi: the first function
        tau: the second function
    """
    # Constants
    cos_theta = np.cos(theta)  # pre-calculated for efficiency
    sin_theta = np.sin(theta)  # "
    n_ang = np.size(theta)  # number of scattering angles
    n_terms = np.rint(alpha + 4.05*alpha**0.3333 + 2).astype(int)  # number of terms in series for coeff calculation

    # Initialize arrays
    pi = np.empty(shape=(n_terms, n_ang))
    pi_prime = np.empty(shape=(n_terms, n_ang))
    tau = np.empty(shape=(n_terms, n_ang))

    # Caluclation of pi and tau via Kerker's method
    for i in range(n_terms):  # loop over rows
        pi[i, :] = -sps.lpmv(1, i + 1, cos_theta) / sin_theta  # generates 1st order Legendre polynomials
        if i == 0:
            pi_prime[i, :] = np.zeros(shape=n_ang)
            tau[i, :] = cos_theta * pi[i, :]
        elif i == 1:
            pi_prime[i, :] = 3 * np.ones(shape=n_ang)
            tau[i, :] = cos_theta * pi[i, :] - 3 * sin_theta ** 2
        else:
            for j in range(n_ang):  # loop over columns
                pi_prime[i, j] = (2 * i + 1) * pi[i - 1, j] + pi_prime[i - 2, j]
            tau[i, :] = cos_theta * pi[i, :] - sin_theta ** 2 * pi_prime[i, :]

    return pi, tau


def mie_a_b(alpha, m):
    """
    Calculate the Mie coefficients for Mie scattering by spheres.

    The size of the arrays is based on the recommended by Barber and Hill in Light Scattering by Particles:
    Computational Methods, 1990. The coefficients are calculated using the method outlined by Kerker in
    The Scattering of Light and Other Electromagnetic Radiation, 1969.

    Args:
        alpha: the size parameter of the spheres (2 * pi * radius / wavelength)
        m: the ratio between the refractive indices of the spheres and the medium

    Returns:
        a: first set of Mie coefficients
        b: second set of Mie coefficients
    """
    # Constants
    beta = m * alpha
    n_terms = np.rint(alpha + 4.05*alpha**0.3333 + 2).astype(int)  # number of terms in series for coeff calculation
    x = np.arange(n_terms + 2) + 0.5  # two extra values for low and high bounds

    # Calculation of coefficients via Kerker's method
    psi = np.empty(shape=(2, n_terms + 2))  # add two extra terms for left and right bounds
    chi = np.empty(shape=(2, n_terms + 2))
    psi[0, :] = (np.pi * alpha / 2) ** 0.5 * sps.jv(x, alpha)  # generates Bessel function values of the first kind
    psi[1, :] = (np.pi * beta / 2) ** 0.5 * sps.jv(x, beta)  # generates Bessel function values of the second kind
    chi[0, :] = -(np.pi * alpha / 2) ** 0.5 * sps.yv(x, alpha)
    chi[1, :] = -(np.pi * beta / 2) ** 0.5 * sps.yv(x, beta)

    zeta = psi + 1j * chi  # conversion to complex numbers

    psi_prime = np.empty(shape=(2, n_terms))
    zeta_prime = np.empty(shape=(2, n_terms), dtype=np.complex64)
    psi_prime[0, :] = 0.5 * psi[0, :-2] - 0.5 * psi[0, 2:] + 0.5 / alpha * psi[0, 1:-1]
    psi_prime[1, :] = 0.5 * psi[1, :-2] - 0.5 * psi[1, 2:] + 0.5 / beta * psi[1, 1:-1]
    zeta_prime[0, :] = 0.5 * zeta[0, :-2] - 0.5 * zeta[0, 2:] + 0.5 / alpha * zeta[0, 1:-1]
    zeta_prime[1, :] = 0.5 * zeta[1, :-2] - 0.5 * zeta[1, 2:] + 0.5 / beta * zeta[1, 1:-1]

    a = (psi_prime[1, :] * psi[0, 1:-1] - m * psi_prime[0, :] * psi[1, 1:-1]) / \
        (psi_prime[1, :] * zeta[0, 1:-1] - m * zeta_prime[0, :] * psi[1, 1:-1])
    b = (m * psi_prime[1, :] * psi[0, 1:-1] - psi_prime[0, :] * psi[1, 1:-1]) / \
        (m * psi_prime[1, :] * zeta[0, 1:-1] - zeta_prime[0, :] * psi[1, 1:-1])

    return a, b


def mie_s1_s2(theta, alpha, m):
    """
    Calculate the scattering amplitude functions for Mie scattering by spheres.

    Args:
        theta: the angles at which scattering amplitudes are calculated
        alpha: the size parameter of the spheres (2 * pi * radius / wavelength)
        m: the ratio between the refractive indices of the spheres and the medium

    Returns:
        s1: perpendicular amplitude function at corresponding theta values
        s2: parallel amplitude function at corresponding theta values
        a: first set of Mie coefficients
        b: second set of Mie coefficients
    """
    # Constants
    n_ang = np.size(theta)
    n_terms = np.rint(alpha + 4.05*alpha**0.3333 + 2).astype(int)  # number of terms in series for coeff calculation

    # Calculation of coefficients using other functions
    pi, tau = mie_pi_tau(theta, alpha)
    a, b = mie_a_b(alpha, m)

    # Calculation of scattering amplitude functions
    s1 = np.reshape(np.tile(a, n_ang), (n_terms, n_ang), order='F') * pi \
        + np.reshape(np.tile(b, n_ang), (n_terms, n_ang), order='F') * tau  # perpendicular amplitude function
    s2 = np.reshape(np.tile(a, n_ang), (n_terms, n_ang), order='F') * tau \
        + np.reshape(np.tile(b, n_ang), (n_terms, n_ang), order='F') * pi  # parallel amplitude function
    temp = np.arange(n_terms) + 1  # factor multiplied by each term (row)
    temp = (2 * temp + 1) / (temp ** 2 + temp)
    s1 = s1 * np.reshape(np.tile(temp, n_ang), (n_terms, n_ang), order='F')
    s2 = s2 * np.reshape(np.tile(temp, n_ang), (n_terms, n_ang), order='F')
    s1 = np.sum(s1, axis=0)  # sum along angle (column)
    s2 = np.sum(s2, axis=0)

    return s1, s2, a, b


def percus_yevick(qa, phi):
    """
    Calculate the structure factor for hard spheres using the Percus-Yevick
    closure of the Ornstein-Zernicke equation.

    Args:
        qa: the scattering vector scaled by the sphere radius at the desired scattering angles
        phi: volume fraction of spheres in solution

    Returns:
        s_q: the structure factor at the corresponding scattering angles
    """
    qa2 = 2 * qa  # twice qa, for efficiency/conciseness
    lambda_1 = (1 + 2 * phi) ** 2 / (1 - phi) ** 4
    lambda_2 = -(1 + phi / 2) ** 2 / (1 - phi) ** 4
    c_1 = lambda_1 / qa2 ** 3 * (np.sin(qa2) - qa2 * np.cos(qa2))
    c_2 = -6 * phi * lambda_2 / qa2 ** 4 * (qa2 ** 2 * np.cos(qa2) - 2 * qa2 * np.sin(qa2) - 2 * np.cos(qa2) + 2)
    c_3 = -phi * lambda_1 / 2 / qa2 ** 6 * \
        (qa2 ** 4 * np.cos(qa2) - 4 * qa2 ** 3 * np.sin(qa2) - 12 * qa2 ** 2 * np.cos(qa2)
         + 24 * qa2 * np.sin(qa2) + 24 * np.cos(qa2) - 24)

    s_q = np.divide(1, 1 + 24 * phi * (c_1 + c_2 + c_3))
    return s_q


def percus_yevick_annulus(q, a_p, a_e, phi):
    """
    Calculate the structure factor for hard spheres with an exluded annulus
    using the Percus-Yevick closure of the Ornstein-Zernicke equation.

    Args:
        q: scattering vector at desired scattering angles
        a_p: sphere radius (excluding annulus)
        a_e: effective sphere radius (including annulus)
        phi: volume fraction of spheres in solution

    Returns:
        s_q: the structure factor at the corresponding scattering angles
    """
    phi_e = phi * (a_e / a_p) ** 3  # effective volume fraction

    qa2 = 2 * q * a_e  # twice qa, for efficiency/conciseness
    lambda_1 = (1 + 2 * phi_e) ** 2 / (1 - phi_e) ** 4
    lambda_2 = -(1 + phi_e / 2) ** 2 / (1 - phi_e) ** 4
    c_1 = lambda_1 / qa2 ** 3 * (np.sin(qa2) - qa2 * np.cos(qa2))
    c_2 = -6 * phi_e * lambda_2 / qa2 ** 4 * (qa2 ** 2 * np.cos(qa2) - 2 * qa2 * np.sin(qa2) - 2 * np.cos(qa2) + 2)
    c_3 = -phi_e * lambda_1 / 2 / qa2 ** 6 * \
        (qa2 ** 4 * np.cos(qa2) - 4 * qa2 ** 3 * np.sin(qa2) - 12 * qa2 ** 2 * np.cos(qa2)
         + 24 * qa2 * np.sin(qa2) + 24 * np.cos(qa2) - 24)

    s_q = np.divide(1, 1 + 24 * phi_e * (c_1 + c_2 + c_3))
    return s_q


def sticky_hard_sphere(q, a_p, phi, tau=0.2, epsilon=0.05):
    """
    Calculate the structure factor for hard spheres with square-well interparticle
    potential. Based on the method of Rao et. al. in "A new interpretation of the
    sticky hard sphere model", J. Chem. Phys. 95(12), 9186-9190 (1991).

    Args:
        q: scattering vector at desired scattering angles
        a_p: sphere radius (excluding annulus)
        phi: volume fraction of spheres in solution
        tau: float, "stickiness" parameter within {0, 1}
        epsilon: float, perturbation parameter within {0, 0.1}

    Returns:
        s_q: the structure factor at the corresponding scattering angles
    """
    eta = phi / (1 - epsilon) ** 3
    eta1 = 1 - eta
    sigma = 2 * a_p
    a = sigma / (1 - epsilon)
    kappa = q * a

    # Quadratic for Lambda
    a_q = eta/12
    b_q = -tau - eta / eta1
    c_q = (1 + eta / 2) / eta1 ** 2

    disc = b_q ** 2 - 4 * a_q * c_q  # discriminant
    if disc < 0:  # error if no real roots
        print('Error: No Real Roots within SHS Model')
        return -1*np.ones_like(q)
    lamb = (-b_q + np.sqrt(disc)) / (2 * a_q)
    lamb2 = (-b_q - np.sqrt(disc)) / (2 * a_q)
    if lamb2 < lamb:
        lamb = lamb2  # take smaller root, larger is unphysical

    # Alpha and Beta
    mu = lamb * eta * eta1
    if mu > 1 + 2 * eta:  # model condition set in paper
        print('Error: Unphysical SHS Model Result')
        return -1*np.ones_like(q)
    alpha = (1 + 2 * eta - mu) / eta1 ** 2
    beta = (-3 * eta + mu) / 2 / eta1 ** 2

    # S_q
    k2 = kappa ** 2
    k3 = kappa ** 3
    ksin = np.sin(kappa)
    kcos = np.cos(kappa)
    A = 1 + 12 * eta * (alpha * (ksin - kappa * kcos) / k3 + beta * (1 - kcos) / k2 - lamb / 12 * ksin / kappa)
    B = 12 * eta * (alpha * (1 / 2 / kappa - ksin / k2 + (1 - kcos) / k3)
                    + beta * (1 / kappa - ksin / k2) - lamb / 12 * (1 - kcos) / kappa)
    s_q = 1 / (A ** 2 + B ** 2)

    return s_q


def mie_scattering(n_p, n_s, a_p, lambda_vac, phi, n_ang=100, struct='none', optional_params=None):
    """
    Calculate intensities and statistics for Mie scattering by spheres.

    Methods and algorithms are those outlined by Kerker in
    The Scattering of Light and Other Electromagnetic Radiation, 1969

    Args:
        n_p: float, refractive index of the spheres
        n_s: float, refractive index of the medium
        a_p: float, sphere radius [m]
        lambda_vac: float, wavelength of light in vacuum [m]
        phi: float, volume fraction of particles in solution
        n_ang: int, number of angles calculated (Default = 100)
        struct: str, choice of equation used for structure factor
            - 'none': S(q) = 1
            - 'PY': S(q) is calculated using the Percus-Yevick approximation
            - 'PYAnnulus': S(q) is calculated using the Percus-Yevick approximation with an excluded annulus
            - 'SHS': S(q) is calculated via a perturbative sticky hard sphere solution of the Percus-Yevick closure
        optional_params: list of optional parameters required for some structure factor calculations
            - 'PYAnnulus': [a_e]
                - a_e: float, effective sphere radius (including excluded annulus)
            - 'SHS': [tau, epsilon]
                - tau: float, "stickiness" parameter within {0, 1}
                - epsilon: float, perturbation parameter within {0, 0.1}

    Returns:
        theta: the angles at which scattering intensities are calculated
        i1: perpendicular scattering intensity at corresponding scattering angles
        i2: parallel scattering intensity
        l_star: photon mean free path
        stats: list of other output scattering statistics
            - q: array of scattering vectors corresponding with theta
            - p_q: form factor at corresponding theta/q values
            - s_q: structure factor at corresponding theta/q values
            - l: scattering mean free path
            - q_scat: scattering efficiency factor
            - q_ext: extinction efficiency factor
            - c_scat: scattering cross section
            - c_ext: extinction cross section
    """
    # Initial Constants
    m = n_p / n_s  # refractive index ratio
    lambda_med = lambda_vac / n_s  # wavelength in medium
    k_0 = 2 * np.pi / lambda_med  # scattering wavevector
    alpha = k_0 * a_p
    rho = phi / 4 * 3 / np.pi / a_p ** 3  # number density of scatterers
    n_terms = np.rint(alpha + 4.05*alpha**0.3333 + 2).astype(int)  # number of terms in series for coeff calculation

    theta = np.linspace(0.001, 0.999 * np.pi, n_ang)  # angles (ends are not absolute to prevent zeroing out)
    qa = a_p * 2 * k_0 * np.sin(theta/2)  # radius times scattering vector

    # Intensity Calculation
    s1, s2, a, b = mie_s1_s2(theta, alpha, m)

    i1 = np.real(s1 * np.conjugate(s1))  # perpendicular intensity
    i2 = np.real(s2 * np.conjugate(s2))  # parallel intensity

    # Scattering Statistics Calculation
    temp = 2 * np.arange(n_terms) + 3  # factor multiplied by each term (row)
    # scattering efficiency factor
    q_scat = 2 / alpha ** 2 * np.sum(temp * np.real(a * np.conjugate(a) + b * np.conjugate(b)))
    q_ext = 2 / alpha ** 2 * np.sum(temp * np.real(a + b))  # extinction efficiency factor
    c_scat = q_scat * 3.14158 * a_p ** 2  # scattering cross section
    c_ext = q_ext * 3.14158 * a_p ** 2  # extinction cross section

    # Structure and Form Factor Calculation
    # structure factor
    if struct == 'none':
        s_q = np.ones_like(theta)
    elif struct == 'PY':
        s_q = percus_yevick(qa, phi)
    elif struct == 'PYAnnulus':
        s_q = percus_yevick_annulus(qa / a_p, a_p, optional_params[0], phi)
    elif struct == 'SHS':
        s_q = sticky_hard_sphere(qa / a_p, a_p, phi, optional_params[0], optional_params[1])
    else:
        sys.exit('Error: Unexpected Input for Choice of Structure Factor Equation')
    p_q = (i1 + i2) / 2  # form factor (average intensity)

    # l Calculation
    integrand = qa * p_q * s_q
    l = (k_0 ** 4 * a_p ** 2) / \
        (2 * np.pi * rho * 0.5 * np.sum((integrand[:-1] + integrand[1:]) * (qa[1:] - qa[:-1])))

    # l* Calculation
    integrand = qa ** 3 * p_q * s_q
    l_star = (k_0 ** 6 * a_p ** 4) / \
        (np.pi * rho * 0.5 * np.sum((integrand[:-1] + integrand[1:]) * (qa[1:] - qa[:-1])))

    # Gather Misc Statistics
    q = 2 * k_0 * np.sin(theta / 2)  # scattering vector
    stats = [q, p_q, s_q, l, q_scat, q_ext, c_scat, c_ext]

    return theta, i1, i2, l_star, stats


def rayleigh_scattering(n_p, n_s, a_p, lambda_vac, phi, n_ang=100, struct='none', optional_params=None):
    """
    Calculate intensities and statistics for Rayleigh scattering by spheres.

    Args:
        n_p: float, refractive index of the spheres
        n_s: float, refractive index of the medium
        a_p: float, sphere radius [m]
        lambda_vac: float, wavelength of light in vacuum [m]
        phi: float, volume fraction of particles in solution
        n_ang: int, number of angles calculated (Default = 100)
        struct: str, choice of equation used for structure factor
            - 'none': S(q) = 1
            - 'PY': S(q) is calculated using the Percus-Yevick approximation
            - 'PYAnnulus': S(q) is calculated using the Percus-Yevick approximation with an excluded annulus
        optional_params: list of optional parameters required for some structure factor calculations
            - 'PYAnnulus': [a_e]
                - a_e: float, effective sphere radius (including excluded annulus)

    Returns:
        theta: the angles at which scattering intensities are calculated
        i1: perpendicular scattering intensity at corresponding scattering angles
        i2: parallel scattering intensity
        l_star: photon mean free path
        stats: list of other output scattering statistics
            - q: array of scattering vectors corresponding with theta
            - p_q: form factor at corresponding theta/q values
            - s_q: structure factor at corresponding theta/q values
            - l: scattering mean free path
    """
    # Initial Constants
    m = n_p / n_s  # refractive index ratio
    lambda_med = lambda_vac / n_s  # wavelength in medium
    k_0 = 2 * np.pi / lambda_med  # scattering wavevector
    rho = phi / 4 * 3 / np.pi / a_p ** 3  # number density of scatterers
    v = 4 / 3 * np.pi * a_p ** 3  # sphere volume

    theta = np.linspace(0.001, 0.999 * np.pi, n_ang)  # angles (ends are not absolute to prevent zeroing out)
    cos_theta = np.cos(theta)
    qa = a_p * 2 * k_0 * np.sin(theta/2)  # radius times scattering vector

    # Intensity Calculation
    i1 = k_0 ** 6 * v ** 2 * ((m - 1) / 2 / np.pi) ** 2 \
        * np.divide(3 * (np.sin(qa) - qa * np.cos(qa)), qa ** 3) ** 2  # perpendicular Rayleigh  intensity
    i2 = i1 * cos_theta ** 2  # parallel Rayleigh  intensity

    # Form and Structure Factor Calculation
    # structure factor
    if struct == 'none':
        s_q = np.ones_like(theta)
    elif struct == 'PY':
        s_q = percus_yevick(qa, phi)
    elif struct == 'PYAnnulus':
        s_q = percus_yevick_annulus(qa / a_p, a_p, optional_params[0], phi)
    else:
        sys.exit('Error: Unexpected Input for Choice of Structure Factor Equation')
    p_q = (i1 + i2) / 2  # form factor (average intensity)

    # l Calculation
    integrand = qa * p_q * s_q
    l = (k_0 ** 4 * a_p ** 2) / \
        (2 * np.pi * rho * 0.5 * np.sum((integrand[:-1] + integrand[1:]) * (qa[1:] - qa[:-1])))

    # l* Calculation
    integrand = qa ** 3 * p_q * s_q
    l_star = (k_0 ** 6 * a_p ** 4) / \
        (np.pi * rho * 0.5 * np.sum((integrand[:-1] + integrand[1:]) * (qa[1:] - qa[:-1])))

    # Gather Misc Statistics
    q = 2 * k_0 * np.sin(theta / 2)  # scattering vector
    stats = [q, p_q, s_q, l]

    return theta, i1, i2, l_star, stats
