import math
import numpy as np
from scipy.stats import linregress
from collections import namedtuple
from scipy.optimize import least_squares


def Gran_F1(mass: list, emf: list, T: float, m0: float):
    # assumes data was not collected below a certain pH/emf
    # TODO maybe put in some good tools to find the optimal data range?
    Gran_data = namedtuple(
        "Gran_data",
        ["F1", "F1_mass", "slope", "intercept", "goodness_of_fit", "indices"],
    )
    k = k_boltz(T)
    F1_all_data = m0 * np.exp(emf / k)
    indices = (0, np.count_nonzero(F1_all_data > 100))
    F1 = F1_all_data[indices[0] : indices[1]]
    F1_mass = mass[indices[0] : indices[1]]
    slope, intercept, r_value, _, _ = linregress(F1_mass, F1)
    goodness_of_fit = r_value**2

    return Gran_data(F1, F1_mass, slope, intercept, goodness_of_fit, indices)


def k_boltz(T: float):
    return 8.31451 * T / 96484.56


def find_data_in_range(low: float, high: float, data: iter):
    """
    Finds all indices in data that is in inclusive range

    Args:
        low (float): _description_
        high (float): _description_
        data (iter): _description_

    Returns:
        iter: _description_
    """
    return [i for i, x in enumerate(data) if low <= x <= high]


def estimate_AT_E0(
    titrant_mass: np.ndarray, emf: np.ndarray, T: float, m0: float, titrant_conc: float
) -> tuple[float, float]:
    gran_data = Gran_F1(titrant_mass, emf, T, m0)
    mass_eq = -gran_data.intercept / gran_data.slope
    k = k_boltz(T)
    titrant_moles = titrant_conc * (titrant_mass - mass_eq)
    E0_est = np.mean(emf - k * np.log(titrant_moles / (m0 + titrant_mass)))
    H_est = np.exp((emf - np.mean(E0_est)) / k)
    AT_est = mass_eq * titrant_conc / m0
    return AT_est, E0_est


def AT_residuals(x, sample, titration):
    ST = sample.ST
    KS = sample.KS
    FT = sample.FT
    KF = sample.KF
    m0 = sample.m0
    KW = sample.KW
    CHCl = titration.titrant.concentration
    m = titration.weight
    H = 10 ** -(titration.pH_est)
    f, AT = x
    Z = 1 + ST / KS
    HSO4 = m0 + ST / (1 + (Z * KS) / (f * H))
    HF = m0 * FT / (1 + KF / (f * H))
    residual = (
        m0 * AT + HSO4 + HF - m * CHCl + (m0 + m) * ((f * H / Z) - (Z * KW / (f * H)))
    )
    return residual
