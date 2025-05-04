import math
import numpy as np
from scipy.stats import linregress
from collections import namedtuple


def Gran_F1(mass: list, emf: list, T: float, m0: float):
    # assumes data was not collected below a certain pH/emf
    # TODO maybe put in some good tools to find the optimal data range?
    Gran_data = namedtuple(
        "Gran_data",
        ["F1", "F1_mass", "slope", "intercept", "goodness_of_fit", "indices"],
    )
    k = 8.31451 * T / 96484.56
    F1_all_data = m0 * np.exp(emf / k)
    indices = (0, np.count_nonzero(F1_all_data > 100))
    F1 = F1_all_data[indices[0] : indices[1]]
    F1_mass = mass[indices[0] : indices[1]]
    slope, intercept, r_value, _, _ = linregress(F1_mass, F1)
    goodness_of_fit = r_value**2

    return Gran_data(F1, F1_mass, slope, intercept, goodness_of_fit, indices)
