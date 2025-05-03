# base class for any solution
from statistics import mean
from math import exp, log
from typing import Union
import numpy as np


# TODO the point of the properties here is that some vlaues can be (re)calculated on the fly when the nsolution changes,
# including temperature or adding stuff so that volume and concentartion and ionic strength/salinity changes
class Solution:

    # everything is estimated at 20 deg. C
    def __init__(
        self,
        type: str = "solution",
        I: float = 0,
        temp: float = 20,
        id: str = None,
    ):
        self.type = type
        self.I = I
        self.id = id
        if temp > 200:
            self.T = temp
            self.t = self.T - 273.15
        else:
            self.t = temp
            self.T = self.t + 273.15

    @property
    def k(self):
        return 8.31451 * self.T / 96484.56

    @property
    def K1(self):
        # i is the ionic strength that is presumed to be a property of the solution using this version of K1
        # by necessity on H+(free)
        m = self.I * 1000 / (1000 - self.I * (22.99 + 35.45))
        A = 35.2911 * m**0.5 + 0.8491 * m - 0.32 * m**1.5 + 0.055 * m**2
        B = -1583.09 * m**0.5
        C = -5.4366 * m**0.5
        pK1 = (
            -402.56788
            + 11656.46 / self.T
            + 72.173 * log(self.T)
            - 0.161325 * self.T
            + 7.5526e-5 * self.T**2
        )
        return 10 ** -(A + B / self.T + C * log(self.T) + pK1)

    @property
    def K2(self, T=293):
        # i is the ionic strength that is presumed to be a property of the solution using this version of K1
        # by necessity on H+(free)
        m = self.I * 1000 / (1000 - self.I * (22.99 + 35.45))
        A = 38.2746 * m**0.5 + 1.6057 * m - 0.647 * m**1.5 + 0.113 * m**2
        B = -1738.16 * m**0.5
        C = -6.0346 * m**0.5
        pK2 = -122.4994 + 5811.18 / self.T + 20.5263 * log(self.T) - 0.0120897 * self.T
        return 10 ** -(A + B / self.T + C * log(self.T) + pK2)


class NaCl(Solution):

    def __init__(self, concentration: float = None):
        self.c = concentration
        self.I = self.c  # calculate from concentration
        super().__init__("NaCl")
        # assume these have not been added
        self.ST = 0
        self.FT = 0
        self.BT = 0

    @property
    def KW(self):
        p00 = 14.83
        p10 = -0.4914
        p01 = -0.0471
        p20 = 0.3917
        p11 = 0.0001724
        p02 = 0.0004381
        p30 = -0.07454
        p21 = 3.419e-06
        p12 = -1.081e-05
        p03 = -4.382e-06
        return 10 ** -(
            p00
            + p10 * self.c
            + p01 * self.t
            + p20 * self.c**2
            + p11 * self.c * self.t
            + p02 * self.t**2
            + p30 * self.c**3
            + p21 * self.c**2 * self.t
            + p12 * self.c * self.t**2
            + p03 * self.t**3
        )


class KCl(Solution):

    def __init__(self, concentration: float = None):
        self.c = concentration
        self.I = self.c  # calculate from concentration
        super().__init__("KCl")
        # assume these have not been added
        self.ST = 0
        self.FT = 0
        self.BT = 0

    @property
    def KW(self):
        p00 = 14.86
        p10 = -1.062
        p01 = -0.0479
        p20 = 2.031
        p11 = 0.0005349
        p02 = 0.000483
        p30 = -1.095
        p21 = -0.0005724
        p12 = -9.838e-06
        p03 = -5.091e-06
        return 10 ** -(
            p00
            + p10 * self.c
            + p01 * self.t
            + p20 * self.c**2
            + p11 * self.c * self.t
            + p02 * self.t**2
            + p30 * self.c**3
            + p21 * self.c**2 * self.t
            + p12 * self.c * self.t**2
            + p03 * self.t**3
        )


class SW(Solution):

    def __init__(self, salinity: float = 35):
        self.S = salinity
        self.I = self.S  # calculate from concentration
        super().__init__("SW")

    @property
    def KW(self):
        # used with H+(tot) --> not the original equation
        return exp(
            -13847.26 / self.T
            + 148.9652
            - 23.652 * log(self.T)
            + (118.67 / self.T - 5.977 + 1.0495 * log(self.T)) * self.S**0.5
            - 0.01615 * self.S
        )

    @property
    def K1(self):
        # used with H+(tot)
        return 10 ** (
            -3633.86 / self.T
            + 61.2172
            - 9.67770 * log(self.T)
            + 0.011555 * self.S
            - 0.0001152 * self.S**2
        )

    @property
    def K2(self):
        # used with H+(tot)
        return 10 ** (
            -471.78 / self.T
            - 25.9290
            + 3.16967 * log(self.T)
            + 0.01781 * self.S
            - 0.0001122 * self.S**2
        )

    @property
    def ST(self):
        return 0.14 / 96.062 * self.S / 1.80655

    @property
    def FT(self):
        return 0.000067 / 18.998 * self.S / 1.80655

    @property
    def BT(self, BTrat: str = "uppstrom", ratio: float = 0.5):
        # ratio will favor Lee, i.e., higher ratio will weight Lee BT/S heavier
        if ratio > 1:
            raise Exception("Ratio for BT/S cannot be above 1.")
        if BTrat.lower() == "uppstrom":
            return self.BT_uppstrom()
        elif BTrat.lower() == "lee":
            return self.BT_lee()
        elif BTrat.lower() == "mix":
            return self.BT_lee() * ratio + self.BT_uppstrom() * (1 - ratio)

    def BT_uppstrom(self):
        return 0.000232 / 10.811 * self.S / 1.80655

    def BT_lee(self):
        return 0.0002414 / 10.811 * self.S / 1.80655

    @property
    def KS(self):
        # used with H+(free)
        return exp(
            -4276.1 / self.T
            + 141.328
            - 23.093 * log(self.T)
            + (-13856 / self.T + 324.57 - 47.986 * log(self.T)) * self.I**0.5
            + (35474 / self.T - 771.54 + 114.723 * log(self.T)) * self.I
            - 2698.0 / self.T * self.I**1.5
            + 1776.0 / self.T * self.I**2
            + log(1 - 0.001005 * self.S)
        )

    @property
    def KF(self):
        # used with H+(tot)
        return exp(874.0 / self.T - 9.68 + 0.111 * self.S**0.5)

    @property
    def KB(self):
        # used with H+(tot)
        return exp(
            (
                -8966.90
                - 2890.53 * self.S**0.5
                - 77.942 * self.S
                + 1.728 * self.S**1.5
                - 0.0996 * self.S**2
            )
            / self.T
            + 148.0248
            + 137.1942 * self.S**0.5
            + 1.62142 * self.S
            + (-24.4344 - 25.085 * self.S**0.5 - 0.2474 * self.S) * log(self.T)
            + 0.053105 * self.S**0.5 * self.T
        )

    @property
    def KSi(self):
        # from Dickson et al. 2007
        return exp(
            -8904.2 / self.T
            + 117.385
            - 19.334 * log(self.T)
            + (-458.79 / self.T + 3.5913) * self.I**0.5
            + (188.74 / self.T - 1.5998) * self.I
            + (-12.1652 / self.T + 0.07871) * self.I**2
            + log(1 - 0.001005 * self.S)
        )

    @property
    def KP1(self):
        # from Dickson et al. 2007
        return exp(
            -4576.752 / self.T
            + 115.525
            - 18.453 * log(self.T)
            + (-106.736 / self.T + 0.69171) * self.S**0.5
            + (-0.65643 / self.T - 0.01844) * self.S
        )

    @property
    def KP2(self):
        # from Dickson et al. 2007
        return exp(
            -8814.715 / self.T
            + 172.0883
            - 27.927 * log(self.T)
            + (-160.340 / self.T + 1.3566) * self.S**0.5
            + (0.37335 / self.T - 0.05778) * self.S
        )

    @property
    def KP3(self):
        # from Dickson et al. 2007
        return exp(
            -3070.75 / self.T
            - 18.141
            + (17.27039 / self.T + 2.81197) * self.S**0.5
            + (-44.99486 / self.T - 0.09984) * self.S
        )


class titrant:

    def __init__(self, name: str, concentration: float):
        self.name = name
        self.c = concentration
        self.weight = list()

    @property
    def density(self):
        # this should read from a datasheet the first time, and save the density equation in memory

        a = 0
        b = 0
        c = 1
        d = 0

        return lambda weight: weight * a**3 + weight * b**2 + weight * c + d


class Titration:
    def __init__(
        self, weight: list[float], emf: list[float], temp: Union[float, list[float]]
    ):
        self.weight = np.array(weight)
        self.emf = np.array(emf)
        self.temp = np.array(temp)
